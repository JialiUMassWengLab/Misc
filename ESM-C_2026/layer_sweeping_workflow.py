#! /usr/bin/env python

import re
import os
import numpy as np
import pandas as pd
import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler


def get_hidden_states(trainingFile):
    model = AutoModelForMaskedLM.from_pretrained(
        "biohub/ESMC-600M",
        dtype="auto",
    ).eval()
    tokenizer = AutoTokenizer.from_pretrained("biohub/ESMC-600M")

    df = pd.read_csv(trainingFile, index_col=0)
    input_seqs = []
    input_ids = []
    for accession, seq in df['sequence'].to_dict().items():
        input_seqs.append(seq)
        input_ids.append(accession)

    embedding_tensors = {}
    batchSize = 10
    for i in range(0, len(input_seqs), batchSize):
        print('Processing proteins %d to %d' % (i+1,i+batchSize))
        inputs = tokenizer(input_seqs[i:i+batchSize], return_tensors="pt", padding=True)
        inputs = {k: v.to(model.device) for k, v in inputs.items()}
        with torch.inference_mode():
            output = model(**inputs, output_hidden_states=True)
            mask_f = inputs['attention_mask'][None, :, :, None].to(dtype=output.hidden_states.dtype)
            summed = (output.hidden_states * mask_f).sum(dim=2)
            counts = mask_f.sum(dim=2).clamp(min=1)
            mean_pooled = summed / counts
            embedding_tensors.update({'batch%d' % (int(i/batchSize) + 1): mean_pooled})

    torch.save(embedding_tensors, 'secreted_training_embedding_tensors.pt')


def layer_sweep_cv(
    embeddings_file='secreted_training_files/secreted_training_embedding_tensors.pt',
    labels_file='secreted_training_files/secreted_training_set.csv',
    n_splits=5,
    random_state=42,
):
    """5-fold CV logistic regression AUROC sweep across all embedding layers."""
    embeddings = torch.load(embeddings_file, map_location='cpu')
    df = pd.read_csv(labels_file, index_col=0)
    labels = df['label'].values

    # embeddings shape per batch: [n_layers, batch_size, hidden_dim]
    n_layers = next(iter(embeddings.values())).shape[0]

    rng = np.random.default_rng(random_state)
    shuffle_idx = rng.permutation(len(labels))
    labels_shuffled = labels[shuffle_idx]

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    layer_aurocs = []
    for layer in range(n_layers):
        # Concatenate all batches for this layer -> [n_proteins, hidden_dim]
        X = torch.cat([v[layer, :, :] for v in embeddings.values()], dim=0).numpy()
        X = X[shuffle_idx]

        fold_aurocs = []
        for train_idx, val_idx in cv.split(X, labels_shuffled):
            X_train, X_val = X[train_idx], X[val_idx]
            y_train, y_val = labels_shuffled[train_idx], labels_shuffled[val_idx]

            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_val = scaler.transform(X_val)

            clf = LogisticRegression(max_iter=1000, C=1e-7, random_state=random_state)
            clf.fit(X_train, y_train)
            proba = clf.predict_proba(X_val)[:, 1]
            fold_aurocs.append(roc_auc_score(y_val, proba))

        mean_auroc = float(np.mean(fold_aurocs))
        std_auroc = float(np.std(fold_aurocs))
        layer_aurocs.append({'layer': layer, 'mean_auroc': mean_auroc, 'std_auroc': std_auroc})
        print(f'Layer {layer:>2d}: AUROC = {mean_auroc:.4f} ± {std_auroc:.4f}')

    results = pd.DataFrame(layer_aurocs)
    best = results.loc[results['mean_auroc'].idxmax()]
    print(f'\nBest layer: {int(best["layer"])} (AUROC = {best["mean_auroc"]:.4f})')
    return results


def plot_layer_sweep(
    layers,
    auc_mean,
    auc_std,
    best_auc_layer,
    last_layer,
    save_path="ec_layer_sweep.png",
):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(layers, auc_mean, "o-", linewidth=2, label="MCC", color="steelblue")
    ax.fill_between(
        layers,
        auc_mean - auc_std,
        auc_mean + auc_std,
        alpha=0.2,
        color="steelblue",
        label="±1 std",
    )
    ax.axvline(
        best_auc_layer + 1,
        color="steelblue",
        linestyle=":",
        linewidth=1.5,
        label=f"Best MCC layer = {best_auc_layer + 1}",
    )
    ax.axvline(
        last_layer + 1,
        color="gray",
        linestyle=":",
        linewidth=1.5,
        label=f"Last layer = {last_layer + 1}",
    )
    ax.set_xlabel("Transformer layer", fontsize=12)
    ax.set_ylabel("CV performance", fontsize=12)
    ax.set_title(
        "Secreted protein classification performance by ESMC layer\n"
        "(mean-pooled embeddings + logistic regression, 5-fold CV)",
        fontsize=13,
    )
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.set_xlim(0.5, last_layer + 1.5)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=10)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to {save_path}")
    plt.show()


def fit_and_save_model(
    layer,
    embeddings_file='secreted_training_files/secreted_training_embedding_tensors.pt',
    labels_file='secreted_training_files/secreted_training_set.csv',
    output_file='secreted_training_files/logistic_regression_model.npz',
    random_state=42,
):
    """Fit logistic regression on the full dataset for a given layer and save weights."""
    embeddings = torch.load(embeddings_file, map_location='cpu')
    df = pd.read_csv(labels_file, index_col=0)
    labels = df['label'].values

    X = torch.cat([v[layer, :, :] for v in embeddings.values()], dim=0).numpy()

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    clf = LogisticRegression(C=0.01, max_iter=1000, random_state=random_state)
    clf.fit(X_scaled, labels)

    np.savez(
        output_file,
        coef=clf.coef_,
        intercept=clf.intercept_,
        scaler_mean=scaler.mean_,
        scaler_scale=scaler.scale_,
        layer=np.array([layer]),
        classes=clf.classes_,
    )
    print(f'Model for layer {layer} saved to {output_file}')


def fit_and_predict_on_test(
    layer,
    target_dir,
    random_state=42,
):
    training_embeddings_file = os.path.join(target_dir, 'secreted_training_embedding_tensors.pt')
    testing_embeddings_file = os.path.join(target_dir, 'secreted_testing_embedding_tensors.pt')
    labels_file = os.path.join(target_dir, 'secreted_training_set.csv')
    output_file = os.path.join(target_dir, 'secreted_testing_predict_prob.csv')

    """Fit logistic regression on the full dataset for a given layer and save weights."""
    embeddings1 = torch.load(training_embeddings_file, map_location='cpu')
    embeddings2 = torch.load(testing_embeddings_file, map_location='cpu')
    df = pd.read_csv(labels_file, index_col=0)
    labels = df['label'].values

    X1 = torch.cat([v[layer, :, :] for v in embeddings1.values()], dim=0).numpy()
    X2 = torch.cat([v[layer, :, :] for v in embeddings2.values()], dim=0).numpy()

    scaler = StandardScaler()
    X1_scaled = scaler.fit_transform(X1)
    X2_scaled = scaler.transform(X2)

    clf = LogisticRegression(C=1e-5, max_iter=1000, random_state=random_state)
    clf.fit(X1_scaled, labels)
    out = pd.Series(clf.predict_proba(X2_scaled)[:, 1], index=embeddings2.keys(), name='prob')

    out.to_csv(output_file)
    print(f'Prediction for layer {layer} saved to {output_file}')


def main():
    #get_hidden_states('secreted_training_set.csv')
    '''
    res = layer_sweep_cv()
    res.to_csv('secreted_training_auroc.csv', index=False)
    res = pd.read_csv('secreted_training_auroc.csv', index_col=0)
    N_LAYERS = res.shape[0]
    plot_layer_sweep(
        np.arange(1, N_LAYERS + 1),
        res['mean_auroc'],
        res['std_auroc'],
        int(res['mean_auroc'].argmax()),
        N_LAYERS - 1,
        save_path='secreted_training_auroc.png',
    )
    '''
    #fit_and_save_model(17)
    fit_and_predict_on_test(17, 'secreted_training_files')

if __name__ == '__main__':
    main()