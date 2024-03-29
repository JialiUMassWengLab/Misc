import re
import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import numpy as np

class proteomicsDataset(Dataset):
    def __init__(self,exprMat,labelDim):
        self.exprMat = exprMat.iloc[:,:-labelDim].astype(float)
        self.labels = exprMat.iloc[:,-labelDim:].astype(float)
        
    def __len__(self):
        return(self.labels.shape[0])
    
    def __getitem__(self,idx):
        #exprVector = np.array(self.exprMat.iloc[idx,:])
        #label = np.array(self.labels.iloc[idx,:])
        exprVector = torch.tensor(self.exprMat.iloc[idx,:], dtype=torch.float64)
        label = torch.tensor(self.labels.iloc[idx,:], dtype=torch.float64)
        return exprVector,label

    
class VariationalAutoEncoder(nn.Module):
    def __init__(self, input_dim, z_dim, h_dim=200):
        super().__init__()
        # encoder
        self.img_2hid = nn.Linear(input_dim, h_dim)

        # one for mu and one for stds, note how we only output
        # diagonal values of covariance matrix. Here we assume
        # the pixels are conditionally independent 
        self.hid_2mu = nn.Linear(h_dim, z_dim)
        self.hid_2sigma = nn.Linear(h_dim, z_dim)

        # decoder
        self.z_2hid = nn.Linear(z_dim, h_dim)
        self.hid_2img = nn.Linear(h_dim, input_dim)
        
        # decoder2 (for computing extra terms involving sigma)
        self.z_2hid2 = nn.Linear(z_dim, h_dim, bias=False)
        self.hid_2img2 = nn.Linear(h_dim, input_dim, bias=False)
        del self.z_2hid2.weight
        del self.hid_2img2.weight
        
        self.double()

    def encode(self, x):
        h = F.relu(self.img_2hid(x))
        mu = self.hid_2mu(h)
        sigma = self.hid_2sigma(h)
        return mu, sigma

    def decode(self, z):
        new_h = self.z_2hid(z)
        x = self.hid_2img(new_h)
        return x
    
    def extra(self, sigma):
        CA = torch.matmul(self.hid_2img.weight,self.z_2hid.weight)
        cov = torch.diag_embed(sigma)
        add_on = torch.diagonal(torch.matmul(torch.matmul(CA,cov),CA.transpose(0,1)),dim1=-2,dim2=-1)
        return add_on

    def forward(self, x):
        mu, sigma = self.encode(x)
        sigma = torch.exp(sigma)

        x_reconst = self.decode(mu)
        add_on = self.extra(sigma)
        
        return x_reconst, add_on, mu, sigma


def readProteomicsNPX():
    npx = pd.read_csv('combined_pheno_forconsortium_v1_NPX.tsv',sep='\t',index_col=0,low_memory=False)
    missingRatioProt = npx.apply(lambda x: x.isna().sum()/npx.shape[0],axis=0)
    npx = npx.loc[:,list(missingRatioProt[missingRatioProt < .1].index)]
    missingRatioSamp = npx.apply(lambda x: x.isna().sum()/npx.shape[1],axis=1)
    npx = npx.loc[list(missingRatioSamp[missingRatioSamp < .2].index),:]

    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info = info.loc[:,map(lambda x: re.search(r'in urine',x)==None,info.columns)]
    info['Sex'] = info.apply(lambda x: 0 if x['Sex']=='F' else 1 if x['Sex']=='M' else np.nan,axis=1)
    print(info.shape)
    df = pd.concat([npx,info],join='inner',axis=1)
    return df,info.shape[1]


def train(train_loader,device,INPUT_DIM,num_epochs,model,optimizer,loss_fn):
    # Start training
    for epoch in range(num_epochs):
        loop = tqdm(enumerate(train_loader))
        for i, (x, y) in loop:
            # Forward pass
            x = x.to(device).view(-1, INPUT_DIM)
            nan_in_x = torch.isnan(x)
            x = torch.nan_to_num(x)
            x_reconst, add_on, mu, sigma = model(x)
            
            # loss, formulas from https://www.youtube.com/watch?v=igP03FXZqgo&t=2182s
            reconst_loss = loss_fn(x_reconst[~nan_in_x], x[~nan_in_x]) + add_on[~nan_in_x].sum()
            kl_div = - torch.sum(1 + torch.log(sigma.pow(2)) - mu.pow(2) - sigma.pow(2))

            # Backprop and optimize
            loss = reconst_loss + kl_div
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            loop.set_postfix(loss=loss.item())


def main():
    df,labelDim = readProteomicsNPX()
    dataset = proteomicsDataset(df,labelDim)
    print(len(dataset))

    # Configuration
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    INPUT_DIM = 2931
    Z_DIM = 20
    H_DIM = 250
    NUM_EPOCHS = 20
    BATCH_SIZE = 80
    LR_RATE = 3e-4

    # Initialize model, optimizer, loss
    model = VariationalAutoEncoder(INPUT_DIM, Z_DIM, H_DIM).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR_RATE)
    loss_fn = nn.MSELoss(reduction="sum")

    # Run training
    train_loader = DataLoader(dataset=dataset, batch_size=BATCH_SIZE, shuffle=True)
    train(train_loader,device,INPUT_DIM,NUM_EPOCHS,model,optimizer,loss_fn)

    # Get Mu's and Sigma's and reconstructed x's for each sample
    mus = []
    sigmas = []
    reconst_x = []
    for vector,label in dataset:
        vector = torch.nan_to_num(vector)
        with torch.no_grad():
            mu, sigma = model.encode(vector)
            sigma = torch.exp(sigma)
            mus.append(mu)
            sigmas.append(sigma)
        
            # Sample from latent distribution from encoder
            epsilon = torch.randn_like(sigma)
            z_reparametrized = mu + sigma*epsilon
            x = model.decode(z_reparametrized)
            reconst_x.append(x)

    embed_array = []
    reconst_array = []
    for i,(vector,label) in enumerate(dataset):
        embed_array.append(pd.Series(np.array(mus[i]),name=df.index[i]))
        reconst_array.append(pd.Series(np.array(reconst_x[i]),name=df.index[i]))

    embedding = pd.concat(embed_array,join='inner',axis=1).transpose()
    embedding.columns = ['latent'+str(i) for i in embedding.columns]
    recon = pd.concat(reconst_array,join='inner',axis=1).transpose()
    recon.columns = df.columns[:-labelDim]

    embedding.to_csv('VAE_noSampling_embeddings.csv')
    recon.to_csv('VAE_noSampling_reconst.csv')


if __name__=='__main__':
    main()
