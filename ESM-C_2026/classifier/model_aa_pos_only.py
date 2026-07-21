"""PositionOnlyClassifier: MLP over the substitution-position embedding vector.

Baseline model that uses only the (1152,) embedding-difference vector at the
amino-acid substitution position, skipping attention pooling entirely.
Motivated by the observation that ||diff[N]|| is 26-107x larger than the mean
||diff[other]|| in this dataset, so most of the LOF-vs-GOF signal is at N.
"""

from __future__ import annotations

import torch
import torch.nn as nn


class PositionOnlyClassifier(nn.Module):
    """MLP: (B, 1152) -> (B, 1) logits. No attention, no pooling, no mask.

    Same MLP head shape as AttentionPoolClassifier (Linear-GELU-Dropout stack)
    so parameter counts differ only by the attention block. Sigmoid is NOT
    applied here — training uses BCEWithLogitsLoss.
    """

    def __init__(
        self,
        embed_dim: int = 1152,
        hidden_dims: tuple[int, ...] = (512, 128),
        dropout: float = 0.2,
    ) -> None:
        super().__init__()
        layers: list[nn.Module] = []
        prev = embed_dim
        for h in hidden_dims:
            layers += [nn.Linear(prev, h), nn.GELU(), nn.Dropout(dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.mlp = nn.Sequential(*layers)

    def forward(
        self,
        embedding: torch.Tensor,
        attention_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        # attention_mask is accepted for signature-compatibility with
        # AttentionPoolClassifier so train.run_epoch can call both uniformly.
        del attention_mask
        return self.mlp(embedding)

    @torch.no_grad()
    def predict_proba(
        self,
        embedding: torch.Tensor,
        attention_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        return torch.sigmoid(self.forward(embedding, attention_mask))
