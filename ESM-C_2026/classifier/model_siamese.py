"""SiameseAttentionClassifier: shared-weight attention pooling on WT and variant
sequences, then MLP on the difference of the pooled vectors.

Contrast with AttentionPoolClassifier (which pools the pre-computed
WT-minus-variant difference tensor):
- diff-input model:      MLP( attention_pool( WT_emb - variant_emb ) )
- this siamese model:    MLP( attention_pool(WT_emb) - attention_pool(variant_emb) )

The attention block is applied twice with the SAME weights, so the two pooled
vectors live in the same latent space and their subtraction is meaningful.

Sigmoid is NOT applied in forward — training uses BCEWithLogitsLoss.
"""

from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class SiameseAttentionClassifier(nn.Module):
    """Two-input model: shared attention encoder + MLP on pooled difference.

    forward(wt_embedding, wt_mask, var_embedding, var_mask) -> logits (B, 1)
    """

    def __init__(
        self,
        embed_dim: int = 1152,
        hidden_dims: tuple[int, ...] = (512, 128),
        dropout: float = 0.2,
        pool_type: str = "projected",
    ) -> None:
        super().__init__()
        if pool_type not in {"projected", "direct"}:
            raise ValueError(f"pool_type must be 'projected' or 'direct', got {pool_type!r}")
        self.embed_dim = embed_dim
        self.pool_type = pool_type

        # Shared attention block: one learnable query, one projection matrix
        # applied to BOTH WT and variant embeddings.
        self.query = nn.Parameter(torch.randn(embed_dim) * 0.02)
        if pool_type == "projected":
            self.attn_proj = nn.Linear(embed_dim, embed_dim)

        mlp_layers: list[nn.Module] = []
        prev = embed_dim
        for h in hidden_dims:
            mlp_layers += [nn.Linear(prev, h), nn.GELU(), nn.Dropout(dropout)]
            prev = h
        mlp_layers.append(nn.Linear(prev, 1))
        self.mlp = nn.Sequential(*mlp_layers)

    def _attention_pool(
        self, x: torch.Tensor, attention_mask: torch.Tensor
    ) -> torch.Tensor:
        # x: (B, seq_len, D)   mask: (B, seq_len) bool, True at valid positions
        if self.pool_type == "projected":
            keys = self.attn_proj(x)
            scores = torch.einsum("bsd,d->bs", keys, self.query) / math.sqrt(self.embed_dim)
        else:
            scores = torch.einsum("bsd,d->bs", x, self.query)
        scores = scores.masked_fill(~attention_mask, float("-inf"))
        weights = F.softmax(scores, dim=-1)
        return torch.einsum("bs,bsd->bd", weights, x)

    def forward(
        self,
        wt_embedding: torch.Tensor,
        wt_mask: torch.Tensor,
        var_embedding: torch.Tensor,
        var_mask: torch.Tensor,
    ) -> torch.Tensor:
        # Attention pool each sequence independently through the SAME weights.
        wt_pooled = self._attention_pool(wt_embedding, wt_mask)
        var_pooled = self._attention_pool(var_embedding, var_mask)
        # Subtract in latent space, then classify.
        return self.mlp(wt_pooled - var_pooled)

    @torch.no_grad()
    def predict_proba(
        self,
        wt_embedding: torch.Tensor,
        wt_mask: torch.Tensor,
        var_embedding: torch.Tensor,
        var_mask: torch.Tensor,
    ) -> torch.Tensor:
        return torch.sigmoid(
            self.forward(wt_embedding, wt_mask, var_embedding, var_mask)
        )
