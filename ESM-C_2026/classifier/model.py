"""AttentionPoolClassifier: learnable-query attention pooling + MLP head.

Two pooling variants selectable via `pool_type`:
- "projected": scores = (Linear(x) @ query) / sqrt(D), the original variant.
- "direct":    scores = x @ w, where w is a learnable (D,) weight vector, no
               projection and no scale factor. Follows spec:
               scores = w @ E,  alpha = softmax(scores),  pooled = E @ alpha^T
               (E has shape (D, seq_len); we operate on the transposed layout).
"""

from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class AttentionPoolClassifier(nn.Module):
    """Pools (B, seq_len, embed_dim) -> (B, embed_dim) via a single learned query,
    then runs an MLP to a single logit.

    Sigmoid is NOT applied in forward — training uses BCEWithLogitsLoss which
    fuses sigmoid + BCE for numerical stability. Use predict_proba() at inference.
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
        if self.pool_type == "projected":
            keys = self.attn_proj(x)
            scores = torch.einsum("bsd,d->bs", keys, self.query) / math.sqrt(self.embed_dim)
        else:
            scores = torch.einsum("bsd,d->bs", x, self.query)
        scores = scores.masked_fill(~attention_mask, float("-inf"))
        weights = F.softmax(scores, dim=-1)
        return torch.einsum("bs,bsd->bd", weights, x)

    def forward(
        self, embedding: torch.Tensor, attention_mask: torch.Tensor
    ) -> torch.Tensor:
        pooled = self._attention_pool(embedding, attention_mask)
        return self.mlp(pooled)

    @torch.no_grad()
    def predict_proba(
        self, embedding: torch.Tensor, attention_mask: torch.Tensor
    ) -> torch.Tensor:
        return torch.sigmoid(self.forward(embedding, attention_mask))
