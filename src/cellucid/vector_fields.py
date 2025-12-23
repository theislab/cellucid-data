"""
Vector field utilities for Cellucid.

Cellucid visualizes "vector fields" (animated arrow/flow overlays) as **per-cell
displacement vectors in embedding space**. This module provides helpers to
derive such vectors from common CellRank outputs (transition matrices) and to
store them in AnnData `obsm` using Cellucid’s naming convention.
"""

from __future__ import annotations

from typing import Optional, Union

import numpy as np
from scipy import sparse


def compute_transition_drift(
    transition_matrix: Union[np.ndarray, sparse.spmatrix],
    embedding: np.ndarray,
    *,
    normalize_rows: bool = True,
) -> np.ndarray:
    """
    Compute a per-cell drift vector field from a transition matrix.

    The drift is defined as:
        drift = E[next_embedding | current_cell] - current_embedding

    Where:
        E[next_embedding] = T @ embedding

    Parameters
    ----------
    transition_matrix:
        A (n_cells, n_cells) matrix. Sparse matrices are recommended for large datasets.
        Typical source: CellRank kernels' transition matrices.
    embedding:
        A (n_cells, dim) array of embedding coordinates (e.g. `adata.obsm['X_umap']`).
    normalize_rows:
        If True (default), divides each row's expectation by its row-sum. This makes
        the drift invariant to a matrix that isn't strictly row-stochastic.

    Returns
    -------
    np.ndarray
        Drift vectors shaped (n_cells, dim) as float32.
    """
    X = np.asarray(embedding, dtype=np.float32)
    if X.ndim != 2:
        raise ValueError(f"embedding must be 2D, got shape {X.shape}")

    if sparse.issparse(transition_matrix):
        T = transition_matrix.tocsr()
        TX = T.dot(X)
        if normalize_rows:
            row_sums = np.asarray(T.sum(axis=1)).reshape(-1)
            safe = row_sums.copy()
            safe[safe == 0] = 1.0
            TX = TX / safe[:, None]
    else:
        T = np.asarray(transition_matrix, dtype=np.float32)
        if T.ndim != 2:
            raise ValueError(f"transition_matrix must be 2D, got shape {T.shape}")
        TX = T @ X
        if normalize_rows:
            row_sums = T.sum(axis=1)
            safe = np.where(row_sums == 0, 1.0, row_sums).astype(np.float32)
            TX = TX / safe[:, None]

    if TX.shape != X.shape:
        raise ValueError(f"T @ embedding produced shape {TX.shape}, expected {X.shape}")

    return (TX - X).astype(np.float32, copy=False)


def _resolve_embedding_key(
    adata: "anndata.AnnData",
    *,
    basis: str,
    dim: int,
) -> str:
    """
    Resolve an AnnData `obsm` key for an embedding basis + dimension.

    Prefers explicit Cellucid-style keys (e.g. `X_umap_2d`) and falls back to
    Scanpy-style `X_umap` if it matches the requested dimension.
    """
    key_explicit = f"X_{basis}_{dim}d"
    if key_explicit in adata.obsm:
        return key_explicit

    key_base = f"X_{basis}"
    if key_base in adata.obsm:
        arr = adata.obsm[key_base]
        shape = arr.shape
        if len(shape) == 2 and shape[1] == dim:
            return key_base
        raise ValueError(f"{key_base} is {shape[1]}D, expected {dim}D")

    raise KeyError(f"Missing embedding in obsm: {key_explicit} (or {key_base})")


def add_transition_drift_to_obsm(
    adata: "anndata.AnnData",
    transition_matrix: Union[np.ndarray, sparse.spmatrix],
    *,
    basis: str = "umap",
    field_prefix: str = "T_fwd",
    dim: Optional[int] = None,
    explicit_dim_suffix: bool = True,
    normalize_rows: bool = True,
    overwrite: bool = False,
) -> str:
    """
    Compute drift vectors from a transition matrix and store them in `adata.obsm`.

    The output key follows Cellucid’s vector field naming convention:
      - Explicit: `<field_prefix>_<basis>_<dim>d`  (default)
      - Implicit: `<field_prefix>_<basis>`

    Examples
    --------
    - Forward drift in UMAP (2D): `T_fwd_umap_2d`
    - Backward drift in UMAP (2D): `T_bwd_umap_2d`

    Parameters
    ----------
    adata:
        AnnData object to modify.
    transition_matrix:
        (n_cells, n_cells) matrix, dense or sparse.
    basis:
        Embedding basis name (default: "umap").
    field_prefix:
        Prefix for the vector field (default: "T_fwd").
        Use e.g. "T_bwd" for backward matrices.
    dim:
        Embedding dimensionality to use. If None, it is inferred from the best
        available embedding in `obsm`:
        `X_{basis}_3d` → `X_{basis}_2d` → `X_{basis}_1d` → `X_{basis}` (if 1D/2D/3D).
    explicit_dim_suffix:
        If True, appends `_{dim}d` to the output key (recommended for clash-safe imports).
    normalize_rows:
        Whether to row-normalize during drift computation (see `compute_transition_drift`).
    overwrite:
        If False (default), raises if the target key already exists.

    Returns
    -------
    str
        The `obsm` key written.
    """
    if dim is None:
        # Prefer explicit Cellucid-style keys (X_<basis>_<dim>d) and fall back to
        # Scanpy-style X_<basis> only if it matches 1D/2D/3D.
        for candidate_dim in (3, 2, 1):
            explicit_key = f"X_{basis}_{candidate_dim}d"
            if explicit_key not in adata.obsm:
                continue
            shape = adata.obsm[explicit_key].shape
            if len(shape) != 2 or shape[1] != candidate_dim:
                raise ValueError(f"{explicit_key} has shape {shape}, expected (n_cells, {candidate_dim})")
            dim = candidate_dim
            break

        if dim is None:
            base_key = f"X_{basis}"
            if base_key not in adata.obsm:
                raise KeyError(
                    f"dim not provided and no embedding found: "
                    f"{base_key} or X_{basis}_{{1,2,3}}d"
                )
            shape = adata.obsm[base_key].shape
            if len(shape) == 1:
                dim = 1
            elif len(shape) == 2 and shape[1] in (1, 2, 3):
                dim = int(shape[1])
            else:
                raise ValueError(f"Cannot infer dim from {base_key} with shape {shape}")

    assert dim is not None
    emb_key = _resolve_embedding_key(adata, basis=basis, dim=dim)
    X = np.asarray(adata.obsm[emb_key], dtype=np.float32)
    drift = compute_transition_drift(transition_matrix, X, normalize_rows=normalize_rows)

    base_id = f"{field_prefix}_{basis}"
    out_key = f"{base_id}_{dim}d" if explicit_dim_suffix else base_id

    if not overwrite and out_key in adata.obsm:
        raise KeyError(f"adata.obsm already contains key '{out_key}' (set overwrite=True to replace)")

    adata.obsm[out_key] = drift
    return out_key
