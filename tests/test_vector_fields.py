import numpy as np
from scipy import sparse

import anndata as ad

from cellucid.vector_fields import compute_transition_drift, add_transition_drift_to_obsm


def test_compute_transition_drift_identity_dense():
    rng = np.random.default_rng(0)
    X = rng.standard_normal((5, 2), dtype=np.float32)
    T = np.eye(5, dtype=np.float32)
    drift = compute_transition_drift(T, X)
    assert drift.dtype == np.float32
    assert drift.shape == X.shape
    assert np.allclose(drift, 0)


def test_compute_transition_drift_row_normalization_dense():
    rng = np.random.default_rng(1)
    X = rng.standard_normal((5, 3), dtype=np.float32)
    T = 2.0 * np.eye(5, dtype=np.float32)  # not row-stochastic, but proportional to identity
    drift = compute_transition_drift(T, X, normalize_rows=True)
    assert np.allclose(drift, 0)


def test_compute_transition_drift_sparse():
    rng = np.random.default_rng(2)
    X = rng.standard_normal((5, 2), dtype=np.float32)
    T = sparse.eye(5, format="csr", dtype=np.float32) * 3.0
    drift = compute_transition_drift(T, X, normalize_rows=True)
    assert np.allclose(drift, 0, atol=1e-6)


def test_add_transition_drift_to_obsm_infers_dim_from_explicit_key():
    adata = ad.AnnData(X=np.zeros((3, 0), dtype=np.float32))
    adata.obsm["X_umap_3d"] = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float32
    )
    T = np.eye(3, dtype=np.float32)
    key = add_transition_drift_to_obsm(adata, T, basis="umap", field_prefix="T_fwd")
    assert key == "T_fwd_umap_3d"
    assert key in adata.obsm
    assert adata.obsm[key].shape == (3, 3)
    assert np.allclose(adata.obsm[key], 0)


def test_add_transition_drift_to_obsm_infers_dim_from_base_key():
    adata = ad.AnnData(X=np.zeros((3, 0), dtype=np.float32))
    adata.obsm["X_umap"] = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    T = np.eye(3, dtype=np.float32)
    key = add_transition_drift_to_obsm(adata, T, basis="umap", field_prefix="T_fwd")
    assert key == "T_fwd_umap_2d"
    assert key in adata.obsm
    assert adata.obsm[key].shape == (3, 2)
    assert np.allclose(adata.obsm[key], 0)


def test_add_transition_drift_to_obsm_prefers_explicit_over_base_key():
    adata = ad.AnnData(X=np.zeros((3, 0), dtype=np.float32))
    adata.obsm["X_umap_2d"] = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    adata.obsm["X_umap"] = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float32
    )
    T = np.eye(3, dtype=np.float32)
    key = add_transition_drift_to_obsm(adata, T, basis="umap", field_prefix="T_fwd")
    assert key == "T_fwd_umap_2d"
    assert adata.obsm[key].shape == (3, 2)
