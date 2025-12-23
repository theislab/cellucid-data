#!/usr/bin/env python3
"""
Export raw dataframes/arrays to files used by the WebGL viewer.

Includes memory/disk optimization features:
- Quantization for continuous data (var/gene expression and obs continuous)
- Auto dtype selection for categorical obs based on category count
- Gzip compression for all binary files

Instead of AnnData, accepts:
- latent_space: (n_cells, n_dims) numpy/sparse array for outlier quantile calculation
- X_umap_1d / X_umap_2d / X_umap_3d: explicit embeddings (at least one required)
- vector_fields: dict[str, array] of per-cell displacement vectors (optional)
- obs: pandas DataFrame with cell metadata columns
- var: pandas DataFrame with gene/feature metadata
- gene_expression: (n_cells, n_genes) numpy/sparse array for gene expression matrix
- var_gene_id_column: column name in var containing gene identifiers, or "index" to use var.index
- connectivities: sparse matrix with KNN connectivities
"""
import gzip
import json
import re
import tqdm
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence, Union, Literal

import numpy as np
import pandas as pd
from scipy import sparse

DEFAULT_EXPORT_DIR = Path.cwd() / "exports"
DEFAULT_OBS_DIRNAME = "obs"
DEFAULT_VAR_DIRNAME = "var"
DEFAULT_CONNECTIVITY_DIRNAME = "connectivity"

# Manifest format version for compact format
MANIFEST_FORMAT_VERSION = "compact_v1"


def _safe_filename_component(name: str) -> str:
    """Return a filesystem-friendly version of a field key."""
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name))
    safe = safe.strip("._")
    return safe or "field"


def _to_dense(arr: Union[np.ndarray, sparse.spmatrix]) -> np.ndarray:
    """Convert sparse matrix to dense numpy array if necessary."""
    if sparse.issparse(arr):
        return np.asarray(arr.toarray())
    return np.asarray(arr)


def _file_exists_skip(path: Path, description: str, force: bool = False) -> bool:
    """Check if file exists and should be skipped. Returns True if should skip."""
    if path.exists() and not force:
        print(f"⚠ Skipping {description}: {path} already exists (use force=True to overwrite)")
        return True
    return False


def _write_binary(
    path: Path,
    data: np.ndarray,
    compression: Optional[int] = None,
) -> Path:
    """
    Write binary data, optionally with gzip compression.
    
    Parameters
    ----------
    path : Path
        Output path. If compression is enabled, '.gz' will be appended.
    data : np.ndarray
        Data to write.
    compression : int or None
        Gzip compression level (1-9). None or 0 means no compression.
        
    Returns
    -------
    Path
        Actual path written (may have .gz suffix).
    """
    if compression and compression > 0:
        gz_path = Path(str(path) + ".gz")
        with gzip.open(gz_path, 'wb', compresslevel=compression) as f:
            f.write(data.tobytes())
        return gz_path
    else:
        data.tofile(path)
        return path


def _quantize_continuous(
    values: np.ndarray,
    bits: int = 8,
    field_name: str = "unknown",
) -> tuple[np.ndarray, float, float, float, dict]:
    """
    Quantize continuous float32 values to uint8 or uint16.
    
    Parameters
    ----------
    values : np.ndarray
        Float32 values to quantize.
    bits : int
        Number of bits for quantization (8 or 16).
    field_name : str
        Name of field for debug messages.
        
    Returns
    -------
    quantized : np.ndarray
        Quantized values as uint8 or uint16.
    min_val : float
        Minimum value (for dequantization).
    max_val : float
        Maximum value (for dequantization).
    scale : float
        Scale factor (for dequantization).
    stats : dict
        Statistics about the quantization for debugging.
    """
    n_total = len(values)
    
    # Identify problematic values
    nan_mask = np.isnan(values)
    inf_mask = np.isinf(values)
    invalid_mask = nan_mask | inf_mask
    valid_mask = ~invalid_mask
    
    n_nan = int(nan_mask.sum())
    n_inf = int(inf_mask.sum())
    n_valid = int(valid_mask.sum())
    
    stats = {
        "n_total": n_total,
        "n_valid": n_valid,
        "n_nan": n_nan,
        "n_inf": n_inf,
    }
    
    if n_valid == 0:
        # All invalid values
        min_val, max_val = 0.0, 1.0
        stats["warning"] = "all_invalid"
    else:
        valid_values = values[valid_mask]
        min_val = float(np.min(valid_values))
        max_val = float(np.max(valid_values))
        stats["data_min"] = min_val
        stats["data_max"] = max_val
    
    # Avoid division by zero
    if max_val == min_val:
        max_val = min_val + 1.0
        stats["constant_value"] = True
    
    if bits == 8:
        max_quant = 254  # Reserve 255 for NaN/Inf
        dtype = np.uint8
        nan_value = 255
    else:  # 16 bits
        max_quant = 65534  # Reserve 65535 for NaN/Inf
        dtype = np.uint16
        nan_value = 65535
    
    scale = max_quant / (max_val - min_val)
    
    # Create output array
    quantized = np.empty(n_total, dtype=dtype)
    
    # Only quantize valid values to avoid numpy warnings
    if n_valid > 0:
        normalized = (values[valid_mask] - min_val) * scale
        quantized[valid_mask] = np.clip(normalized, 0, max_quant).astype(dtype)
    
    # Set invalid values to reserved marker
    quantized[invalid_mask] = nan_value
    
    return quantized, min_val, max_val, scale, stats


def _select_category_dtype(n_categories: int) -> tuple[np.dtype, int]:
    """
    Select optimal dtype for category codes based on number of categories.
    
    Parameters
    ----------
    n_categories : int
        Number of unique categories (not counting missing).
        
    Returns
    -------
    dtype : np.dtype
        Optimal dtype (uint8 or uint16).
    missing_value : int
        Value to use for missing/NaN codes.
    """
    if n_categories <= 254:
        # uint8 can hold 0-254 for categories, 255 for missing
        return np.uint8, 255
    else:
        # uint16 can hold 0-65534 for categories, 65535 for missing
        return np.uint16, 65535


def _report_quantization_stats(field_name: str, stats: dict, field_type: str = "field") -> None:
    """
    Report quantization statistics and warnings.
    
    Parameters
    ----------
    field_name : str
        Name of the field being quantized.
    stats : dict
        Statistics from _quantize_continuous.
    field_type : str
        Type of field for messages (e.g., "obs field", "gene", "outlier quantiles").
    """
    n_nan = stats.get("n_nan", 0)
    n_inf = stats.get("n_inf", 0)
    n_total = stats.get("n_total", 0)
    n_valid = stats.get("n_valid", 0)
    
    issues = []
    
    if n_nan > 0:
        pct = 100 * n_nan / n_total if n_total > 0 else 0
        issues.append(f"{n_nan:,} NaN values ({pct:.1f}%)")
    
    if n_inf > 0:
        pct = 100 * n_inf / n_total if n_total > 0 else 0
        issues.append(f"{n_inf:,} Inf values ({pct:.1f}%)")
    
    if stats.get("warning") == "all_invalid":
        print(f"  ⚠ WARNING: {field_type} '{field_name}' has NO valid values (all NaN/Inf)")
    elif issues:
        print(f"  ⚠ {field_type} '{field_name}': {', '.join(issues)} → mapped to missing marker")
    
    if stats.get("constant_value"):
        print(f"  ℹ {field_type} '{field_name}' has constant value (min=max)")


def _check_data_quality(values: np.ndarray, field_name: str, field_type: str = "field") -> None:
    """
    Check data quality and print warnings for common issues.
    
    Parameters
    ----------
    values : np.ndarray
        Values to check.
    field_name : str
        Name of the field.
    field_type : str
        Type of field for messages.
    """
    n_total = len(values)
    if n_total == 0:
        print(f"  ⚠ WARNING: {field_type} '{field_name}' is empty")
        return
    
    n_nan = int(np.isnan(values).sum())
    n_inf = int(np.isinf(values).sum())
    n_neg_inf = int(np.isneginf(values).sum())
    n_pos_inf = int(np.isposinf(values).sum())
    
    issues = []
    if n_nan > 0:
        pct = 100 * n_nan / n_total
        issues.append(f"{n_nan:,} NaN ({pct:.1f}%)")
    if n_neg_inf > 0:
        issues.append(f"{n_neg_inf:,} -Inf")
    if n_pos_inf > 0:
        issues.append(f"{n_pos_inf:,} +Inf")
    
    if issues:
        print(f"  ⚠ {field_type} '{field_name}' contains: {', '.join(issues)}")


def _compute_centroids_for_field(
    coords: np.ndarray,
    codes: np.ndarray,
    categories: list[str],
    outlier_quantile: float = 0.95,
    min_points: int = 10,
) -> list[dict]:
    """
    Compute centroids per category with outlier removal (based on embedding coords for display).

    Works with any dimensionality (1D, 2D, 3D, etc.) - the position will have the same
    number of dimensions as the input coords.
    """
    if coords.shape[0] != codes.shape[0]:
        raise ValueError("coords and codes must have the same length.")

    centroids: list[dict] = []

    if not (0.5 < outlier_quantile < 1.0):
        outlier_quantile = 0.95

    for code, label in enumerate(categories):
        mask = codes == code
        idx = np.nonzero(mask)[0]
        n = idx.size
        if n < min_points:
            continue

        pts = coords[idx, :]  # (n, ndim)
        center = pts.mean(axis=0)

        if n > min_points:
            dists = np.linalg.norm(pts - center, axis=1)
            thr = float(np.quantile(dists, outlier_quantile))
            inlier_mask = dists <= thr
            n_in = int(inlier_mask.sum())
            if n_in >= min_points:
                pts_in = pts[inlier_mask, :]
                center = pts_in.mean(axis=0)
                used_count = n_in
            else:
                used_count = n
        else:
            used_count = n

        centroids.append(
            {
                "category": str(label),
                "position": center.astype(float).tolist(),
                "n_points": int(used_count),
            }
        )

    return centroids


def _compute_centroids_for_all_dimensions(
    embeddings: dict[int, np.ndarray],
    codes: np.ndarray,
    categories: list[str],
    outlier_quantile: float = 0.95,
    min_points: int = 10,
) -> dict[int, list[dict]]:
    """
    Compute centroids for each available dimension.

    Returns a dictionary keyed by dimension (1, 2, 3) with centroid lists for each.
    """
    centroids_by_dim = {}
    for dim, coords in embeddings.items():
        centroids_by_dim[dim] = _compute_centroids_for_field(
            coords, codes, categories, outlier_quantile, min_points
        )
    return centroids_by_dim


def _compute_latent_space_quantiles(
    latent: np.ndarray,
    codes: np.ndarray,
    categories: list[str],
    min_points: int = 10,
) -> np.ndarray:
    """
    Compute per-cell outlier quantiles based on latent space distances to category centroids.
    """
    n_cells = latent.shape[0]
    quantiles = np.full(n_cells, np.nan, dtype=np.float32)
    
    for code, label in enumerate(categories):
        mask = codes == code
        idx = np.nonzero(mask)[0]
        n = idx.size
        
        if n < min_points:
            continue
        
        pts = latent[idx, :]
        centroid = pts.mean(axis=0)
        dists = np.linalg.norm(pts - centroid, axis=1)
        sorted_dists = np.sort(dists)
        ranks = np.searchsorted(sorted_dists, dists, side='right')
        cell_quantiles = ranks.astype(np.float32) / n
        quantiles[idx] = cell_quantiles
    
    return quantiles


def prepare(
    latent_space: Optional[Union[np.ndarray, sparse.spmatrix]] = None,
    obs: Optional[pd.DataFrame] = None,
    var: Optional[pd.DataFrame] = None,
    gene_expression: Optional[Union[np.ndarray, sparse.spmatrix]] = None,
    var_gene_id_column: str = "index",
    gene_identifiers: Optional[Sequence[str]] = None,
    connectivities: Optional[sparse.spmatrix] = None,
    out_dir: Path | str = DEFAULT_EXPORT_DIR,
    obs_keys: Optional[Sequence[str]] = None,
    centroid_outlier_quantile: float = 0.95,
    centroid_min_points: int = 10,
    obs_manifest_filename: str = "obs_manifest.json",
    obs_binary_dirname: str = DEFAULT_OBS_DIRNAME,
    var_manifest_filename: str = "var_manifest.json",
    var_binary_dirname: str = DEFAULT_VAR_DIRNAME,
    connectivity_manifest_filename: str = "connectivity_manifest.json",
    connectivity_binary_dirname: str = DEFAULT_CONNECTIVITY_DIRNAME,
    force: bool = False,
    # Optimization parameters (all disabled by default)
    var_quantization: Optional[int] = None,
    obs_continuous_quantization: Optional[int] = None,
    obs_categorical_dtype: Literal["auto", "uint8", "uint16"] = "auto",
    compression: Optional[int] = None,
    # Dataset metadata parameters (for dataset_identity.json)
    dataset_name: Optional[str] = None,
    dataset_description: Optional[str] = None,
    dataset_id: Optional[str] = None,
    source_name: Optional[str] = None,
    source_url: Optional[str] = None,
    source_citation: Optional[str] = None,
    # Multi-dimensional embedding parameters (at least one required)
    X_umap_1d: Optional[np.ndarray] = None,
    X_umap_2d: Optional[np.ndarray] = None,
    X_umap_3d: Optional[np.ndarray] = None,
    X_umap_4d: Optional[np.ndarray] = None,
    # Optional per-cell vector fields aligned to the embedding(s)
    # (e.g. scVelo velocity, CellRank drift vectors)
    vector_fields: Optional[dict[str, Union[np.ndarray, sparse.spmatrix]]] = None,
) -> None:
    """
    Export raw data arrays to files used by the WebGL viewer.

    Memory/Disk Optimization Options
    --------------------------------
    var_quantization : int or None
        Bits for gene expression quantization (8, 16, or None for full float32).
        8-bit reduces file size by 4x with minimal visual impact for colormapping.
    obs_continuous_quantization : int or None
        Bits for continuous obs field quantization (8, 16, or None for full float32).
    obs_categorical_dtype : 'auto', 'uint8', or 'uint16'
        - 'auto': Select based on number of categories (uint8 if ≤254, else uint16)
        - 'uint8': Force uint8 (max 254 categories)
        - 'uint16': Force uint16 (max 65534 categories)
    compression : int or None
        Gzip compression level (1-9). None or 0 disables compression.
        Level 6 is a good balance of speed and size. Files get .gz extension.

    Multi-Dimensional Embeddings
    ----------------------------
    At least one dimensional embedding must be provided. The viewer supports
    switching between different dimensionalities of the same data at runtime.
    All embeddings must have the same number of cells (rows) but different
    column counts matching their dimensionality.

    IMPORTANT: Each embedding is normalized independently to fit within the
    [-1, 1] coordinate range. Within each dimension, the same scale factor is
    used for all axes to preserve aspect ratios. This ensures each dimension
    fills the viewing area optimally without requiring manual zoom adjustment.

    X_umap_1d : np.ndarray, optional
        1D embedding coordinates, shape (n_cells, 1). Stored as points_1d.bin.
    X_umap_2d : np.ndarray, optional
        2D embedding coordinates, shape (n_cells, 2). Stored as points_2d.bin.
    X_umap_3d : np.ndarray, optional
        3D embedding coordinates, shape (n_cells, 3). Stored as points_3d.bin.
        This is the primary visualization and is used for centroid computation.
    X_umap_4d : np.ndarray, optional
        4D embedding coordinates, shape (n_cells, 4). Stored as points_4d.bin.
        NOTE: 4D visualization is not yet implemented in the viewer.

    vector_fields : dict[str, np.ndarray] or None
        Optional per-cell displacement vectors aligned to the embedding space.
        Keys follow the same naming convention as AnnData `obsm`:
        - Explicit: `<field>_umap_<dim>d` (e.g. `velocity_umap_2d`, `T_fwd_umap_3d`)
        - Implicit: `<field>_umap` with shape `(n_cells, 1|2|3)`
          (used only if the explicit key for that dim is not provided)

        Each value must be shaped `(n_cells, dim)` (or `(n_cells,)` for 1D).
        Vectors are scaled by the SAME per-dimension normalization scale as points.

    Standard Parameters
    -------------------
    latent_space : np.ndarray or sparse matrix
        Latent space for outlier quantile calculation, shape (n_cells, n_dims).
    obs : pd.DataFrame
        Cell metadata, shape (n_cells, n_obs_columns).
    var : pd.DataFrame, optional
        Gene/feature metadata. Required if gene_expression is provided.
    gene_expression : np.ndarray or sparse matrix, optional
        Gene expression matrix, shape (n_cells, n_genes).
    var_gene_id_column : str
        Column name in var containing gene identifiers, or "index" to use var.index.
    gene_identifiers : sequence of str, optional
        Which genes to export. If None, all genes are exported.
    connectivities : sparse matrix, optional
        KNN connectivity matrix from scanpy (n_cells, n_cells).
    out_dir : Path or str
        Output directory (default: exports/ under the current working directory).
    obs_keys : sequence of str or None
        Which obs columns to export. If None, all columns are exported.
    centroid_outlier_quantile : float
        Quantile of distances to keep as inliers when computing centroids.
    centroid_min_points : int
        Minimum number of points in a category to compute a centroid.
    force : bool
        If True, overwrite existing files. If False, skip files that already exist.

    Dataset Metadata Parameters
    ---------------------------
    dataset_name : str, optional
        Human-readable name for the dataset (e.g., "Human Lung Cell Atlas").
        If not provided, defaults to the output directory name.
    dataset_description : str, optional
        Description of the dataset.
    dataset_id : str, optional
        Unique identifier for the dataset. If not provided, a filesystem-safe
        version of the dataset_name is used.
    source_name : str, optional
        Name of the data source (e.g., "HLCA Consortium").
    source_url : str, optional
        URL to the data source.
    source_citation : str, optional
        Citation text for the data source.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    obs_binary_dir = out_dir / obs_binary_dirname
    obs_binary_dir.mkdir(parents=True, exist_ok=True)

    # Normalize compression parameter
    if compression is not None and compression <= 0:
        compression = None

    # =========================================================================
    # MULTI-DIMENSIONAL EMBEDDING VALIDATION & PROCESSING
    # =========================================================================
    # Collect all provided embeddings
    embeddings: dict[int, np.ndarray] = {}
    if X_umap_1d is not None:
        embeddings[1] = np.asarray(X_umap_1d, dtype=np.float32)
    if X_umap_2d is not None:
        embeddings[2] = np.asarray(X_umap_2d, dtype=np.float32)
    if X_umap_3d is not None:
        embeddings[3] = np.asarray(X_umap_3d, dtype=np.float32)
    if X_umap_4d is not None:
        # 4D is a hook for future development - raise error for now
        raise NotImplementedError(
            "4D visualization is not yet implemented. "
            "The X_umap_4d parameter is reserved for future development. "
            "Please use X_umap_1d, X_umap_2d, or X_umap_3d for now."
        )

    if not embeddings:
        raise ValueError(
            "At least one dimensional embedding must be provided. "
            "Use X_umap_1d, X_umap_2d, or X_umap_3d."
        )

    # Validate each embedding has correct dimensions
    n_cells = None
    for dim, arr in embeddings.items():
        if arr.ndim != 2:
            raise ValueError(
                f"X_umap_{dim}d must be a 2D array, got shape {arr.shape}."
            )
        if arr.shape[1] != dim:
            raise ValueError(
                f"X_umap_{dim}d must have exactly {dim} columns, got {arr.shape[1]}. "
                f"Shape is {arr.shape}."
            )
        if n_cells is None:
            n_cells = arr.shape[0]
        elif arr.shape[0] != n_cells:
            raise ValueError(
                f"All embeddings must have the same number of cells. "
                f"First embedding has {n_cells} cells, but X_umap_{dim}d has {arr.shape[0]} cells."
            )

    # =========================================================================
    # NORMALIZE EACH EMBEDDING INDEPENDENTLY TO FIT WITHIN [-1, 1] RANGE
    # =========================================================================
    # Each dimensional embedding (1D, 2D, 3D) is normalized independently so that
    # it fills the viewing area optimally. Within each dimension, we use the same
    # scale factor for all axes to preserve aspect ratios.
    #
    # This ensures that switching between dimensions doesn't require manual zoom
    # adjustments - each dimension will fill the view appropriately.

    normalization_info = {}
    for dim, arr in embeddings.items():
        # Find min/max for each axis
        axis_mins = arr.min(axis=0)
        axis_maxs = arr.max(axis=0)
        axis_ranges = axis_maxs - axis_mins

        # Use the maximum range across all axes to preserve aspect ratio
        max_range = float(axis_ranges.max())
        if max_range < 1e-8:
            max_range = 1.0  # Avoid division by zero for degenerate data

        # Center of the bounding box
        center = (axis_mins + axis_maxs) / 2

        # Scale to fit in [-1, 1] based on the max range (preserves aspect ratio)
        scale_factor = 2.0 / max_range
        embeddings[dim] = ((arr - center) * scale_factor).astype(np.float32)

        # Store info for logging
        normalization_info[dim] = {
            'original_range': max_range,
            'center': center.tolist(),
            'scale_factor': scale_factor,
        }

    # Determine primary 3D coords for centroid computation
    # Priority: 3D > 2D (padded) > 1D (padded)
    # Note: 4D support is reserved for future development
    if 3 in embeddings:
        coords3d = embeddings[3]
    elif 2 in embeddings:
        # Pad 2D to 3D with zeros in Z
        coords3d = np.hstack([embeddings[2], np.zeros((n_cells, 1), dtype=np.float32)])
    elif 1 in embeddings:
        # Pad 1D to 3D with zeros in Y and Z
        coords3d = np.hstack([embeddings[1], np.zeros((n_cells, 2), dtype=np.float32)])
    else:
        raise ValueError("No usable embedding for 3D coordinates.")

    # Track available dimensions for metadata
    available_dimensions = sorted(embeddings.keys())

    # Determine default dimension (priority: 3D > 2D > 1D)
    default_dimension = 3 if 3 in embeddings else (2 if 2 in embeddings else 1)

    # Print export settings summary
    print("=" * 60)
    print("Export Settings:")
    print(f"  Output directory: {out_dir}")
    print(f"  Compression: {'gzip level ' + str(compression) if compression else 'disabled'}")
    print(f"  Var (gene) quantization: {str(var_quantization) + '-bit' if var_quantization else 'disabled (float32)'}")
    print(f"  Obs continuous quantization: {str(obs_continuous_quantization) + '-bit' if obs_continuous_quantization else 'disabled (float32)'}")
    print(f"  Obs categorical dtype: {obs_categorical_dtype}")
    print(f"  Available dimensions: {available_dimensions}")
    print(f"  Default dimension: {default_dimension}D")
    print(f"  Coordinate normalization (per-dimension, aspect-ratio preserved):")
    for dim in sorted(normalization_info.keys()):
        info = normalization_info[dim]
        print(f"    {dim}D: range {info['original_range']:.2f} → [-1, 1]")
    print("=" * 60)

    # Validate and convert latent space
    if latent_space is None:
        raise ValueError("latent_space is required for outlier quantile calculation.")
    latent = _to_dense(latent_space).astype(np.float32)
    if latent.shape[0] != n_cells:
        raise ValueError(
            f"Latent space has {latent.shape[0]} cells, but embeddings have {n_cells} cells."
        )

    # Validate obs
    if obs is None:
        raise ValueError("obs DataFrame is required.")
    if len(obs) != n_cells:
        raise ValueError(
            f"obs has {len(obs)} rows, but embeddings have {n_cells} cells."
        )

    # =========================================================================
    # SAVE DIMENSIONAL EMBEDDING FILES
    # =========================================================================
    for dim, arr in embeddings.items():
        dim_filename = f"points_{dim}d.bin"
        dim_path = out_dir / dim_filename
        check_path = Path(str(dim_path) + ".gz") if compression and compression > 0 else dim_path
        if _file_exists_skip(check_path, check_path.name, force):
            pass
        else:
            actual_path = _write_binary(dim_path, arr, compression)
            suffix = " (gzip)" if compression else ""
            print(f"✓ Wrote {dim}D positions ({arr.shape[0]:,} cells × {dim} dims) to {actual_path}{suffix}")

    # =========================================================================
    # SAVE VECTOR FIELDS (OPTIONAL)
    # =========================================================================
    vector_fields_identity: Optional[dict] = None
    if vector_fields:
        if not isinstance(vector_fields, dict):
            raise TypeError("vector_fields must be a dict[str, array]")

        vectors_dir = out_dir / "vectors"
        vectors_dir.mkdir(parents=True, exist_ok=True)

        def _infer_vector_shape(arr: Union[np.ndarray, sparse.spmatrix], name: str) -> tuple[int, np.ndarray]:
            dense = _to_dense(arr)
            dense = np.asarray(dense)
            if dense.ndim == 1:
                dense = dense.reshape(-1, 1)
            if dense.ndim != 2:
                raise ValueError(f"Vector field '{name}' must be 1D or 2D array, got shape {dense.shape}")
            if dense.shape[1] not in (1, 2, 3):
                raise ValueError(f"Vector field '{name}' must have 1/2/3 components, got shape {dense.shape}")
            return int(dense.shape[1]), dense

        def _label_for(field_id: str) -> str:
            base = field_id[:-5] if field_id.endswith("_umap") else field_id
            base = base.replace("_", " ").strip()
            titled = (base[:1].upper() + base[1:]) if base else field_id
            return f"{titled} (UMAP)" if field_id.endswith("_umap") else titled

        # First pass: collect explicit keys (<field>_<dim>d), then fill gaps with implicit keys.
        grouped: dict[str, dict[int, Union[np.ndarray, sparse.spmatrix]]] = {}
        explicit_dims: dict[str, set[int]] = {}

        suffix_re = re.compile(r"^(.*)_([123])d$")
        for name, arr in vector_fields.items():
            if arr is None:
                continue
            key = str(name)
            match = suffix_re.match(key)
            if not match:
                continue
            field_id = match.group(1)
            dim = int(match.group(2))
            explicit_dims.setdefault(field_id, set()).add(dim)
            grouped.setdefault(field_id, {})[dim] = arr

        for name, arr in vector_fields.items():
            if arr is None:
                continue
            key = str(name)
            if suffix_re.match(key):
                continue

            inferred_dim, _dense = _infer_vector_shape(arr, key)
            field_id = key
            if inferred_dim in explicit_dims.get(field_id, set()):
                continue  # clash-safe: explicit wins
            grouped.setdefault(field_id, {})[inferred_dim] = arr

        fields_meta: dict[str, dict] = {}
        gz_suffix = ".gz" if compression else ""

        for field_id, by_dim in grouped.items():
            safe_id = _safe_filename_component(field_id)
            if safe_id != field_id:
                raise ValueError(
                    f"Vector field id '{field_id}' contains unsupported characters. "
                    f"Use '{safe_id}' instead."
                )

            files: dict[str, str] = {}
            dims: list[int] = []

            for dim in sorted(by_dim.keys()):
                if dim not in embeddings:
                    print(
                        f"  ⚠ Skipping vector field '{field_id}' {dim}D: "
                        f"embedding points_{dim}d not provided"
                    )
                    continue

                inferred_dim, dense = _infer_vector_shape(by_dim[dim], f"{field_id}_{dim}d")
                if inferred_dim != dim:
                    raise ValueError(
                        f"Vector field '{field_id}' declared as {dim}D but has shape {dense.shape}"
                    )
                if dense.shape[0] != n_cells:
                    raise ValueError(
                        f"Vector field '{field_id}' {dim}D has {dense.shape[0]} rows, expected {n_cells}"
                    )

                vec = np.asarray(dense, dtype=np.float32)
                scale_factor = float(normalization_info.get(dim, {}).get("scale_factor", 1.0))
                if scale_factor != 1.0:
                    vec *= scale_factor

                filename = f"{field_id}_{dim}d.bin"
                path = vectors_dir / filename
                check_path = Path(str(path) + ".gz") if compression and compression > 0 else path
                if _file_exists_skip(check_path, check_path.name, force):
                    pass
                else:
                    actual_path = _write_binary(path, vec, compression)
                    suffix = " (gzip)" if compression else ""
                    print(
                        f"✓ Wrote vector field '{field_id}' {dim}D "
                        f"({vec.shape[0]:,} cells × {dim} comps) to {actual_path}{suffix}"
                    )

                files[f"{dim}d"] = f"vectors/{filename}{gz_suffix}"
                dims.append(dim)

            if not dims:
                continue

            entry = {
                "label": _label_for(field_id),
                "available_dimensions": dims,
                "default_dimension": max(dims),
                "files": files,
            }
            if field_id.endswith("_umap"):
                entry["basis"] = "umap"
            fields_meta[field_id] = entry

        if fields_meta:
            default_field = "velocity_umap" if "velocity_umap" in fields_meta else sorted(fields_meta.keys())[0]
            vector_fields_identity = {
                "default_field": default_field,
                "fields": fields_meta,
            }


    # Decide which obs columns to export
    if obs_keys is None:
        obs_keys = list(obs.columns)
    else:
        obs_keys = list(obs_keys)
        missing = [k for k in obs_keys if k not in obs.columns]
        if missing:
            raise KeyError(
                f"obs_keys contain columns not in obs: {missing}. "
                f"Available columns: {list(obs.columns)}"
            )

    # Collect lightweight metadata for obs fields (used in dataset identity)
    obs_field_summaries: list[dict] = []
    for key in obs_keys:
        s = obs[key]
        if pd.api.types.is_categorical_dtype(s):
            kind = "category"
        elif pd.api.types.is_bool_dtype(s):
            kind = "category"
        elif pd.api.types.is_numeric_dtype(s):
            kind = "continuous"
        else:
            kind = "category"

        if kind == "continuous":
            if obs_continuous_quantization is None:
                dtype_str = "float32"
            elif obs_continuous_quantization == 8:
                dtype_str = "uint8"
            else:
                dtype_str = "uint16"
            obs_field_summaries.append(
                {
                    "key": str(key),
                    "kind": "continuous",
                    "quantized": obs_continuous_quantization is not None,
                    "quantization_bits": int(obs_continuous_quantization)
                    if obs_continuous_quantization is not None
                    else None,
                    "dtype": dtype_str,
                }
            )
        else:
            cat = s.astype("category")
            categories = [str(c) for c in cat.cat.categories]
            n_categories = len(categories)
            if obs_categorical_dtype == "auto":
                dtype_str = "uint8" if n_categories <= 254 else "uint16"
            elif obs_categorical_dtype == "uint8":
                dtype_str = "uint8"
            else:
                dtype_str = "uint16"
            obs_field_summaries.append(
                {
                    "key": str(key),
                    "kind": "category",
                    "category_count": n_categories,
                    "codes_dtype": dtype_str,
                    "outlier_quantized": obs_continuous_quantization is not None,
                    "outlier_quantization_bits": int(obs_continuous_quantization)
                    if obs_continuous_quantization is not None
                    else None,
                }
            )

    # Check if obs manifest already exists
    obs_manifest_path = out_dir / obs_manifest_filename
    if _file_exists_skip(obs_manifest_path, "obs manifest", force):
        pass
    else:
        # Compact format: separate lists for continuous and categorical fields
        obs_continuous_fields: list = []
        obs_categorical_fields: list = []
        # Track dtype info for schema (will use first encountered)
        continuous_dtype_info: dict = {}
        categorical_dtype_info: dict = {}

        for key in obs_keys:
            s = obs[key]
            safe_key = _safe_filename_component(key)
            
            # Decide kind: continuous vs categorical
            if pd.api.types.is_categorical_dtype(s):
                kind = "category"
            elif pd.api.types.is_bool_dtype(s):
                kind = "category"
            elif pd.api.types.is_numeric_dtype(s):
                kind = "continuous"
            else:
                kind = "category"

            if kind == "continuous":
                values = pd.to_numeric(s, errors="coerce").to_numpy(dtype=np.float32, copy=False)
                if values.shape[0] != n_cells:
                    raise ValueError(
                        f"Length mismatch for obs['{key}']: {values.shape[0]} vs {n_cells}"
                    )

                # Apply quantization if requested
                if obs_continuous_quantization is not None:
                    quantized, min_val, max_val, scale, stats = _quantize_continuous(
                        values, bits=obs_continuous_quantization, field_name=key
                    )
                    _report_quantization_stats(key, stats, "obs continuous")
                    
                    if obs_continuous_quantization == 8:
                        dtype_str = "uint8"
                        ext = "u8"
                    else:
                        dtype_str = "uint16"
                        ext = "u16"
                    
                    value_fname = f"{safe_key}.values.{ext}"
                    value_path = obs_binary_dir / value_fname
                    actual_path = _write_binary(value_path, quantized, compression)
                    
                    # Adjust path in manifest if compressed
                    manifest_path = f"{obs_binary_dirname}/{value_fname}"
                    if compression:
                        manifest_path += ".gz"
                    
                    # Compact format: [key, minValue, maxValue]
                    obs_continuous_fields.append([key, min_val, max_val])
                    if not continuous_dtype_info:
                        continuous_dtype_info["ext"] = ext
                        continuous_dtype_info["dtype"] = dtype_str
                        continuous_dtype_info["quantized"] = True
                        continuous_dtype_info["quantizationBits"] = obs_continuous_quantization
                else:
                    # Full precision
                    value_fname = f"{safe_key}.values.f32"
                    value_path = obs_binary_dir / value_fname
                    actual_path = _write_binary(value_path, values, compression)

                    # Compact format: [key]
                    obs_continuous_fields.append([key])
                    if not continuous_dtype_info:
                        continuous_dtype_info["ext"] = "f32"
                        continuous_dtype_info["dtype"] = "float32"
                        continuous_dtype_info["quantized"] = False

            else:
                # Categorical
                cat = s.astype("category")
                categories = [str(c) for c in cat.cat.categories]
                codes = cat.cat.codes.to_numpy(dtype=np.int32)  # -1 for NaN

                if codes.shape[0] != n_cells:
                    raise ValueError(
                        f"Length mismatch for obs['{key}']: {codes.shape[0]} vs {n_cells}"
                    )

                n_categories = len(categories)
                
                # Select dtype based on settings
                if obs_categorical_dtype == "auto":
                    dtype, missing_value = _select_category_dtype(n_categories)
                elif obs_categorical_dtype == "uint8":
                    if n_categories > 254:
                        raise ValueError(
                            f"Field '{key}' has {n_categories} categories, "
                            f"but uint8 can only hold 254. Use 'auto' or 'uint16'."
                        )
                    dtype, missing_value = np.uint8, 255
                else:  # uint16
                    dtype, missing_value = np.uint16, 65535
                
                codes_typed = np.full(n_cells, missing_value, dtype=dtype)
                valid_mask = codes >= 0
                codes_typed[valid_mask] = codes[valid_mask].astype(dtype)
                
                if dtype == np.uint8:
                    codes_fname = f"{safe_key}.codes.u8"
                    dtype_str = "uint8"
                else:
                    codes_fname = f"{safe_key}.codes.u16"
                    dtype_str = "uint16"
                
                codes_path = obs_binary_dir / codes_fname
                actual_path = _write_binary(codes_path, codes_typed, compression)
                
                manifest_codes_path = f"{obs_binary_dirname}/{codes_fname}"
                if compression:
                    manifest_codes_path += ".gz"

                # Compute centroids for all available dimensions
                if centroid_outlier_quantile is None:
                    centroids_by_dim = {dim: [] for dim in embeddings.keys()}
                else:
                    centroids_by_dim = _compute_centroids_for_all_dimensions(
                        embeddings,
                        codes,
                        categories,
                        outlier_quantile=centroid_outlier_quantile,
                        min_points=centroid_min_points,
                    )

                # Compute per-cell outlier quantiles based on latent space
                outlier_quantiles = _compute_latent_space_quantiles(
                    latent=latent,
                    codes=codes,
                    categories=categories,
                    min_points=centroid_min_points,
                )
                
                # Quantize outlier quantiles (they're always 0-1)
                if obs_continuous_quantization is not None:
                    oq_quantized, oq_min, oq_max, oq_scale, oq_stats = _quantize_continuous(
                        outlier_quantiles, bits=obs_continuous_quantization, field_name=f"{key}_outliers"
                    )
                    _report_quantization_stats(f"{key}_outliers", oq_stats, "outlier quantiles")
                    
                    if obs_continuous_quantization == 8:
                        oq_dtype_str = "uint8"
                        oq_ext = "u8"
                    else:
                        oq_dtype_str = "uint16"
                        oq_ext = "u16"
                    
                    outlier_fname = f"{safe_key}.outliers.{oq_ext}"
                    outlier_path = obs_binary_dir / outlier_fname
                    _write_binary(outlier_path, oq_quantized, compression)
                    
                    manifest_outlier_path = f"{obs_binary_dirname}/{outlier_fname}"
                    if compression:
                        manifest_outlier_path += ".gz"
                    
                    # Compact format: [key, categories, codesDtype, codesMissingValue, centroidsByDim, outlierMinValue, outlierMaxValue]
                    # centroidsByDim is a dict keyed by dimension: {"1": [...], "2": [...], "3": [...]}
                    centroids_serializable = {str(dim): cents for dim, cents in centroids_by_dim.items()}
                    obs_categorical_fields.append([
                        key, categories, dtype_str, int(missing_value), centroids_serializable, oq_min, oq_max
                    ])
                    if not categorical_dtype_info:
                        categorical_dtype_info["codesExt"] = "u8" if dtype == np.uint8 else "u16"
                        categorical_dtype_info["outlierExt"] = oq_ext
                        categorical_dtype_info["outlierDtype"] = oq_dtype_str
                        categorical_dtype_info["outlierQuantized"] = True
                else:
                    # Full precision outliers
                    outlier_fname = f"{safe_key}.outliers.f32"
                    outlier_path = obs_binary_dir / outlier_fname
                    _write_binary(outlier_path, outlier_quantiles.astype(np.float32), compression)

                    # Compact format: [key, categories, codesDtype, codesMissingValue, centroidsByDim]
                    # centroidsByDim is a dict keyed by dimension: {"1": [...], "2": [...], "3": [...]}
                    centroids_serializable = {str(dim): cents for dim, cents in centroids_by_dim.items()}
                    obs_categorical_fields.append([
                        key, categories, dtype_str, int(missing_value), centroids_serializable
                    ])
                    if not categorical_dtype_info:
                        categorical_dtype_info["codesExt"] = "u8" if dtype == np.uint8 else "u16"
                        categorical_dtype_info["outlierExt"] = "f32"
                        categorical_dtype_info["outlierDtype"] = "float32"
                        categorical_dtype_info["outlierQuantized"] = False

        # Build compact manifest with schemas
        gz_suffix = ".gz" if compression else ""

        obs_schemas = {}
        if continuous_dtype_info:
            obs_schemas["continuous"] = {
                "pathPattern": f"{obs_binary_dirname}/{{key}}.values.{continuous_dtype_info['ext']}{gz_suffix}",
                "ext": continuous_dtype_info["ext"],
                "dtype": continuous_dtype_info["dtype"],
                "quantized": continuous_dtype_info.get("quantized", False),
            }
            if continuous_dtype_info.get("quantized"):
                obs_schemas["continuous"]["quantizationBits"] = continuous_dtype_info["quantizationBits"]

        if categorical_dtype_info:
            obs_schemas["categorical"] = {
                "codesPathPattern": f"{obs_binary_dirname}/{{key}}.codes.{{ext}}{gz_suffix}",
                "outlierPathPattern": f"{obs_binary_dirname}/{{key}}.outliers.{categorical_dtype_info['outlierExt']}{gz_suffix}",
                "outlierExt": categorical_dtype_info["outlierExt"],
                "outlierDtype": categorical_dtype_info["outlierDtype"],
                "outlierQuantized": categorical_dtype_info.get("outlierQuantized", False),
            }

        obs_manifest_payload = {
            "_format": MANIFEST_FORMAT_VERSION,
            "n_points": int(n_cells),
            "centroid_outlier_quantile": float(centroid_outlier_quantile)
            if centroid_outlier_quantile is not None
            else None,
            "latent_key": "latent_space",
            "compression": compression if compression else None,
            "_obsSchemas": obs_schemas,
            "_continuousFields": obs_continuous_fields,
            "_categoricalFields": obs_categorical_fields,
        }
        obs_manifest_path.write_text(json.dumps(obs_manifest_payload), encoding="utf-8")

        total_fields = len(obs_continuous_fields) + len(obs_categorical_fields)
        print(
            f"✓ Wrote obs manifest ({total_fields} fields: {len(obs_continuous_fields)} continuous, "
            f"{len(obs_categorical_fields)} categorical) to {obs_manifest_path} "
            f"with binaries in {obs_binary_dirname}/"
        )

    # Process gene expression if provided
    genes_to_export: list[str] = []
    if gene_expression is not None:
        if var is None:
            raise ValueError("var DataFrame must be provided when gene_expression is given.")

        gene_expr_is_sparse = sparse.issparse(gene_expression)
        gene_expression_for_export = (
            gene_expression.tocsc() if gene_expr_is_sparse and not sparse.isspmatrix_csc(gene_expression) else gene_expression
        )

        n_expr_cells = gene_expression_for_export.shape[0]
        n_genes = gene_expression_for_export.shape[1]

        if n_expr_cells != n_cells:
            raise ValueError(
                f"gene_expression has {n_expr_cells} cells, but embeddings have {n_cells} cells."
            )

        if len(var) != n_genes:
            raise ValueError(
                f"var has {len(var)} rows, but gene_expression has {n_genes} genes."
            )

        if var_gene_id_column == "index" or var_gene_id_column is None:
            all_gene_ids = var.index.astype(str).tolist()
        else:
            if var_gene_id_column not in var.columns:
                raise KeyError(
                    f"var_gene_id_column '{var_gene_id_column}' not found in var. "
                    f"Available columns: {list(var.columns)}"
                )
            all_gene_ids = var[var_gene_id_column].astype(str).tolist()

        gene_id_to_idx = {gid: idx for idx, gid in enumerate(all_gene_ids)}

        if gene_identifiers is None:
            genes_to_export = all_gene_ids
        else:
            genes_to_export = list(gene_identifiers)
            missing_genes = [g for g in genes_to_export if g not in gene_id_to_idx]
            if missing_genes:
                print(f"⚠ Warning: {len(missing_genes)} gene identifiers not found in var: {missing_genes[:5]}...")
            genes_to_export = [g for g in genes_to_export if g in gene_id_to_idx]

        var_manifest_path = out_dir / var_manifest_filename
        if _file_exists_skip(var_manifest_path, "var manifest", force):
            pass
        else:
            var_binary_dir = out_dir / var_binary_dirname
            var_binary_dir.mkdir(parents=True, exist_ok=True)

            var_manifest_fields: list[dict] = []

            # Track problematic genes for aggregated reporting
            genes_with_nan: list[str] = []
            genes_with_inf: list[str] = []
            genes_all_invalid: list[str] = []

            for gene_id in tqdm.tqdm(genes_to_export, desc="Exporting genes"):
                gene_idx = gene_id_to_idx[gene_id]
                safe_gene_id = _safe_filename_component(gene_id)

                if gene_expr_is_sparse:
                    col = gene_expression_for_export.getcol(gene_idx).toarray().flatten()
                else:
                    col = gene_expression_for_export[:, gene_idx]

                values = np.asarray(col, dtype=np.float32)

                if values.shape[0] != n_cells:
                    raise ValueError(
                        f"Gene '{gene_id}' expression length mismatch: {values.shape[0]} vs {n_cells}"
                    )

                # Apply quantization if requested
                if var_quantization is not None:
                    quantized, min_val, max_val, scale, stats = _quantize_continuous(
                        values, bits=var_quantization, field_name=gene_id
                    )

                    # Track aggregated stats for genes (don't spam per-gene)
                    if stats.get("n_nan", 0) > 0:
                        genes_with_nan.append(gene_id)
                    if stats.get("n_inf", 0) > 0:
                        genes_with_inf.append(gene_id)
                    if stats.get("warning") == "all_invalid":
                        genes_all_invalid.append(gene_id)

                    if var_quantization == 8:
                        dtype_str = "uint8"
                        ext = "u8"
                    else:
                        dtype_str = "uint16"
                        ext = "u16"

                    value_fname = f"{safe_gene_id}.values.{ext}"
                    value_path = var_binary_dir / value_fname
                    _write_binary(value_path, quantized, compression)

                    manifest_path = f"{var_binary_dirname}/{value_fname}"
                    if compression:
                        manifest_path += ".gz"

                    # Compact format: [key, minValue, maxValue]
                    var_manifest_fields.append([gene_id, min_val, max_val])
                else:
                    # Full precision
                    value_fname = f"{safe_gene_id}.values.f32"
                    value_path = var_binary_dir / value_fname
                    _write_binary(value_path, values, compression)

                    # Compact format: [key] for non-quantized
                    var_manifest_fields.append([gene_id])

            # Report aggregated gene stats
            if genes_with_nan:
                print(f"  ⚠ {len(genes_with_nan)} genes contain NaN values (mapped to missing marker)")
                if len(genes_with_nan) <= 10:
                    print(f"    Genes: {', '.join(genes_with_nan)}")
                else:
                    print(f"    First 10: {', '.join(genes_with_nan[:10])}...")
            if genes_with_inf:
                print(f"  ⚠ {len(genes_with_inf)} genes contain Inf values (mapped to missing marker)")
                if len(genes_with_inf) <= 10:
                    print(f"    Genes: {', '.join(genes_with_inf)}")
            if genes_all_invalid:
                print(f"  ⚠ WARNING: {len(genes_all_invalid)} genes have NO valid values (all NaN/Inf)")
                if len(genes_all_invalid) <= 10:
                    print(f"    Genes: {', '.join(genes_all_invalid)}")

            # Build compact manifest with schema
            gz_suffix = ".gz" if compression else ""
            if var_quantization is not None:
                ext = "u8" if var_quantization == 8 else "u16"
                dtype_str = "uint8" if var_quantization == 8 else "uint16"
                var_schema = {
                    "kind": "continuous",
                    "pathPattern": f"{var_binary_dirname}/{{key}}.values.{ext}{gz_suffix}",
                    "ext": ext,
                    "dtype": dtype_str,
                    "quantized": True,
                    "quantizationBits": var_quantization,
                }
            else:
                var_schema = {
                    "kind": "continuous",
                    "pathPattern": f"{var_binary_dirname}/{{key}}.values.f32{gz_suffix}",
                    "ext": "f32",
                    "dtype": "float32",
                    "quantized": False,
                }

            var_manifest_payload = {
                "_format": MANIFEST_FORMAT_VERSION,
                "n_points": int(n_cells),
                "var_gene_id_column": var_gene_id_column,
                "compression": compression if compression else None,
                "quantization": var_quantization,
                "_varSchema": var_schema,
                "fields": var_manifest_fields,
            }
            var_manifest_path.write_text(json.dumps(var_manifest_payload), encoding="utf-8")

            compression_info = f", gzip level {compression}" if compression else ""
            quant_info = f", {var_quantization}-bit quantized" if var_quantization else ""
            print(
                f"✓ Wrote var manifest ({len(var_manifest_fields)} genes{quant_info}{compression_info}) "
                f"to {var_manifest_path}"
            )
    else:
        print("INFO: No gene expression data provided, skipping var export.")

    # Process connectivity data if provided
    # GPU-optimized edge format for instanced rendering
    connectivity_meta = {
        "n_edges": None,
        "max_neighbors": None,
        "index_dtype": None,
    }
    if connectivities is not None:
        connectivity_manifest_path = out_dir / connectivity_manifest_filename
        if connectivities.shape[0] != n_cells or connectivities.shape[1] != n_cells:
            raise ValueError(
                f"Connectivity matrix shape {connectivities.shape} does not match "
                f"number of cells {n_cells}."
            )

        # Determine optimal dtype based on cell count up front (used in hash even if we skip writing)
        if n_cells <= 65535:
            index_dtype = np.uint16
            index_dtype_str = "uint16"
            index_bytes = 2
        elif n_cells <= 4294967295:
            index_dtype = np.uint32
            index_dtype_str = "uint32"
            index_bytes = 4
        else:
            index_dtype = np.uint64
            index_dtype_str = "uint64"
            index_bytes = 8

        approx_edges = None
        if sparse.issparse(connectivities):
            approx_edges = int(connectivities.nnz // 2)
        else:
            approx_edges = int(np.count_nonzero(connectivities) // 2)
        connectivity_meta["n_edges"] = approx_edges
        connectivity_meta["index_dtype"] = index_dtype_str

        if _file_exists_skip(connectivity_manifest_path, "connectivity manifest", force):
            pass
        else:
            connectivity_binary_dir = out_dir / connectivity_binary_dirname
            connectivity_binary_dir.mkdir(parents=True, exist_ok=True)

            if not sparse.isspmatrix_csr(connectivities):
                connectivities = sparse.csr_matrix(connectivities)

            # Symmetrize and binarize the connectivity matrix
            connectivities_sym = connectivities + connectivities.T
            connectivities_sym.data[:] = 1
            connectivities_csr = connectivities_sym.tocsr()

            indptr = connectivities_csr.indptr
            indices = connectivities_csr.indices

            # Extract unique edges (src < dst to avoid duplicates)
            # Using vectorized operations for speed with large datasets
            print(f"  Extracting unique edges from {n_cells:,} cells...")

            edge_sources = []
            edge_destinations = []
            max_neighbors_found = 0

            # Process in chunks for memory efficiency with very large datasets
            chunk_size = 100000
            for chunk_start in range(0, n_cells, chunk_size):
                chunk_end = min(chunk_start + chunk_size, n_cells)
                for cell_idx in range(chunk_start, chunk_end):
                    start = indptr[cell_idx]
                    end = indptr[cell_idx + 1]
                    neighbor_count = end - start
                    if neighbor_count > max_neighbors_found:
                        max_neighbors_found = neighbor_count

                    for j in range(start, end):
                        neighbor_idx = indices[j]
                        # Only keep edges where src < dst (avoid duplicates)
                        if cell_idx < neighbor_idx:
                            edge_sources.append(cell_idx)
                            edge_destinations.append(neighbor_idx)

            edge_sources = np.array(edge_sources, dtype=index_dtype)
            edge_destinations = np.array(edge_destinations, dtype=index_dtype)
            n_unique_edges = len(edge_sources)
            connectivity_meta["n_edges"] = int(n_unique_edges)
            connectivity_meta["max_neighbors"] = int(max_neighbors_found)
            connectivity_meta["index_dtype"] = index_dtype_str

            print(f"  Found {n_unique_edges:,} unique edges, max {max_neighbors_found} neighbors/cell")

            # Sort edges by source, then by destination for optimal gzip compression
            print(f"  Sorting edges for optimal compression...")
            sort_idx = np.lexsort((edge_destinations, edge_sources))
            edge_sources = edge_sources[sort_idx]
            edge_destinations = edge_destinations[sort_idx]

            # Write binary files (column-separated for better compression)
            sources_fname = "edges.src.bin"
            dests_fname = "edges.dst.bin"
            sources_path = connectivity_binary_dir / sources_fname
            dests_path = connectivity_binary_dir / dests_fname

            _write_binary(sources_path, edge_sources, compression)
            _write_binary(dests_path, edge_destinations, compression)

            manifest_sources = f"{connectivity_binary_dirname}/{sources_fname}"
            manifest_dests = f"{connectivity_binary_dirname}/{dests_fname}"
            if compression:
                manifest_sources += ".gz"
                manifest_dests += ".gz"

            # Write manifest
            connectivity_manifest_payload = {
                "format": "edge_pairs",
                "n_cells": int(n_cells),
                "n_edges": int(n_unique_edges),
                "max_neighbors": int(max_neighbors_found),
                "index_bytes": index_bytes,
                "index_dtype": index_dtype_str,
                "sourcesPath": manifest_sources,
                "destinationsPath": manifest_dests,
                "compression": compression if compression else None,
            }

            connectivity_manifest_path.write_text(
                json.dumps(connectivity_manifest_payload), encoding="utf-8"
            )

            print(
                f"✓ Wrote connectivity ({n_unique_edges:,} edges, "
                f"max {max_neighbors_found} neighbors/cell, {index_dtype_str}) "
                f"to {connectivity_binary_dir}"
            )
    else:
        print("INFO: No connectivity data provided, skipping connectivity export.")

    # =========================================================================
    # Generate dataset_identity.json (metadata for multi-dataset support)
    # =========================================================================
    identity_path = out_dir / "dataset_identity.json"

    # Try to import version from package
    try:
        from cellucid import __version__ as cellucid_version
    except ImportError:
        cellucid_version = "unknown"

    # Determine dataset ID and name
    if dataset_id is None:
        if dataset_name:
            dataset_id = _safe_filename_component(dataset_name)
        else:
            dataset_id = _safe_filename_component(out_dir.name)

    if dataset_name is None:
        dataset_name = out_dir.name

    # Count genes
    n_genes = 0
    if gene_expression is not None and var is not None:
        if sparse.issparse(gene_expression):
            n_genes = gene_expression.shape[1]
        else:
            n_genes = np.asarray(gene_expression).shape[1]

    # Count obs fields from the obs_field_summaries list populated earlier (lines 528-582)
    # This works regardless of whether the obs manifest was written or skipped
    n_obs_fields = len(obs_field_summaries)
    n_categorical_fields = sum(1 for f in obs_field_summaries if f.get("kind") == "category")
    n_continuous_fields = sum(1 for f in obs_field_summaries if f.get("kind") == "continuous")

    # Validate that field counts add up (sanity check)
    if n_obs_fields != n_categorical_fields + n_continuous_fields:
        print(
            f"  ⚠ WARNING: Field count mismatch! "
            f"n_obs_fields={n_obs_fields}, categorical={n_categorical_fields}, continuous={n_continuous_fields}"
        )
        # Debug: show fields with unexpected kinds
        unexpected = [f for f in obs_field_summaries if f.get("kind") not in ("category", "continuous")]
        if unexpected:
            print(f"    Fields with unexpected kind: {[f.get('key') for f in unexpected]}")

    # Convert obs_field_summaries to the format expected in dataset_identity.json
    # (simplified version with just key, kind, and n_categories for categorical)
    identity_obs_fields = []
    for field_info in obs_field_summaries:
        entry = {
            "key": field_info["key"],
            "kind": field_info["kind"]
        }
        if field_info["kind"] == "category" and "category_count" in field_info:
            entry["n_categories"] = field_info["category_count"]
        identity_obs_fields.append(entry)

    # Build source info if provided
    source_info = None
    if source_name or source_url or source_citation:
        source_info = {}
        if source_name:
            source_info["name"] = source_name
        if source_url:
            source_info["url"] = source_url
        if source_citation:
            source_info["citation"] = source_citation

    # Build export settings
    export_settings = {
        "compression": compression if compression else None,
        "var_quantization": var_quantization,
        "obs_continuous_quantization": obs_continuous_quantization,
        "obs_categorical_dtype": obs_categorical_dtype
    }

    # Build embeddings metadata
    gz_suffix = ".gz" if compression else ""
    embeddings_meta = {
        "available_dimensions": available_dimensions,
        "default_dimension": default_dimension,
        "files": {}
    }
    for dim in available_dimensions:
        embeddings_meta["files"][f"{dim}d"] = f"points_{dim}d.bin{gz_suffix}"

    # Build identity payload
    identity_payload = {
        "version": 2,  # Bumped version for multi-dimensional support
        "id": dataset_id,
        "name": dataset_name,
        "description": dataset_description or "",
        "created_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "cellucid_data_version": cellucid_version,
        "stats": {
            "n_cells": int(n_cells),
            "n_genes": int(n_genes),
            "n_obs_fields": int(n_obs_fields),
            "n_categorical_fields": int(n_categorical_fields),
            "n_continuous_fields": int(n_continuous_fields),
            "has_connectivity": connectivity_meta.get("n_edges") is not None,
            "n_edges": connectivity_meta.get("n_edges")
        },
        "embeddings": embeddings_meta,
        "obs_fields": identity_obs_fields,
        "export_settings": export_settings
    }

    if source_info:
        identity_payload["source"] = source_info

    if vector_fields_identity:
        identity_payload["vector_fields"] = vector_fields_identity

    identity_path.write_text(json.dumps(identity_payload, indent=2), encoding="utf-8")
    print(f"✓ Wrote dataset identity to {identity_path}")


def generate_datasets_manifest(
    exports_dir: Union[str, Path] = DEFAULT_EXPORT_DIR,
    output_filename: str = "datasets.json",
    default_dataset: Optional[str] = None,
) -> Optional[Path]:
    """
    Scan an exports directory for dataset_identity.json files and generate a datasets.json manifest.

    This utility helps maintain the datasets.json manifest file that the frontend uses
    to discover available demo datasets. Run this after adding or removing datasets.

    Parameters
    ----------
    exports_dir : Path or str
        Directory containing dataset subdirectories (default: exports/).
    output_filename : str
        Name of the output manifest file (default: datasets.json).
    default_dataset : str or None
        ID of the default dataset. If None, uses the first dataset found.

    Returns
    -------
    Path
        Path to the generated datasets.json file.

    Example
    -------
    >>> from cellucid.prepare_data import generate_datasets_manifest
    >>> generate_datasets_manifest("./exports", default_dataset="my_dataset")
    """
    exports_dir = Path(exports_dir)
    if not exports_dir.exists():
        raise FileNotFoundError(f"Exports directory not found: {exports_dir}")

    print(f"Scanning {exports_dir} for datasets...")

    datasets = []

    # Scan subdirectories for dataset_identity.json
    for subdir in sorted(exports_dir.iterdir()):
        if not subdir.is_dir():
            continue

        identity_file = subdir / "dataset_identity.json"
        if not identity_file.exists():
            print(f"  ⚠ Skipping {subdir.name}: no dataset_identity.json")
            continue

        try:
            identity = json.loads(identity_file.read_text(encoding="utf-8"))

            dataset_entry = {
                "id": identity.get("id", subdir.name),
                "path": f"{subdir.name}/",
                "name": identity.get("name", subdir.name),
            }

            # Include quick stats for display in dropdown
            stats = identity.get("stats", {})
            if stats.get("n_cells"):
                dataset_entry["n_cells"] = stats["n_cells"]
            if stats.get("n_genes"):
                dataset_entry["n_genes"] = stats["n_genes"]

            datasets.append(dataset_entry)
            print(f"  ✓ Found dataset: {dataset_entry['name']} ({dataset_entry['id']})")

        except json.JSONDecodeError as e:
            print(f"  ⚠ Skipping {subdir.name}: invalid JSON - {e}")
        except Exception as e:
            print(f"  ⚠ Skipping {subdir.name}: {e}")

    if not datasets:
        print("  No datasets found!")
        return None

    # Determine default dataset
    if default_dataset:
        if not any(d["id"] == default_dataset for d in datasets):
            print(f"  ⚠ Warning: default_dataset '{default_dataset}' not found, using first dataset")
            default_dataset = datasets[0]["id"]
    else:
        default_dataset = datasets[0]["id"]

    # Build manifest
    manifest = {
        "version": 1,
        "default": default_dataset,
        "datasets": datasets
    }

    # Write manifest
    manifest_path = exports_dir / output_filename
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"✓ Wrote datasets manifest with {len(datasets)} datasets to {manifest_path}")
    print(f"  Default dataset: {default_dataset}")

    return manifest_path
