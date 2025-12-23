"""
AnnData Adapter for Cellucid

Provides a server-side adapter that reads AnnData objects (in memory, backed h5ad,
or zarr stores) and serves data in the same format that the Cellucid web viewer
expects. This allows users to visualize AnnData directly without running
prepare first.

Key features:
- Works with in-memory AnnData, backed h5ad files, and zarr stores
- Handles sparse matrices (CSR/CSC) transparently with automatic format conversion
- Auto-detects UMAP dimensions from explicit obsm keys (X_umap_1d, X_umap_2d, X_umap_3d) or X_umap
- Computes centroids and outlier quantiles on-demand with caching
- Lazy loading: gene expression and obs data are loaded on-demand
- No quantization (full float32 precision, gzip compression for network transfer)

Lazy Loading Behavior:
----------------------
- h5ad (backed='r'): Uses HDF5 memory-mapping. Arrays are loaded only when accessed.
  Gene expression columns are fetched individually, minimizing memory usage.
- zarr: Arrays are stored as separate files. Reading is inherently lazy - only
  requested array chunks are loaded into memory. Note: zarr does not support
  backed mode like h5ad; the AnnData structure is loaded but data arrays are lazy.
- in-memory: All data is already in RAM. No lazy loading.

Memory Management:
-----------------
The adapter maintains several caches for performance:
- Embedding cache: Normalized UMAP coordinates (one per dimension)
- Centroid cache: Computed label centroids
- Gene expression LRU cache: Recently accessed gene columns (max 100 genes)
- CSC cache: For CSR matrices, a CSC copy is created for efficient column access

IMPORTANT: Always use the context manager or call close() when done:
    with AnnDataAdapter.from_file("data.h5ad") as adapter:
        # use adapter
    # resources are automatically released

Usage:
    from cellucid.anndata_adapter import AnnDataAdapter

    # From h5ad file (lazy loading via backed mode)
    with AnnDataAdapter.from_file("/path/to/data.h5ad") as adapter:
        manifest = adapter.get_obs_manifest()

    # From zarr store (lazy array access)
    adapter = AnnDataAdapter.from_file("/path/to/data.zarr")

    # From in-memory AnnData
    adapter = AnnDataAdapter(adata)

    # Use with server
    from cellucid import show_anndata
    show_anndata(adata)  # or show_anndata("/path/to/data.h5ad")
"""

from __future__ import annotations

import gzip
import json
import logging
import re
import weakref
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Optional, Union

import numpy as np
import pandas as pd
from scipy import sparse

if TYPE_CHECKING:
    import anndata

logger = logging.getLogger("cellucid.anndata_adapter")

# Import shared utilities from prepare_data to avoid duplication
# If running standalone, define locally
try:
    from .prepare_data import _safe_filename_component
except ImportError:
    def _safe_filename_component(name: str) -> str:
        """Return a filesystem-friendly version of a field key."""
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name))
        safe = safe.strip("._")
        return safe or "field"


class LRUCache:
    """
    Simple LRU cache with O(1) operations using OrderedDict.

    WARNING: This class is NOT thread-safe. The AnnDataAdapter is designed
    for single-threaded use within an HTTP server (each request is handled
    sequentially). If you need thread-safe caching, wrap operations in
    external locks or use a thread-safe cache implementation.

    The HTTP server uses a single adapter instance per server, and Python's
    GIL provides some level of safety for simple operations, but concurrent
    access from multiple threads may cause race conditions.
    """

    def __init__(self, max_size: int = 100):
        if max_size < 1:
            raise ValueError("max_size must be at least 1")
        self._cache: OrderedDict[str, Any] = OrderedDict()
        self._max_size = max_size

    def get(self, key: str) -> Optional[Any]:
        """Get item and move to end (most recently used)."""
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        return None

    def put(self, key: str, value: Any) -> None:
        """Add item, evicting oldest if at capacity."""
        if key in self._cache:
            self._cache.move_to_end(key)
        else:
            if len(self._cache) >= self._max_size:
                self._cache.popitem(last=False)  # Remove oldest
        self._cache[key] = value

    def __contains__(self, key: str) -> bool:
        return key in self._cache

    def clear(self) -> None:
        """Clear all cached items and release memory."""
        self._cache.clear()

    def __len__(self) -> int:
        return len(self._cache)

    @property
    def max_size(self) -> int:
        """Maximum cache size."""
        return self._max_size


def _to_dense_1d(arr: Union[np.ndarray, sparse.spmatrix]) -> np.ndarray:
    """Convert sparse vector to dense numpy array."""
    if sparse.issparse(arr):
        return np.asarray(arr.toarray()).flatten()
    arr = np.asarray(arr)
    if arr.ndim > 1:
        return arr.flatten()
    return arr


def _to_dense_2d(arr: Union[np.ndarray, sparse.spmatrix]) -> np.ndarray:
    """Convert sparse matrix to dense numpy array."""
    if sparse.issparse(arr):
        return np.asarray(arr.toarray())
    return np.asarray(arr)


class AnnDataAdapter:
    """
    Adapter that wraps AnnData and provides data in Cellucid format.

    This adapter generates all the data that would normally be created by
    prepare, but reads directly from AnnData without creating
    intermediate files. This is slower but more convenient for interactive use.
    """

    def __init__(
        self,
        adata: "anndata.AnnData",
        latent_key: Optional[str] = None,
        gene_id_column: str = "index",
        normalize_embeddings: bool = True,
        centroid_outlier_quantile: float = 0.95,
        centroid_min_points: int = 10,
        dataset_name: Optional[str] = None,
        dataset_id: Optional[str] = None,
    ):
        """
        Initialize the adapter.

        Parameters
        ----------
        adata : AnnData
            AnnData object to adapt. Can be in-memory or backed (h5ad file).
        latent_key : str, optional
            Key in obsm for latent space used for outlier quantile calculation.
            If None, attempts to find: 'X_pca', 'X_scvi', 'scanvi', 'scvi', or first obsm.
        gene_id_column : str
            Column in var for gene identifiers. Use "index" for var.index.
        normalize_embeddings : bool
            If True, normalize embeddings to [-1, 1] range (recommended).
        centroid_outlier_quantile : float
            Quantile for outlier removal in centroid computation.
        centroid_min_points : int
            Minimum points per category for centroid computation.
        dataset_name : str, optional
            Human-readable dataset name.
        dataset_id : str, optional
            Dataset identifier.
        """
        self.adata = adata
        self.latent_key = latent_key
        self.gene_id_column = gene_id_column
        self.normalize_embeddings = normalize_embeddings

        # Validate centroid parameters
        if not (0.5 < centroid_outlier_quantile < 1.0):
            logger.warning(
                f"centroid_outlier_quantile={centroid_outlier_quantile} is outside "
                f"recommended range (0.5, 1.0). Using default 0.95."
            )
            centroid_outlier_quantile = 0.95
        self.centroid_outlier_quantile = centroid_outlier_quantile

        if centroid_min_points < 1:
            logger.warning(
                f"centroid_min_points={centroid_min_points} must be positive. Using default 10."
            )
            centroid_min_points = 10
        self.centroid_min_points = centroid_min_points

        # Dataset metadata
        self.dataset_name = dataset_name or "AnnData Dataset"
        self.dataset_id = dataset_id or _safe_filename_component(self.dataset_name)

        # Source tracking (set by from_file, defaults for in-memory)
        self._source_path: Optional[str] = None
        self._source_type: str = "memory"
        self._is_lazy: bool = False  # Set by from_file for backed h5ad/zarr

        # Caches for computed data
        self._embedding_cache: dict[int, np.ndarray] = {}
        # Per-dimension normalization info for embeddings (center + scale_factor).
        # Used to scale vector fields into the same normalized space as points.
        self._embedding_norm: dict[int, dict[str, Any]] = {}

        # Optional per-cell vector fields (e.g. velocity, CellRank drift).
        # Structure matches dataset_identity.json: {"default_field": ..., "fields": {...}}.
        self._vector_fields_metadata: Optional[dict[str, Any]] = None
        self._vector_field_cache: dict[tuple[str, int], np.ndarray] = {}
        self._centroid_cache: dict[str, dict] = {}
        self._outlier_quantile_cache: dict[str, np.ndarray] = {}
        self._latent_space: Optional[np.ndarray] = None
        self._gene_ids_cache: Optional[list[str]] = None
        self._gene_id_to_idx_cache: Optional[dict[str, int]] = None

        # Connectivity cache: (sources, destinations, n_edges, max_neighbors)
        self._connectivity_cache: Optional[tuple[np.ndarray, np.ndarray, int, int]] = None

        # CSC cache for efficient column access on CSR matrices
        # Note: This doubles memory for sparse matrices, but is necessary for O(1) gene access
        # Set to None initially; created lazily only when needed
        self._X_csc_cache: Optional[sparse.csc_matrix] = None

        # LRU cache for gene expression values using O(1) OrderedDict
        self._gene_expression_cache = LRUCache(max_size=100)

        # UMAP embedding key resolution (dimension -> obsm key)
        self._umap_embedding_key_by_dim: dict[int, str] = {}
        # Optional metadata for UI notifications about X_umap aliasing/clashes
        self._umap_resolution: Optional[dict[str, Any]] = None

        # Track if adapter has been closed
        self._closed = False

        # Auto-detect available dimensions
        self._available_dimensions = self._detect_dimensions()
        self._default_dimension = self._select_default_dimension()

        # Detect optional per-cell vector fields aligned to UMAP.
        self._vector_fields_metadata = self._detect_vector_fields()

        # Auto-detect latent space
        if self.latent_key is None:
            self.latent_key = self._detect_latent_key()

        logger.info(f"AnnDataAdapter initialized: {self.n_cells:,} cells, "
                   f"{self.n_genes:,} genes, dimensions: {self._available_dimensions}")

    @classmethod
    def from_file(
        cls,
        path: Union[str, Path],
        backed: Union[bool, Literal["r", "r+"]] = "r",
        **kwargs
    ) -> "AnnDataAdapter":
        """
        Create adapter from h5ad file or zarr store with lazy loading.

        Supports both:
        - .h5ad files: HDF5-based, supports true backed mode with memory-mapping
        - .zarr directories: Directory-based, arrays are loaded on-demand

        Lazy Loading Behavior
        ---------------------

        h5ad (backed mode):
            When backed='r', the file is memory-mapped. Only accessed data is
            loaded into RAM. This is ideal for large datasets. The X matrix and
            layer arrays support lazy column/row access.

        zarr:
            Zarr stores individual arrays as separate files on disk. While
            anndata.read_zarr() loads the AnnData structure (obs, var metadata),
            the actual X matrix data is loaded lazily when accessed. This is
            because zarr's internal chunking mechanism defers loading until
            data is requested. Note: zarr does not support the same backed mode
            API as h5ad, but achieves similar lazy behavior through its design.

        Parameters
        ----------
        path : str or Path
            Path to h5ad file or zarr directory.
            - For h5ad: path/to/file.h5ad
            - For zarr: path/to/store.zarr (must be a directory)
        backed : bool or 'r' or 'r+'
            For h5ad only:
            - 'r': Read-only backed mode (recommended for visualization)
            - 'r+': Read-write backed mode
            - True: Same as 'r'
            - False: Load entire file into memory
            For zarr: This parameter is ignored (zarr is always lazy).
        **kwargs
            Additional arguments passed to AnnDataAdapter.__init__:
            - latent_key: Key in obsm for latent space
            - gene_id_column: Column in var for gene IDs
            - normalize_embeddings: Normalize UMAP to [-1,1]
            - dataset_name: Human-readable name

        Returns
        -------
        AnnDataAdapter
            Adapter instance wrapping the loaded data.

        Raises
        ------
        FileNotFoundError
            If the path does not exist.
        ValueError
            If the path is not a valid h5ad or zarr store.

        Examples
        --------
        >>> # Load h5ad with lazy loading (default)
        >>> adapter = AnnDataAdapter.from_file("data.h5ad")
        >>> # Load h5ad fully into memory
        >>> adapter = AnnDataAdapter.from_file("data.h5ad", backed=False)
        >>> # Load zarr store
        >>> adapter = AnnDataAdapter.from_file("data.zarr")
        """
        import anndata as ad

        path = Path(path).resolve()

        if not path.exists():
            raise FileNotFoundError(f"Path not found: {path}")

        # Detect format based on path structure
        # Zarr stores are directories containing .zgroup or .zattrs files
        is_zarr = (
            path.suffix == ".zarr"
            or (path.is_dir() and (path / ".zgroup").exists())
            or (path.is_dir() and (path / ".zattrs").exists())
        )

        # Additional check: if it's a directory without h5ad extension, likely zarr
        if path.is_dir() and not is_zarr:
            # Check for zarr indicators
            if any((path / f).exists() for f in [".zgroup", ".zattrs", ".zarray"]):
                is_zarr = True
            elif path.suffix != ".h5ad":
                # If directory doesn't look like zarr and isn't h5ad, warn
                logger.warning(
                    f"Directory '{path}' does not appear to be a zarr store. "
                    f"Expected .zgroup or .zattrs file."
                )

        if is_zarr:
            # Zarr store - data arrays are loaded lazily by zarr's chunking
            if not path.is_dir():
                raise ValueError(
                    f"Zarr path must be a directory: {path}. "
                    f"Zarr stores are directories, not single files."
                )
            adata = ad.read_zarr(path)
            source_type = "zarr"
            is_lazy = True  # Zarr provides lazy array access
            logger.info(f"Loaded zarr store: {path} (lazy array access)")
        else:
            # h5ad file
            if not path.is_file():
                raise ValueError(f"h5ad path must be a file: {path}")
            if backed:
                backed_mode = backed if isinstance(backed, str) else "r"
                adata = ad.read_h5ad(path, backed=backed_mode)
                is_lazy = True
            else:
                adata = ad.read_h5ad(path)
                is_lazy = False
            source_type = "h5ad"
            logger.info(f"Loaded h5ad file: {path} (backed={backed}, lazy={is_lazy})")

        # Use filename/dirname as default name
        if "dataset_name" not in kwargs:
            kwargs["dataset_name"] = path.stem if path.is_file() else path.name

        # Create adapter and store source info
        adapter = cls(adata, **kwargs)
        adapter._source_path = str(path)
        adapter._source_type = source_type
        adapter._is_lazy = is_lazy
        return adapter

    def _check_closed(self) -> None:
        """Raise error if adapter has been closed."""
        if self._closed:
            raise RuntimeError(
                "AnnDataAdapter has been closed. Create a new adapter instance."
            )

    @property
    def n_cells(self) -> int:
        """Number of cells."""
        self._check_closed()
        return self.adata.n_obs

    @property
    def n_genes(self) -> int:
        """Number of genes."""
        self._check_closed()
        return self.adata.n_vars

    @property
    def is_backed(self) -> bool:
        """
        Whether the AnnData is backed (lazy loading from disk).

        Returns False if the adapter is closed or if adata is None.
        """
        if self._closed or self.adata is None:
            return False
        try:
            return self.adata.isbacked
        except Exception:
            return False

    # =========================================================================
    # DIMENSION DETECTION
    # =========================================================================

    def _detect_dimensions(self) -> list[int]:
        """
        Detect available UMAP dimensions from obsm keys.

        Looks for embeddings in this order:
        1. Explicit dimension keys: X_umap_1d, X_umap_2d, X_umap_3d
        2. Generic key: X_umap (interpreted as 1D/2D/3D if shape matches and no clash)

        Returns
        -------
        list[int]
            Sorted list of available dimensions (e.g., [2, 3])

        Raises
        ------
        ValueError
            If no valid UMAP embeddings are found in obsm.
        """
        # Edge case: empty AnnData
        if self.n_cells == 0:
            logger.warning("AnnData has 0 cells. Visualization may not work correctly.")

        # Reset resolution metadata on each detection pass
        self._umap_embedding_key_by_dim = {}
        self._umap_resolution = None

        available: list[int] = []

        def _shape_of(arr: Any) -> tuple[int, ...]:
            if sparse.issparse(arr):
                return tuple(arr.shape)
            try:
                shape = arr.shape  # type: ignore[attr-defined]
                return tuple(shape)
            except Exception:
                return tuple(np.asarray(arr).shape)

        # Check for explicit dimension keys first
        for dim in [1, 2, 3]:
            key = f"X_umap_{dim}d"
            if key in self.adata.obsm:
                arr = self.adata.obsm[key]
                # Handle sparse obsm (rare but possible)
                arr_shape = _shape_of(arr)

                if len(arr_shape) != 2:
                    logger.warning(f"Key '{key}' is not 2D (shape: {arr_shape}), skipping")
                    continue
                if arr_shape[0] != self.n_cells:
                    logger.warning(
                        f"Key '{key}' has {arr_shape[0]} rows but AnnData has {self.n_cells} cells"
                    )
                    continue
                if arr_shape[1] == dim:
                    available.append(dim)
                    self._umap_embedding_key_by_dim[dim] = key
                else:
                    logger.warning(
                        f"Key '{key}' has {arr_shape[1]} columns, expected {dim}"
                    )

        # Optional fallback: interpret X_umap as a dimensioned embedding when unambiguous.
        if "X_umap" in self.adata.obsm:
            arr = self.adata.obsm["X_umap"]
            arr_shape = _shape_of(arr)

            if len(arr_shape) != 2:
                self._umap_resolution = {
                    "source_key": "X_umap",
                    "action": "ignored",
                    "reason": "invalid_shape",
                    "shape": list(arr_shape),
                }
            elif arr_shape[0] != self.n_cells:
                self._umap_resolution = {
                    "source_key": "X_umap",
                    "action": "ignored",
                    "reason": "row_count_mismatch",
                    "shape": list(arr_shape),
                    "expected_rows": self.n_cells,
                }
            else:
                n_dims = int(arr_shape[1])
                alias_key = f"X_umap_{n_dims}d"
                if n_dims not in (1, 2, 3):
                    self._umap_resolution = {
                        "source_key": "X_umap",
                        "action": "ignored",
                        "reason": "unsupported_dim",
                        "n_dims": n_dims,
                    }
                elif n_dims in self._umap_embedding_key_by_dim:
                    # Clash: explicit key for this dimension exists; ignore X_umap.
                    self._umap_resolution = {
                        "source_key": "X_umap",
                        "action": "ignored",
                        "reason": "explicit_key_present",
                        "dim": n_dims,
                        "alias_key": alias_key,
                    }
                else:
                    available.append(n_dims)
                    self._umap_embedding_key_by_dim[n_dims] = "X_umap"
                    self._umap_resolution = {
                        "source_key": "X_umap",
                        "action": "used_as",
                        "dim": n_dims,
                        "alias_key": alias_key,
                    }

        if not available:
            raise ValueError(
                "No valid UMAP embeddings found in adata.obsm. "
                "Expected keys: 'X_umap_1d', 'X_umap_2d', 'X_umap_3d', or 'X_umap' (1D/2D/3D). "
                f"Available obsm keys: {list(self.adata.obsm.keys())}. "
                "Please add UMAP coordinates to your AnnData object using scanpy.tl.umap() "
                "and then store explicit embeddings under X_umap_1d / X_umap_2d / X_umap_3d."
            )

        # Ensure uniqueness and stable ordering.
        return sorted(set(available))

    def _select_default_dimension(self) -> int:
        """Select default dimension (prefer 3D > 2D > 1D)."""
        for dim in [3, 2, 1]:
            if dim in self._available_dimensions:
                return dim
        return self._available_dimensions[0]

    def _detect_vector_fields(self) -> Optional[dict[str, Any]]:
        """
        Detect per-cell vector fields in `adata.obsm` aligned to UMAP.

        Naming convention (UMAP basis):
        - Explicit: `<field>_umap_<dim>d` (preferred)
          - Examples: `velocity_umap_2d`, `T_fwd_umap_3d`
        - Implicit: `<field>_umap` with shape `(n_cells, 1|2|3)`
          - Used only if the explicit key for that dim is not present (clash-safe)

        Returns a structure suitable for dataset_identity.json under `vector_fields`,
        including per-field `files` entries pointing at server-served binaries:
        `vectors/<fieldId>_<dim>d.bin`.
        """
        keys = list(self.adata.obsm.keys())
        if not keys:
            return None

        def _shape_of(arr: Any) -> tuple[int, ...]:
            if sparse.issparse(arr):
                return tuple(arr.shape)
            try:
                shape = arr.shape  # type: ignore[attr-defined]
                return tuple(shape)
            except Exception:
                return tuple(np.asarray(arr).shape)

        suffix_re = re.compile(r"^(.*)_([123])d$")

        # Internal working structure: fieldId -> {dims:set, explicit:set, obsm_keys:dict[str,str]}
        fields: dict[str, dict[str, Any]] = {}

        def _ensure(field_id: str) -> dict[str, Any]:
            if field_id not in fields:
                fields[field_id] = {"dims": set(), "explicit": set(), "obsm_keys": {}}
            return fields[field_id]

        # 1) Explicit per-dimension keys.
        for key in keys:
            match = suffix_re.match(key)
            if not match:
                continue
            field_id = match.group(1)
            dim = int(match.group(2))

            if not field_id.endswith("_umap"):
                continue
            if field_id.startswith("X_"):
                continue  # reserved for embeddings
            if dim not in (1, 2, 3):
                continue
            if dim not in self._available_dimensions:
                continue  # can't render a dim that doesn't exist in points

            arr = self.adata.obsm.get(key)
            if arr is None:
                continue
            shape = _shape_of(arr)
            if len(shape) != 2 or shape[0] != self.n_cells or shape[1] != dim:
                continue

            entry = _ensure(field_id)
            entry["obsm_keys"][f"{dim}d"] = key
            entry["dims"].add(dim)
            entry["explicit"].add(dim)

        # 2) Implicit keys (<field>_umap), shape-derived, clash-safe.
        for key in keys:
            field_id = key
            if not field_id.endswith("_umap"):
                continue
            if field_id.startswith("X_"):
                continue
            if suffix_re.match(field_id):
                continue  # skip explicit keys

            arr = self.adata.obsm.get(field_id)
            if arr is None:
                continue
            shape = _shape_of(arr)
            if len(shape) == 1:
                if shape[0] != self.n_cells:
                    continue
                dim = 1
            elif len(shape) == 2:
                if shape[0] != self.n_cells:
                    continue
                dim = int(shape[1])
            else:
                continue

            if dim not in (1, 2, 3):
                continue
            if dim not in self._available_dimensions:
                continue

            entry = _ensure(field_id)
            if dim in entry["explicit"]:
                continue  # clash-safe: explicit wins
            slot = f"{dim}d"
            if slot not in entry["obsm_keys"]:
                entry["obsm_keys"][slot] = field_id
                entry["dims"].add(dim)

        if not fields:
            return None

        def _label_for(field_id: str) -> str:
            base = field_id[:-5] if field_id.endswith("_umap") else field_id
            base = base.replace("_", " ").strip()
            titled = (base[:1].upper() + base[1:]) if base else field_id
            return f"{titled} (UMAP)" if field_id.endswith("_umap") else titled

        fields_out: dict[str, Any] = {}
        for field_id, entry in fields.items():
            dims = sorted(entry["dims"])
            if not dims:
                continue
            files = {f"{d}d": f"vectors/{field_id}_{d}d.bin" for d in dims}
            fields_out[field_id] = {
                "label": _label_for(field_id),
                "basis": "umap",
                "available_dimensions": dims,
                "default_dimension": max(dims),
                "files": files,
                "obsm_keys": entry["obsm_keys"],
            }

        if not fields_out:
            return None

        default_field = "velocity_umap" if "velocity_umap" in fields_out else sorted(fields_out.keys())[0]
        return {"default_field": default_field, "fields": fields_out}

    def _detect_latent_key(self) -> Optional[str]:
        """Auto-detect latent space key from obsm."""
        # Priority order for latent space detection
        candidates = ["scanvi", "scvi", "X_scvi", "X_pca", "harmony", "X_harmony"]

        for key in candidates:
            if key in self.adata.obsm:
                logger.info(f"Auto-detected latent space: '{key}'")
                return key

        # Fall back to first available obsm that's not a UMAP
        for key in self.adata.obsm.keys():
            if not key.startswith("X_umap"):
                logger.info(f"Using '{key}' as latent space (first non-UMAP obsm)")
                return key

        logger.warning("No latent space found in obsm, centroid outlier computation will be skipped")
        return None

    # =========================================================================
    # EMBEDDING DATA
    # =========================================================================

    def get_embedding(self, dim: int) -> np.ndarray:
        """
        Get embedding coordinates for a dimension.

        Returns normalized Float32 array of shape (n_cells, dim).
        """
        self._check_closed()

        if dim not in self._available_dimensions:
            raise ValueError(f"Dimension {dim}D not available. Available: {self._available_dimensions}")

        if dim in self._embedding_cache:
            return self._embedding_cache[dim]

        # Get raw embedding
        key = self._umap_embedding_key_by_dim.get(dim)
        if not key or key not in self.adata.obsm:
            raise ValueError(
                f"Could not find embedding for dimension {dim}. "
                f"Expected an embedding resolved for {dim}D. Available obsm keys: {list(self.adata.obsm.keys())}"
            )
        raw = _to_dense_2d(self.adata.obsm[key])

        raw = raw.astype(np.float32)

        if self.normalize_embeddings:
            # Normalize to [-1, 1] range (same as prepare)
            axis_mins = raw.min(axis=0)
            axis_maxs = raw.max(axis=0)
            axis_ranges = axis_maxs - axis_mins
            max_range = float(axis_ranges.max())
            if max_range < 1e-8:
                max_range = 1.0
            center = (axis_mins + axis_maxs) / 2
            scale_factor = 2.0 / max_range
            raw = ((raw - center) * scale_factor).astype(np.float32)
            self._embedding_norm[dim] = {
                "center": center.astype(np.float32),
                "scale_factor": float(scale_factor),
            }
        else:
            self._embedding_norm[dim] = {
                "center": np.zeros((dim,), dtype=np.float32),
                "scale_factor": 1.0,
            }

        self._embedding_cache[dim] = raw
        return raw

    def get_embedding_3d(self, dim: int) -> np.ndarray:
        """
        Get embedding padded to 3D for WebGL rendering.

        1D -> (x, 0, 0)
        2D -> (x, y, 0)
        3D -> (x, y, z)
        """
        embedding = self.get_embedding(dim)
        n_cells = embedding.shape[0]

        if dim == 3:
            return embedding
        elif dim == 2:
            result = np.zeros((n_cells, 3), dtype=np.float32)
            result[:, :2] = embedding
            return result
        elif dim == 1:
            result = np.zeros((n_cells, 3), dtype=np.float32)
            result[:, 0] = embedding[:, 0]
            return result
        else:
            raise ValueError(f"Unsupported dimension: {dim}")

    def get_points_binary(self, dim: int, compress: bool = False) -> bytes:
        """Get embedding as binary data (for HTTP response)."""
        embedding = self.get_embedding(dim)
        data = embedding.astype(np.float32).tobytes()

        if compress:
            return gzip.compress(data, compresslevel=6)
        return data

    def get_vector_field_binary(self, field_id: str, dim: int, compress: bool = False) -> bytes:
        """
        Get a per-cell vector field (displacement vectors) as binary float32 data.

        Vector fields are scaled by the SAME per-dimension normalization scale as
        the embedding points, so they are in the same normalized space as the
        `points_{dim}d.bin` responses.
        """
        self._check_closed()

        field = str(field_id or "")
        d = int(dim)
        if d not in (1, 2, 3):
            raise ValueError(f"Invalid dimension {dim}. Expected 1, 2, or 3.")

        meta = self._vector_fields_metadata or {}
        fields = meta.get("fields") or {}
        entry = fields.get(field)
        if not entry:
            raise ValueError(f"Vector field '{field}' not available")

        obsm_keys = entry.get("obsm_keys") or {}
        obsm_key = obsm_keys.get(f"{d}d")
        if not obsm_key or obsm_key not in self.adata.obsm:
            raise ValueError(f"Vector field '{field}' does not provide {d}D data")

        cache_key = (field, d)
        if cache_key not in self._vector_field_cache:
            # Ensure embedding normalization scale is computed.
            self.get_embedding(d)
            scale_factor = float(self._embedding_norm.get(d, {}).get("scale_factor", 1.0))

            raw = _to_dense_2d(self.adata.obsm[obsm_key]).astype(np.float32)
            if raw.ndim == 1:
                raw = raw.reshape(-1, 1)
            if raw.ndim != 2 or raw.shape[0] != self.n_cells or raw.shape[1] != d:
                raise ValueError(
                    f"Vector field '{obsm_key}' has shape {tuple(raw.shape)}, expected ({self.n_cells}, {d})"
                )

            if self.normalize_embeddings and scale_factor != 1.0:
                raw *= scale_factor

            self._vector_field_cache[cache_key] = raw

        data = self._vector_field_cache[cache_key].tobytes()
        if compress:
            return gzip.compress(data, compresslevel=6)
        return data

    # =========================================================================
    # LATENT SPACE FOR OUTLIER COMPUTATION
    # =========================================================================

    def _get_latent_space(self) -> Optional[np.ndarray]:
        """Get latent space for outlier quantile computation."""
        if self._latent_space is not None:
            return self._latent_space

        if self.latent_key is None or self.latent_key not in self.adata.obsm:
            return None

        self._latent_space = _to_dense_2d(self.adata.obsm[self.latent_key]).astype(np.float32)
        return self._latent_space

    # =========================================================================
    # OBS DATA (CELL METADATA)
    # =========================================================================

    def get_obs_keys(self) -> list[str]:
        """Get list of obs column names."""
        self._check_closed()
        return list(self.adata.obs.columns)

    def get_obs_field_kind(self, key: str) -> Literal["continuous", "category"]:
        """
        Determine if an obs field is continuous or categorical.

        Classification rules:
        - Categorical dtype → category
        - Boolean dtype → category
        - Numeric dtype → continuous
        - String/object → category (treated as labels)
        - Empty column → category (safe default)
        """
        self._check_closed()

        if key not in self.adata.obs.columns:
            raise KeyError(f"Field '{key}' not found in obs. Available: {list(self.adata.obs.columns)}")

        s = self.adata.obs[key]

        # Handle empty series
        if len(s) == 0:
            logger.warning(f"obs field '{key}' is empty, treating as categorical")
            return "category"

        # Check categorical dtype (pandas 2.0+ compatible)
        try:
            is_cat = isinstance(s.dtype, pd.CategoricalDtype)
        except AttributeError:
            # Fallback for older pandas versions
            is_cat = hasattr(s, 'cat')

        if is_cat:
            return "category"
        elif pd.api.types.is_bool_dtype(s):
            return "category"
        elif pd.api.types.is_numeric_dtype(s):
            # Check if all values are NaN - still treat as continuous
            return "continuous"
        elif pd.api.types.is_string_dtype(s) or pd.api.types.is_object_dtype(s):
            return "category"
        else:
            # Unknown dtype - safe default is categorical
            logger.warning(f"obs field '{key}' has unknown dtype {s.dtype}, treating as categorical")
            return "category"

    def get_obs_continuous_values(self, key: str, compress: bool = False) -> bytes:
        """
        Get continuous obs field as binary float32 data.

        NaN/Inf values are preserved in the output (client handles visualization).

        Raises
        ------
        KeyError
            If the field is not found in obs.
        """
        self._check_closed()

        if key not in self.adata.obs.columns:
            raise KeyError(
                f"obs field '{key}' not found. Available fields: {list(self.adata.obs.columns)}"
            )

        s = self.adata.obs[key]
        values = pd.to_numeric(s, errors="coerce").to_numpy(dtype=np.float32)

        # Log warning for problematic data
        n_nan = np.isnan(values).sum()
        n_inf = np.isinf(values).sum()
        if n_nan > 0 or n_inf > 0:
            logger.debug(f"obs field '{key}': {n_nan} NaN, {n_inf} Inf values")

        data = values.tobytes()

        if compress:
            return gzip.compress(data, compresslevel=6)
        return data

    def get_obs_categorical_codes(self, key: str, compress: bool = False) -> tuple[bytes, list[str], int]:
        """
        Get categorical obs field as binary codes.

        Returns:
            (binary_codes, category_list, missing_value)

        Categories are assigned codes 0 to n-1.
        Missing values (NaN) are encoded as the missing_value sentinel.

        Raises
        ------
        KeyError
            If the field is not found in obs.
        """
        self._check_closed()

        if key not in self.adata.obs.columns:
            raise KeyError(
                f"obs field '{key}' not found. Available fields: {list(self.adata.obs.columns)}"
            )

        s = self.adata.obs[key]
        cat = s.astype("category")
        categories = [str(c) for c in cat.cat.categories]
        codes = cat.cat.codes.to_numpy(dtype=np.int32)

        # Handle edge case: empty categories
        n_categories = len(categories)
        if n_categories == 0:
            logger.warning(f"obs field '{key}' has no categories (all values are NaN/missing)")
            # Return all missing values
            dtype = np.uint8
            missing_value = 255
            codes_typed = np.full(self.n_cells, missing_value, dtype=dtype)
            data = codes_typed.tobytes()
            if compress:
                data = gzip.compress(data, compresslevel=6)
            return data, [], int(missing_value)

        # Select dtype based on number of categories
        if n_categories <= 254:
            dtype = np.uint8
            missing_value = 255
        else:
            dtype = np.uint16
            missing_value = 65535

        # Convert codes (-1 for NaN -> missing_value)
        codes_typed = np.full(self.n_cells, missing_value, dtype=dtype)
        valid_mask = codes >= 0
        codes_typed[valid_mask] = codes[valid_mask].astype(dtype)

        # Log missing value count
        n_missing = (~valid_mask).sum()
        if n_missing > 0:
            logger.debug(f"obs field '{key}': {n_missing} missing values out of {self.n_cells}")

        data = codes_typed.tobytes()
        if compress:
            data = gzip.compress(data, compresslevel=6)

        return data, categories, int(missing_value)

    # =========================================================================
    # CENTROIDS AND OUTLIER QUANTILES
    # =========================================================================

    def _compute_centroids_for_field(
        self,
        key: str,
        dim: int,
    ) -> list[dict]:
        """Compute centroids for a categorical field at a given dimension."""
        cache_key = f"{key}_{dim}d"
        if cache_key in self._centroid_cache:
            return self._centroid_cache[cache_key]

        coords = self.get_embedding(dim)
        s = self.adata.obs[key]
        cat = s.astype("category")
        categories = [str(c) for c in cat.cat.categories]
        codes = cat.cat.codes.to_numpy(dtype=np.int32)

        centroids = []
        for code, label in enumerate(categories):
            mask = codes == code
            idx = np.nonzero(mask)[0]
            n = idx.size

            if n < self.centroid_min_points:
                continue

            pts = coords[idx, :]
            center = pts.mean(axis=0)

            # Remove outliers for centroid computation
            if n > self.centroid_min_points:
                dists = np.linalg.norm(pts - center, axis=1)
                thr = float(np.quantile(dists, self.centroid_outlier_quantile))
                inlier_mask = dists <= thr
                n_in = int(inlier_mask.sum())
                if n_in >= self.centroid_min_points:
                    pts_in = pts[inlier_mask, :]
                    center = pts_in.mean(axis=0)
                    used_count = n_in
                else:
                    used_count = n
            else:
                used_count = n

            centroids.append({
                "category": str(label),
                "position": center.astype(float).tolist(),
                "n_points": int(used_count),
            })

        self._centroid_cache[cache_key] = centroids
        return centroids

    def get_centroids_for_field(self, key: str) -> dict[str, list[dict]]:
        """Get centroids for all available dimensions."""
        result = {}
        for dim in self._available_dimensions:
            centroids = self._compute_centroids_for_field(key, dim)
            result[str(dim)] = centroids
        return result

    def _compute_outlier_quantiles(self, key: str) -> np.ndarray:
        """Compute per-cell outlier quantiles for a categorical field."""
        if key in self._outlier_quantile_cache:
            return self._outlier_quantile_cache[key]

        latent = self._get_latent_space()
        if latent is None:
            # Return all zeros if no latent space
            result = np.zeros(self.n_cells, dtype=np.float32)
            self._outlier_quantile_cache[key] = result
            return result

        s = self.adata.obs[key]
        cat = s.astype("category")
        categories = [str(c) for c in cat.cat.categories]
        codes = cat.cat.codes.to_numpy(dtype=np.int32)

        quantiles = np.full(self.n_cells, np.nan, dtype=np.float32)

        for code, label in enumerate(categories):
            mask = codes == code
            idx = np.nonzero(mask)[0]
            n = idx.size

            if n < self.centroid_min_points:
                continue

            pts = latent[idx, :]
            centroid = pts.mean(axis=0)
            dists = np.linalg.norm(pts - centroid, axis=1)
            sorted_dists = np.sort(dists)
            ranks = np.searchsorted(sorted_dists, dists, side='right')
            cell_quantiles = ranks.astype(np.float32) / n
            quantiles[idx] = cell_quantiles

        self._outlier_quantile_cache[key] = quantiles
        return quantiles

    def get_obs_outlier_quantiles(self, key: str, compress: bool = False) -> bytes:
        """Get outlier quantiles as binary float32 data."""
        quantiles = self._compute_outlier_quantiles(key)
        data = quantiles.tobytes()

        if compress:
            return gzip.compress(data, compresslevel=6)
        return data

    # =========================================================================
    # VAR DATA (GENE EXPRESSION)
    # =========================================================================

    def get_gene_ids(self) -> list[str]:
        """Get list of gene identifiers."""
        self._check_closed()

        if self._gene_ids_cache is not None:
            return self._gene_ids_cache

        if self.gene_id_column == "index" or self.gene_id_column is None:
            self._gene_ids_cache = self.adata.var.index.astype(str).tolist()
        else:
            if self.gene_id_column not in self.adata.var.columns:
                raise KeyError(
                    f"gene_id_column '{self.gene_id_column}' not found in var. "
                    f"Available columns: {list(self.adata.var.columns)}"
                )
            self._gene_ids_cache = self.adata.var[self.gene_id_column].astype(str).tolist()

        # Build index for O(1) lookup
        self._gene_id_to_idx_cache = {gid: idx for idx, gid in enumerate(self._gene_ids_cache)}
        return self._gene_ids_cache

    def _get_gene_idx(self, gene_id: str) -> int:
        """Get gene index with O(1) lookup."""
        # Ensure cache is populated
        if self._gene_id_to_idx_cache is None:
            self.get_gene_ids()

        if gene_id not in self._gene_id_to_idx_cache:
            raise KeyError(f"Gene '{gene_id}' not found in var")

        return self._gene_id_to_idx_cache[gene_id]

    def _get_gene_column(self, gene_idx: int) -> np.ndarray:
        """
        Extract a single gene column from the expression matrix.

        Handles multiple data formats efficiently:
        - Dense numpy arrays: Direct column access O(1)
        - CSC sparse matrices: Efficient column access O(nnz/n_cols)
        - CSR sparse matrices: Convert to CSC and cache for repeated access
        - Other sparse formats (COO, DOK, etc.): Convert to CSC and cache
        - None: Returns zeros (edge case for AnnData without X matrix)
        - Backed/chunked arrays: Direct slicing (no caching to preserve lazy loading)

        Memory-Speed Tradeoff:
        For in-memory CSR matrices, the CSC cache doubles memory usage but enables
        O(nnz/n_cols) column access instead of O(nnz). This pays off after
        accessing more than a few genes.

        IMPORTANT: For backed h5ad files, we do NOT create a CSC cache because:
        1. It would force loading the entire X matrix into memory
        2. It defeats the purpose of lazy loading for large datasets
        Instead, we accept slower column access to preserve memory efficiency.

        Parameters
        ----------
        gene_idx : int
            Column index of the gene to extract.

        Returns
        -------
        np.ndarray
            1D array of expression values for all cells.
        """
        self._check_closed()

        X = self.adata.X

        # Edge case: no expression matrix
        if X is None:
            logger.warning("AnnData has no X matrix. Returning zeros for gene expression.")
            return np.zeros(self.n_cells, dtype=np.float32)

        # Edge case: invalid gene index
        if gene_idx < 0 or gene_idx >= self.n_genes:
            raise IndexError(f"Gene index {gene_idx} out of range [0, {self.n_genes})")

        # CRITICAL: For backed h5ad files, avoid creating CSC cache to preserve
        # lazy loading behavior. The memory cost of CSC cache would defeat the
        # purpose of backed mode for large datasets.
        if self._is_lazy and self._source_type == "h5ad":
            # For backed h5ad, use direct slicing (slower but memory-efficient)
            # Note: backed h5ad matrices support slicing but may be slower
            col = X[:, gene_idx]
            if sparse.issparse(col):
                return np.asarray(col.toarray()).flatten()
            return np.asarray(col).flatten()

        if sparse.issparse(X):
            # For sparse matrix, CSC is efficient for column access
            if sparse.isspmatrix_csc(X):
                return np.asarray(X.getcol(gene_idx).toarray()).flatten()
            elif sparse.isspmatrix_csr(X):
                # CSR is inefficient for column access - convert to CSC and cache
                # This is a one-time cost that pays off with repeated gene queries
                # Only do this for in-memory data (not backed files)
                if self._X_csc_cache is None:
                    logger.info(
                        f"Converting CSR matrix ({X.shape[0]:,}×{X.shape[1]:,}, "
                        f"{X.nnz:,} non-zeros) to CSC for efficient column access"
                    )
                    self._X_csc_cache = X.tocsc()
                return np.asarray(self._X_csc_cache.getcol(gene_idx).toarray()).flatten()
            else:
                # Other sparse format (COO, LIL, DOK, BSR) - convert to CSC
                if self._X_csc_cache is None:
                    logger.info(f"Converting {type(X).__name__} to CSC for column access")
                    self._X_csc_cache = sparse.csc_matrix(X)
                return np.asarray(self._X_csc_cache.getcol(gene_idx).toarray()).flatten()
        else:
            # Dense matrix - direct column access
            return np.asarray(X[:, gene_idx]).flatten()

    def get_gene_expression(self, gene_id: str, compress: bool = False) -> bytes:
        """Get expression values for a single gene as binary float32."""
        # Check LRU cache first (O(1) operations)
        values = self._gene_expression_cache.get(gene_id)
        if values is None:
            # Fetch from data
            gene_idx = self._get_gene_idx(gene_id)
            col = self._get_gene_column(gene_idx)
            values = col.astype(np.float32)

            # Add to LRU cache (handles eviction automatically)
            self._gene_expression_cache.put(gene_id, values)

        data = values.tobytes()
        if compress:
            return gzip.compress(data, compresslevel=6)
        return data

    def get_gene_min_max(self, gene_id: str) -> tuple[float, float]:
        """Get min/max values for a gene (for colormap scaling)."""
        gene_idx = self._get_gene_idx(gene_id)
        col = self._get_gene_column(gene_idx)
        valid = col[np.isfinite(col)]
        if len(valid) == 0:
            return 0.0, 1.0
        return float(valid.min()), float(valid.max())

    # =========================================================================
    # CONNECTIVITY DATA
    # =========================================================================

    def has_connectivity(self) -> bool:
        """Check if connectivity data is available."""
        if self._closed:
            return False
        return "connectivities" in self.adata.obsp

    def _compute_connectivity_edges(self) -> tuple[np.ndarray, np.ndarray, int, int]:
        """
        Compute connectivity edges (cached).

        Returns:
            (sources_array, destinations_array, n_edges, max_neighbors)
        """
        self._check_closed()

        if self._connectivity_cache is not None:
            return self._connectivity_cache

        if not self.has_connectivity():
            raise ValueError("No connectivity data in adata.obsp['connectivities']")

        connectivities = self.adata.obsp["connectivities"]

        if not sparse.isspmatrix_csr(connectivities):
            connectivities = sparse.csr_matrix(connectivities)

        # Symmetrize and binarize (make a copy to avoid modifying original)
        connectivities_sym = (connectivities + connectivities.T).tocsr()
        connectivities_sym.data = np.ones_like(connectivities_sym.data)

        indptr = connectivities_sym.indptr
        indices = connectivities_sym.indices

        # Pre-allocate arrays for edges (upper bound: nnz/2)
        max_edges = connectivities_sym.nnz // 2 + 1
        edge_sources = np.empty(max_edges, dtype=np.int64)
        edge_destinations = np.empty(max_edges, dtype=np.int64)
        edge_count = 0
        max_neighbors = 0

        for cell_idx in range(self.n_cells):
            start = indptr[cell_idx]
            end = indptr[cell_idx + 1]
            neighbor_count = end - start
            if neighbor_count > max_neighbors:
                max_neighbors = neighbor_count

            for j in range(start, end):
                neighbor_idx = indices[j]
                if cell_idx < neighbor_idx:
                    edge_sources[edge_count] = cell_idx
                    edge_destinations[edge_count] = neighbor_idx
                    edge_count += 1

        # Trim to actual size
        edge_sources = edge_sources[:edge_count]
        edge_destinations = edge_destinations[:edge_count]

        # Select dtype based on cell count
        if self.n_cells <= 65535:
            dtype = np.uint16
        elif self.n_cells <= 4294967295:
            dtype = np.uint32
        else:
            dtype = np.uint64

        sources = edge_sources.astype(dtype)
        destinations = edge_destinations.astype(dtype)

        # Sort for better compression
        sort_idx = np.lexsort((destinations, sources))
        sources = sources[sort_idx]
        destinations = destinations[sort_idx]

        self._connectivity_cache = (sources, destinations, edge_count, int(max_neighbors))
        return self._connectivity_cache

    def get_connectivity_edges(self, compress: bool = False) -> tuple[bytes, bytes, int, int]:
        """
        Get connectivity edges as binary data.

        Returns:
            (sources_binary, destinations_binary, n_edges, max_neighbors)
        """
        sources, destinations, n_edges, max_neighbors = self._compute_connectivity_edges()

        sources_data = sources.tobytes()
        destinations_data = destinations.tobytes()

        if compress:
            sources_data = gzip.compress(sources_data, compresslevel=6)
            destinations_data = gzip.compress(destinations_data, compresslevel=6)

        return sources_data, destinations_data, n_edges, max_neighbors

    # =========================================================================
    # MANIFEST GENERATION
    # =========================================================================

    def get_dataset_identity(self) -> dict:
        """Generate dataset_identity.json content."""
        # Count obs fields
        obs_keys = self.get_obs_keys()
        n_categorical = sum(1 for k in obs_keys if self.get_obs_field_kind(k) == "category")
        n_continuous = sum(1 for k in obs_keys if self.get_obs_field_kind(k) == "continuous")

        # Check connectivity
        has_conn = self.has_connectivity()
        n_edges = None
        if has_conn:
            try:
                _, _, n_edges, _ = self.get_connectivity_edges(compress=False)
            except Exception:
                pass

        # Build embeddings metadata
        embeddings_meta = {
            "available_dimensions": self._available_dimensions,
            "default_dimension": self._default_dimension,
            "files": {f"{dim}d": f"points_{dim}d.bin" for dim in self._available_dimensions}
        }
        if self._umap_resolution is not None:
            embeddings_meta["umap_resolution"] = self._umap_resolution

        # Build obs field summaries
        obs_fields = []
        for key in obs_keys:
            kind = self.get_obs_field_kind(key)
            entry = {"key": key, "kind": kind}
            if kind == "category":
                s = self.adata.obs[key]
                cat = s.astype("category")
                entry["n_categories"] = len(cat.cat.categories)
            obs_fields.append(entry)

        identity = {
            "version": 2,
            "id": self.dataset_id,
            "name": self.dataset_name,
            "description": "Loaded directly from AnnData",
            "created_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "cellucid_data_version": "anndata_adapter",
            "stats": {
                "n_cells": self.n_cells,
                "n_genes": self.n_genes,
                "n_obs_fields": len(obs_keys),
                "n_categorical_fields": n_categorical,
                "n_continuous_fields": n_continuous,
                "has_connectivity": has_conn,
                "n_edges": n_edges,
            },
            "embeddings": embeddings_meta,
            "obs_fields": obs_fields,
            "export_settings": {
                # No file compression (data is served dynamically)
                # Network compression (gzip) is applied transparently by the server
                "compression": None,
                # No quantization - full float32 precision for all values
                # This ensures no data loss but increases network transfer size
                "var_quantization": None,
                "obs_continuous_quantization": None,
                "obs_categorical_dtype": "auto",
            },
            "source": {
                "type": "anndata_adapter",
                "format": self._source_type,  # "h5ad", "zarr", or "memory"
                "path": self._source_path,
                "is_backed": self.is_backed,
            },
            # Special flag to indicate this is from AnnData adapter
            "_anndata_adapter": True,
        }

        if self._vector_fields_metadata:
            identity["vector_fields"] = self._vector_fields_metadata

        return identity

    def get_obs_manifest(self) -> dict:
        """Generate obs_manifest.json content."""
        obs_keys = self.get_obs_keys()

        continuous_fields = []
        categorical_fields = []

        for key in obs_keys:
            safe_key = _safe_filename_component(key)
            kind = self.get_obs_field_kind(key)

            if kind == "continuous":
                s = self.adata.obs[key]
                values = pd.to_numeric(s, errors="coerce").to_numpy(dtype=np.float32)
                valid = values[np.isfinite(values)]
                min_val = float(valid.min()) if len(valid) > 0 else 0.0
                max_val = float(valid.max()) if len(valid) > 0 else 1.0
                continuous_fields.append([key, min_val, max_val])
            else:
                # Categorical
                s = self.adata.obs[key]
                cat = s.astype("category")
                categories = [str(c) for c in cat.cat.categories]
                n_categories = len(categories)

                dtype_str = "uint8" if n_categories <= 254 else "uint16"
                missing_value = 255 if n_categories <= 254 else 65535

                # Get centroids for all dimensions
                centroids_by_dim = self.get_centroids_for_field(key)

                # Outlier quantile min/max (always 0-1 since it's a quantile)
                categorical_fields.append([
                    key, categories, dtype_str, missing_value, centroids_by_dim, 0.0, 1.0
                ])

        # Build compact manifest format
        # Path patterns must match prepare format for consistency
        # Note: We don't use compression (.gz) in adapter mode - gzip is applied
        # at the HTTP level when client supports Accept-Encoding: gzip
        obs_schemas = {
            "continuous": {
                "pathPattern": "obs/{key}.values.f32",
                "ext": "f32",
                "dtype": "float32",
                "quantized": False,
            },
            "categorical": {
                "codesPathPattern": "obs/{key}.codes.{ext}",
                "outlierPathPattern": "obs/{key}.outliers.f32",
                "outlierExt": "f32",
                "outlierDtype": "float32",
                "outlierQuantized": False,
            }
        }

        return {
            "_format": "compact_v1",
            "n_points": self.n_cells,
            "centroid_outlier_quantile": self.centroid_outlier_quantile,
            "latent_key": self.latent_key,
            "compression": None,
            "_obsSchemas": obs_schemas,
            "_continuousFields": continuous_fields,
            "_categoricalFields": categorical_fields,
        }

    def get_var_manifest(self) -> dict:
        """Generate var_manifest.json content."""
        gene_ids = self.get_gene_ids()

        # Build field list with min/max values
        # Note: Computing min/max for all genes upfront would be slow
        # We'll compute on-demand and cache, or use a placeholder
        fields = [[gid] for gid in gene_ids]  # No min/max by default

        # Path pattern must match prepare format for consistency
        # Note: We don't use compression (.gz) in adapter mode - gzip is applied
        # at the HTTP level when client supports Accept-Encoding: gzip
        var_schema = {
            "kind": "continuous",
            "pathPattern": "var/{key}.values.f32",
            "ext": "f32",
            "dtype": "float32",
            "quantized": False,
        }

        return {
            "_format": "compact_v1",
            "n_points": self.n_cells,
            "var_gene_id_column": self.gene_id_column,
            "compression": None,
            "quantization": None,
            "_varSchema": var_schema,
            "fields": fields,
        }

    def get_connectivity_manifest(self) -> Optional[dict]:
        """Generate connectivity_manifest.json content."""
        if not self.has_connectivity():
            return None

        try:
            _, _, n_edges, max_neighbors = self.get_connectivity_edges(compress=False)
        except Exception as e:
            logger.error(f"Failed to compute connectivity: {e}")
            return None

        # Determine index dtype
        if self.n_cells <= 65535:
            index_dtype = "uint16"
            index_bytes = 2
        elif self.n_cells <= 4294967295:
            index_dtype = "uint32"
            index_bytes = 4
        else:
            index_dtype = "uint64"
            index_bytes = 8

        return {
            "format": "edge_pairs",
            "n_cells": self.n_cells,
            "n_edges": n_edges,
            "max_neighbors": max_neighbors,
            "index_bytes": index_bytes,
            "index_dtype": index_dtype,
            "sourcesPath": "connectivity/edges.src.bin",
            "destinationsPath": "connectivity/edges.dst.bin",
            "compression": None,
        }

    # =========================================================================
    # CLEANUP AND CONTEXT MANAGER
    # =========================================================================

    def close(self) -> None:
        """
        Close the adapter and release all resources.

        This method:
        1. Clears all caches to free memory (embedding, centroid, CSC, gene expression)
        2. Closes the underlying file handle for backed h5ad files
        3. Marks the adapter as closed to prevent further operations

        Safe to call multiple times. Always call this method when done with the
        adapter, or use the context manager::

            with AnnDataAdapter.from_file("data.h5ad") as adapter:
                # use adapter
            # automatically cleaned up

        Memory Released
        ---------------

        - Embedding cache (normalized UMAP coordinates)
        - Centroid cache (computed label centroids)
        - Outlier quantile cache
        - Gene expression LRU cache (up to 100 gene columns)
        - CSC matrix cache (for CSR->CSC converted matrices)
        - Latent space array
        - Gene ID lookup indices
        """
        if self._closed:
            return

        # Check if backed BEFORE setting _closed flag
        # (is_backed property returns False after _closed is True)
        was_backed = False
        file_to_close = None
        try:
            was_backed = self.adata.isbacked
            if was_backed and hasattr(self.adata, 'file'):
                file_to_close = self.adata.file
        except Exception:
            pass  # adata may already be in invalid state

        self._closed = True

        # Clear all caches to free memory
        self._embedding_cache.clear()
        self._embedding_norm.clear()
        self._vector_field_cache.clear()
        self._vector_fields_metadata = None
        self._centroid_cache.clear()
        self._outlier_quantile_cache.clear()
        self._gene_expression_cache.clear()
        self._latent_space = None
        self._connectivity_cache = None
        self._gene_ids_cache = None
        self._gene_id_to_idx_cache = None

        # Clear CSC cache (can be large for sparse matrices)
        # This is critical for memory management - CSC cache can double memory usage
        if self._X_csc_cache is not None:
            try:
                # For very large sparse matrices, explicit deletion helps GC
                del self._X_csc_cache
            except Exception:
                pass
            self._X_csc_cache = None

        # Close backed file handle if applicable
        if file_to_close is not None:
            try:
                file_to_close.close()
                logger.debug(f"Closed backed file handle for {self.dataset_name}")
            except Exception as e:
                logger.warning(f"Error closing backed AnnData file: {e}")

        # Clear reference to adata to help garbage collection
        # (do this last since we needed it for is_backed check)
        self.adata = None

        logger.debug(f"AnnDataAdapter closed: {self.dataset_name}")

    def __enter__(self) -> "AnnDataAdapter":
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit - ensures cleanup."""
        self.close()

    def __del__(self):
        """
        Cleanup on garbage collection.

        Note: Using __del__ is an anti-pattern but kept for safety.
        Users should prefer:
        - Using context manager: `with AnnDataAdapter(...) as adapter:`
        - Explicitly calling close(): `adapter.close()`
        """
        try:
            self.close()
        except Exception:
            # Ignore errors during GC - object may be partially destroyed
            pass

    def __repr__(self) -> str:
        if self._closed:
            return f"AnnDataAdapter(closed)"
        backed_str = " (backed)" if self.is_backed else ""
        source_str = f" from {self._source_type}" if self._source_type != "memory" else ""
        return (f"AnnDataAdapter({self.n_cells:,} cells, {self.n_genes:,} genes, "
                f"dims={self._available_dimensions}{backed_str}{source_str})")
