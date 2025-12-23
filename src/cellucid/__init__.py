"""
Cellucid - 3D UMAP WebGL visualization toolkit.

Utilities for exporting AnnData/array data for the 3D UMAP web viewer,
serving data locally or remotely, and integrating with Jupyter notebooks.

Main exports:
    - prepare: Export data to cellucid format
    - serve: Serve a dataset directory via HTTP
    - CellucidServer: Server class for more control
    - CellucidViewer: Jupyter notebook integration
    - show_anndata: Directly visualize AnnData (no export needed)
    - serve_anndata: Serve AnnData via HTTP (no export needed)
    - AnnDataAdapter: Adapter for reading AnnData in cellucid format
"""

from importlib import metadata as _metadata

try:
    __version__ = _metadata.version("cellucid")
except _metadata.PackageNotFoundError:  # pragma: no cover - handled in dev installs
    __version__ = "0.0.0"

# All imports are lazy to ensure fast CLI startup.
# Heavy dependencies (numpy, scipy, pandas) are only loaded when needed.
def __getattr__(name):
    # Core export functionality (imports numpy, pandas, scipy)
    if name == "prepare":
        from .prepare_data import prepare
        return prepare

    # Server for pre-exported data (lightweight, no heavy deps)
    if name == "serve":
        from .server import serve
        return serve
    if name == "CellucidServer":
        from .server import CellucidServer
        return CellucidServer

    # Jupyter integration
    if name == "CellucidViewer":
        from .jupyter import CellucidViewer
        return CellucidViewer
    if name == "show":
        from .jupyter import show
        return show

    # AnnData direct visualization
    if name == "show_anndata":
        from .jupyter import show_anndata
        return show_anndata
    if name == "AnnDataViewer":
        from .jupyter import AnnDataViewer
        return AnnDataViewer

    # AnnData server and adapter
    if name == "serve_anndata":
        from .anndata_server import serve_anndata
        return serve_anndata
    if name == "AnnDataServer":
        from .anndata_server import AnnDataServer
        return AnnDataServer
    if name == "AnnDataAdapter":
        from .anndata_adapter import AnnDataAdapter
        return AnnDataAdapter

    # Vector field utilities (CellRank drift, etc.)
    if name == "compute_transition_drift":
        from .vector_fields import compute_transition_drift
        return compute_transition_drift
    if name == "add_transition_drift_to_obsm":
        from .vector_fields import add_transition_drift_to_obsm
        return add_transition_drift_to_obsm

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    # Core export functionality
    "prepare",

    # Serve pre-exported data
    "serve",
    "CellucidServer",

    # Jupyter integration (pre-exported data)
    "CellucidViewer",
    "show",

    # Direct AnnData visualization (no export needed)
    "show_anndata",
    "AnnDataViewer",
    "serve_anndata",
    "AnnDataServer",
    "AnnDataAdapter",

    # Vector field utilities
    "compute_transition_drift",
    "add_transition_drift_to_obsm",

    # Version
    "__version__",
]
