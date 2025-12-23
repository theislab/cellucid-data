# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.0.1a2] - 2025

### Added

- Initial alpha release
- AnnData visualization with UMAP embeddings (1D, 2D, 3D)
- Gene expression overlays with sparse matrix support
- Cell metadata coloring (categorical and continuous)
- Interactive filtering and cell selection
- KNN connectivity visualization
- Multiple deployment modes:
  - Local demo (browser-only)
  - Browser file picker (h5ad and exported formats)
  - Server CLI (`cellucid serve`)
  - Python API (`serve()`, `serve_anndata()`)
  - Jupyter integration (`show()`, `show_anndata()`)
- Export functionality for web deployment

[Unreleased]: https://github.com/theislab/cellucid-python/compare/v0.0.1a2...HEAD
[0.0.1a2]: https://github.com/theislab/cellucid-python/releases/tag/v0.0.1a2
