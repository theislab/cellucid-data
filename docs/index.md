# Cellucid

**Interactive Single-Cell Data Visualization**

Cellucid renders millions of cells in real-time 3D, enabling interactive exploration of UMAP/tSNE embeddings with gene expression overlays, filtering, and KNN connectivity visualization.

::::{grid} 1 1 3 3
:gutter: 2

:::{grid-item}
```{button-link} https://cellucid.com
:color: primary
:expand:

{octicon}`play;1em` Live Demo
```
:::

:::{grid-item}
```{button-link} https://github.com/theislab/cellucid-python
:color: secondary
:expand:

{octicon}`mark-github;1em` GitHub
```
:::

:::{grid-item}
```{button-link} https://pypi.org/project/cellucid/
:color: secondary
:expand:

{octicon}`package;1em` PyPI
```
:::

::::

---

## Installation

```bash
pip install cellucid
```

## Quick Start

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Jupyter Notebook
```python
from cellucid import show_anndata

# In-memory, h5ad, or zarr
show_anndata(adata)
show_anndata("data.h5ad")
show_anndata("data.zarr")
```
:::

:::{grid-item-card} Pre-exported Data
```python
from cellucid import prepare, show

# Export once (recommended for large datasets / sharing)
X_umap = adata.obsm["X_umap"]  # shape: (n_cells, 2) or (n_cells, 3)

prepare(
    latent_space=adata.obsm.get("X_pca", X_umap),
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp.get("connectivities"),
    X_umap_2d=adata.obsm.get("X_umap_2d", X_umap if X_umap.shape[1] == 2 else None),
    X_umap_3d=adata.obsm.get("X_umap_3d", X_umap if X_umap.shape[1] == 3 else None),
    out_dir="./export",
    compression=6,
    var_quantization=8,
    obs_continuous_quantization=8,
)

show("./export")
```
:::
::::

---

## Documentation

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {octicon}`book;1.5em;sd-mr-1` User Guide
:link: user_guide/index
:link-type: doc

Step-by-step tutorials covering all 14 loading options (h5ad, zarr, pre-exported) from local demos to Jupyter integration.
:::

:::{grid-item-card} {octicon}`code;1.5em;sd-mr-1` API Reference
:link: api/index
:link-type: doc

Complete reference for all functions and classes in the Cellucid Python package.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Examples
:link: examples/index
:link-type: doc

Real-world examples with various single-cell datasets ready for visualization.
:::

:::{grid-item-card} {octicon}`people;1.5em;sd-mr-1` Contributing
:link: contributing
:link-type: doc

Guidelines for contributing to Cellucid development.
:::

:::{grid-item-card} {octicon}`history;1.5em;sd-mr-1` Changelog
:link: changelog
:link-type: doc

Version history and release notes.
:::

:::{grid-item-card} {octicon}`mark-github;1.5em;sd-mr-1` GitHub
:link: https://github.com/theislab/cellucid-python

Source code, issues, and discussions.
:::
::::

---

## Features

- **Real-time 3D rendering** - Millions of cells with adaptive LOD
- **Gene expression overlays** - Efficient sparse matrix handling
- **Cell metadata coloring** - Categorical and continuous fields
- **Interactive filtering** - Cell selection and hiding
- **KNN connectivity** - Edge visualization
- **Multi-dimensional** - 1D timelines, 2D, 3D support
- **14 loading options** - h5ad, zarr, pre-exported binary formats

## AnnData Requirements

Your AnnData needs UMAP coordinates in `obsm`:

```python
# Required: one of these
adata.obsm['X_umap_3d']  # shape (n_cells, 3) - recommended
adata.obsm['X_umap_2d']  # shape (n_cells, 2)
adata.obsm['X_umap']     # shape (n_cells, 2 or 3)

# Optional
adata.obs                # Cell metadata
adata.X                  # Gene expression (dense or sparse)
adata.obsp['connectivities']  # KNN edges
```

## License

BSD-3-Clause

```{toctree}
:maxdepth: 2
:hidden:
:caption: Documentation

user_guide/index
api/index
examples/index
contributing
changelog
```
