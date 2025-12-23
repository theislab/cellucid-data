# Jupyter Functions

```{eval-rst}
.. currentmodule:: cellucid
```

Functions for interactive visualization directly in Jupyter notebooks.

---

## show_anndata

Visualize an AnnData object directly without pre-exporting. Supports in-memory objects, `.h5ad` files, and `.zarr` stores.

```python
from cellucid import show_anndata

# In-memory AnnData
show_anndata(adata)

# From file path
show_anndata("data.h5ad")
show_anndata("data.zarr")
```

```{eval-rst}
.. autofunction:: show_anndata
```

---

## show

Display a pre-exported Cellucid dataset in Jupyter. Use this when you've already run {func}`~cellucid.prepare`.

```python
from cellucid import prepare, show

# Export once (see `cellucid.prepare` docs for full options)
X_umap = adata.obsm["X_umap"]
prepare(
    latent_space=adata.obsm.get("X_pca", X_umap),
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    X_umap_2d=adata.obsm.get("X_umap_2d", X_umap if X_umap.shape[1] == 2 else None),
    X_umap_3d=adata.obsm.get("X_umap_3d", X_umap if X_umap.shape[1] == 3 else None),
    out_dir="./my_export",
    compression=6,
)
show("./my_export")
```

```{eval-rst}
.. autofunction:: show
```

---

## See Also

- {class}`~cellucid.AnnDataViewer` - Class-based interface for `show_anndata`
- {class}`~cellucid.CellucidViewer` - Class-based interface for `show`
- {func}`~cellucid.serve_anndata` - Serve AnnData via HTTP instead of embedding
