# Server Functions

```{eval-rst}
.. currentmodule:: cellucid
```

Functions and classes for running a local HTTP visualization server.

---

## CLI Command

The unified `cellucid serve` command auto-detects the data format:

```bash
# Serve any data - format auto-detected
cellucid serve /path/to/data.h5ad      # h5ad file
cellucid serve /path/to/data.zarr      # zarr store
cellucid serve /path/to/export         # pre-exported data

# With options
cellucid serve data.h5ad --port 9000 --no-browser

# Show help
cellucid serve --help
```

### CLI Options

| Option | Description | Default |
|--------|-------------|---------|
| `--port, -p` | Port to serve on | 8765 |
| `--host, -H` | Host to bind to | 127.0.0.1 |
| `--no-browser` | Don't auto-open browser | False |
| `--quiet, -q` | Suppress info messages | False |
| `--no-backed` | Load entire file into memory (h5ad/zarr) | False |
| `--latent-key` | Key in obsm for latent space | Auto-detected |

---

## Functions

### serve_anndata

Serve an AnnData object via HTTP for viewing in a browser. Supports in-memory objects, `.h5ad` files, and `.zarr` stores.

```python
from cellucid import serve_anndata

# Starts server and opens browser
serve_anndata(adata)
serve_anndata("data.h5ad")
serve_anndata("data.zarr")
```

```{eval-rst}
.. autofunction:: serve_anndata
```

---

### serve

Serve a pre-exported Cellucid dataset via HTTP. Use this when you've already run {func}`~cellucid.prepare`.

```python
from cellucid import prepare, serve

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
serve("./my_export")
```

```{eval-rst}
.. autofunction:: serve
```

---

## Classes

### AnnDataServer

Server class for AnnData with fine-grained control over the server lifecycle.

```python
from cellucid import AnnDataServer

server = AnnDataServer(adata)
server.start()
# ... do other work ...
server.stop()
```

```{eval-rst}
.. autoclass:: AnnDataServer
   :members:
   :show-inheritance:
```

---

### CellucidServer

Server class for pre-exported data with fine-grained control.

```python
from cellucid import CellucidServer

server = CellucidServer("./my_export")
server.start()
# ... do other work ...
server.stop()
```

```{eval-rst}
.. autoclass:: CellucidServer
   :members:
   :show-inheritance:
```

---

## See Also

- {func}`~cellucid.show_anndata` - Display in Jupyter instead of browser
- {func}`~cellucid.prepare` - Export data for static hosting
