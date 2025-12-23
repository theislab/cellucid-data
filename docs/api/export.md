# Data Export

```{eval-rst}
.. currentmodule:: cellucid
```

Export arrays (typically sourced from an AnnData object) to the Cellucid on-disk format used by the WebGL viewer and the 14 loading paths (static hosting, server, file picker).

---

## prepare

Write a dataset directory containing:

- `dataset_identity.json` (metadata + embedding + optional vector field manifests)
- `points_1d.bin(.gz)`, `points_2d.bin(.gz)`, `points_3d.bin(.gz)` (any/all dims you provide)
- `obs_manifest.json` + `obs/` binaries
- `var_manifest.json` + `var/` binaries (if `gene_expression` is provided)
- `connectivity_manifest.json` + `connectivity/` binaries (if `connectivities` is provided)
- `vectors/` binaries (if `vector_fields` is provided)

```python
from cellucid import prepare

X_umap = adata.obsm["X_umap"]  # shape: (n_cells, 2) or (n_cells, 3)

vector_fields = {
    # Optional: per-cell displacement vectors in embedding space.
    # Preferred naming: <field>_umap_<dim>d (e.g. velocity_umap_2d, T_fwd_umap_3d)
    "velocity_umap_2d": adata.obsm.get("velocity_umap_2d"),
    "velocity_umap_3d": adata.obsm.get("velocity_umap_3d"),
    "T_fwd_umap_2d": adata.obsm.get("T_fwd_umap_2d"),
}

prepare(
    latent_space=adata.obsm.get("X_pca", X_umap),
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp.get("connectivities"),

    # Embeddings (provide any/all of 1D/2D/3D)
    X_umap_2d=adata.obsm.get("X_umap_2d", X_umap if X_umap.shape[1] == 2 else None),
    X_umap_3d=adata.obsm.get("X_umap_3d", X_umap if X_umap.shape[1] == 3 else None),

    # Optional: vector fields for the animated overlay
    vector_fields=vector_fields,

    out_dir="./my_export",
    compression=6,
    var_quantization=8,
    obs_continuous_quantization=8,
)
```

```{eval-rst}
.. autofunction:: prepare
```

---

## Output Format

```
my_export/
├── dataset_identity.json
├── obs_manifest.json
├── var_manifest.json                 # optional (gene expression)
├── connectivity_manifest.json        # optional (KNN edges)
├── points_1d.bin.gz                  # optional
├── points_2d.bin.gz                  # optional
├── points_3d.bin.gz                  # optional
├── obs/                              # obs field binaries
├── var/                              # gene expression binaries
├── connectivity/                     # KNN edge binaries
└── vectors/                          # optional (vector field binaries)
```

---

## See Also

- {func}`~cellucid.show` - Display exported data in Jupyter
- {func}`~cellucid.serve` - Serve exported data via HTTP
