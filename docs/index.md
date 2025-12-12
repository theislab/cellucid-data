# cellucid

Utilities to export AnnData/numpy payloads for the 3D UMAP WebGL viewer.

## Quickstart
- Install: `pip install -e python[docs,analysis]`
- Export data: `from cellucid import export_data_for_web`
- Default outputs land in `exports/` under your working directory; copy them into the viewer’s `assets/exports/` folder when hosting the site, or override `out_dir` to point there directly.

Companion viewer: static site in `theislab/cellucid` expects the exported `assets/exports/` folder.

## API
```{autofunction} cellucid.export_data_for_web
```
```{automodule} cellucid.anndata_variations
   :members:
   :undoc-members:
   :show-inheritance:
```

## Local docs build
```bash
cd python
pip install -e .[docs]
sphinx-build -b html docs docs/_build/html
```
