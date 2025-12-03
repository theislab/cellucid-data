<p>
  <img src="cellucid-logo.svg" alt="Cellucid logo" width="320">
</p>

[![PyPI version](https://img.shields.io/pypi/v/cellucid.svg)](https://pypi.org/project/cellucid/)
[![Python versions](https://img.shields.io/pypi/pyversions/cellucid.svg)](https://pypi.org/project/cellucid/)
[![Documentation Status](https://readthedocs.org/projects/cellucid/badge/?version=latest)](https://cellucid.readthedocs.io/en/latest/)

# cellucid (Python)

Python utilities for exporting AnnData/numpy payloads into the manifest + binary format consumed by the Cellucid WebGL viewer. Everything here is independent of the website; generate data wherever you like, then copy the resulting files into the viewer’s `assets/exports/` directory when you host the site.

## Layout (repo: theislab/cellucid-data)
- `src/cellucid/prepare_data.py`: `export_data_for_web` writes points/obs/var/connectivity binaries into an `exports/` directory (configurable).
- `src/cellucid/anndata_variations.py`: UMAP sweep helper for experimentation (uses Scanpy).
- `src/cellucid/repo_to_text.py`: utility to dump repository contents to a single text file.
- `notebooks/`: exploratory notebook that demonstrates exporting viewer assets.
- `data/`: example/raw inputs (h5ad, pickles) and experiment outputs live here; this folder stays untracked by default.

## Installation
Requires Python 3.10+.

From the repo root:

```bash
pip install -e python
# or with optional extras
pip install -e python[analysis,docs]
```

The package name is `cellucid` and imports under `cellucid`.

## Usage
```python
import scanpy as sc
from cellucid.prepare_data import export_data_for_web

adata = sc.read_h5ad("path/to/adata.h5ad")
export_data_for_web(
    X_umap=adata.obsm["X_umap"],
    latent_space=adata.X,
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp.get("connectivities"),
    out_dir="exports",  # default; change as needed
    var_quantization=8,
    obs_continuous_quantization=8,
    compression=6,
)
```

The exporter defaults to `exports/` under your current working directory. After generation, copy or sync that folder into the website’s `assets/exports/` directory so the viewer can load the files.

To run the UMAP sweep helper:

```bash
python -m cellucid.anndata_variations
```

Ensure `python/data/raw` contains the expected h5ad/pickle files or override the module constants before calling. Paths in the notebook and sweep helper assume the `data/` directory inside this `python/` project; point them to your own locations as needed.

## Companion web viewer (theislab/cellucid)
The static viewer lives in a separate repository (`theislab/cellucid`). After running `export_data_for_web`, copy or sync the resulting `exports/` directory into that repo’s `assets/exports/` folder before hosting. The manifest/binary schema is shared; no build step is required on the web side.

## Documentation
- Build locally: `pip install -e python[docs] && sphinx-build -b html docs docs/_build/html`
- Read the Docs: `python/.readthedocs.yaml` targets `docs/conf.py` with the `docs` extra.

## Publishing
Install packaging tools if needed: `pip install build twine`

1. From `python/`, bump the version in `pyproject.toml`.
2. Build: `python -m build`
3. Upload to PyPI (or TestPyPI) with `twine upload dist/*`.

Update the `license`, `authors`, and `Repository` URL in `pyproject.toml` before publishing if needed.
