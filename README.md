# TOAD_CLEWS

TOAD-based abrupt-shift detection and clustering pipeline for the Amazon
rainforest, driven by CMIP6 `1pctCO2` model output.

This project was created for the *Earth System Science and the Anthropocene*
lecture (Master's program: Climate, Earth, Water, Sustainability -- CLEWS)
at the University of Potsdam, taught by Johan Rockström and Jonathan Donges.
This analysis project was supervised by Jakob Harteg of the Earth Resilience
Science Unit ([ERSU](https://www.pik-potsdam.de/en/institute/futurelabs-science-units/ersu)),
Potsdam Institute for Climate Impact Research. Thank you, Jakob!

This pipeline is built on [TOAD](https://github.com/tipmip-methods/toad):

> Harteg, J., Roehrich, L., De Maeyer, K., Garbe, J., Sakschewski, B., Klose,
> A. K., Donges, J., Winkelmann, R., and Loriani, S.: TOAD: Tipping and Other
> Abrupt events Detector, Zenodo, https://doi.org/10.5281/zenodo.18316437, 2026.

## Setup

```bash
pip install -r requirements.txt
```

(Optionally, install into a virtual environment first: `python -m venv venv && source venv/bin/activate`.)

Then open the `.ipynb` notebooks directly in VS Code (Jupyter extension required).

## Verification

After installation, verify the environment works:

```python
python -c "import toad; import xarray; import geopandas; print('All imports successful!')"
```

## Project Structure

```
TOAD_CLEWS/
├── Amazon/
│   ├── configs/              # Configuration JSON files (model/variable settings)
│   ├── processed_data/       # Preprocessed NetCDF files
│   ├── plots/                # Visualization outputs
│   ├── results_toad/         # TOAD shift detection & clustering results
│   ├── docs/
│   │   ├── analysis_narrative.ipynb      # Main analysis notebook
│   │   └── TOAD_Amazon_interactive.ipynb # Interactive exploration
│   ├── preprocess.py         # Preprocessing pipeline
│   ├── toad_runner.py        # TOAD orchestration
│   └── run_pipeline.py       # Main orchestrator
├── data/                     # Third-party input data (not redistributed, see below)
└── README.md                 # This file
```

## Running the Pipeline

`Amazon/processed_data/` already ships preprocessed `.nc` files for every
model/variable used in this analysis, so you can go straight to
`--stage toad-only` without needing raw CMIP6 data at all.

```bash
cd Amazon

# Run complete pipeline (preprocess + TOAD)
python run_pipeline.py

# Run only preprocessing
python run_pipeline.py --stage preprocess-only

# Run only TOAD (skip preprocessing if NC files exist)
python run_pipeline.py --stage toad-only

# Process specific model and variable
python run_pipeline.py --model GFDL-ESM4 --variable cVeg --stage toad-only
```

### Extending the Pipeline

Adding a new variable or model requires no code changes -- just edit the
JSON config to add a new entry under `variable_settings`, or drop a new
`config_NEWMODEL.json` into `configs/`. `run_pipeline.py` automatically
discovers every config file and processes every variable defined in it.

## Data & Licensing

- **Code** (this repository): [MIT License](LICENSE)
- **TOAD package**: [BSD-2-Clause](https://github.com/tipmip-methods/toad) (pip dependency, not vendored)
- **CMIP6 model output** (GFDL-ESM4, TaiESM1): CC BY-SA 4.0, not redistributed
  raw -- see [`data/cmip6/README.md`](data/cmip6/README.md) for citation.
  `Amazon/processed_data/` ships the cropped derivative under the same terms.
- **Amazon biome/basin boundary shapefiles**: third-party, not redistributed
  -- see [`data/biome/README.md`](data/biome/README.md) and
  [`data/basin/README.md`](data/basin/README.md).
