# When, Where, and How Much Does the Amazon Rainforest Change Under Global Warming? : Applying TOAD

TOAD-based abrupt-shift detection and clustering pipeline for the Amazon
rainforest, driven by CMIP6 `1pctCO2` model output.

**Authors**: Wonjin Kim, Bianca Pinto, Sekar Ayu Kinasih

This project was developed as part of the Earth System Science and the Anthropocene course in the Master's programme Climate, Earth, Water, Sustainability (CLEWS) at the University of Potsdam, taught by Johan Rockström and Jonathan Donges.

The project was supervised by Jakob Harteg from the Earth Resilience Science Unit (ERSU) at the Potsdam Institute for Climate Impact Research (PIK). Many thanks to Jakob for his guidance and support throughout the project!

This pipeline is built on [TOAD](https://github.com/tipmip-methods/toad):

> Harteg, J., Roehrich, L., De Maeyer, K., Garbe, J., Sakschewski, B., Klose,
> A. K., Donges, J., Winkelmann, R., and Loriani, S.: TOAD: Tipping and Other
> Abrupt events Detector, Zenodo, https://doi.org/10.5281/zenodo.18316437, 2026.

## Setup

```bash
pip install -r requirements.txt
```

(Optionally, install into a virtual environment first: `python -m venv venv && source venv/bin/activate`.)

Then open the `.ipynb` notebooks directly in VS Code (Jupyter extension required):

- [`Amazon/docs/analysis_narrative.ipynb`](Amazon/docs/analysis_narrative.ipynb) -- the full analysis narrative and results. Start here to see what we found!
- [`Amazon/docs/TOAD_Amazon_interactive.ipynb`](Amazon/docs/TOAD_Amazon_interactive.ipynb) -- an interactive sandbox to run TOAD yourself and tune shift/clustering parameters.

## Data

This repository does not redistribute third-party input data -- download what
your run needs. Each entry below covers what it's for, when you need it, and
its license in one place.

- **Shapefiles** (`data/biome/`, `data/basin/`) -- required for *any* run,
  including `--stage toad-only` against the `.nc` files already shipped in
  `Amazon/processed_data/`, since result plots draw the biome/basin boundary
  as an overlay. Third-party, not redistributed here -- see
  [`data/biome/README.md`](data/biome/README.md) and
  [`data/basin/README.md`](data/basin/README.md) for download links, file
  layout, and license terms.
- **Raw CMIP6 data** (`data/cmip6/`) -- `Amazon/processed_data/` already
  ships the preprocessed (cropped) `.nc` files used in this analysis, so you
  only need this if you want to run preprocessing yourself from raw model
  output (the default full pipeline, or `--stage preprocess-only`). CC
  BY-SA 4.0 (GFDL-ESM4, TaiESM1), not redistributed raw here -- see
  [`data/cmip6/README.md`](data/cmip6/README.md) for the download source,
  file naming, and citation.

## Running the Pipeline

See [Data](#data) above for what each `--stage` below needs.

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

Adding a new variable or model requires no code changes:

- **New variable**: add an entry under `variable_settings` in the relevant `config_{MODEL}.json`.
- **New model**: drop a new `config_NEWMODEL.json` into `Amazon/configs/`, and place its raw CMIP6 data under `data/cmip6/NEWMODEL/` (see [Data](#data) and [`data/cmip6/README.md`](data/cmip6/README.md) for the expected folder layout).

`run_pipeline.py` automatically discovers every config file and processes every variable defined in it.

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
├── data/                     # Third-party input data (not redistributed, see Data)
├── requirements.txt          # Python dependencies
├── LICENSE                   # MIT License (code)
└── README.md                 # This file
```

## License

- **Code** (this repository): [MIT License](LICENSE)
- **TOAD package**: [BSD-2-Clause](https://github.com/tipmip-methods/toad) (pip dependency, not vendored)

Third-party input data licensing is covered in [Data](#data) above.
