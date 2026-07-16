# Installation Guide

## Setup

```bash
# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

Then open the `.ipynb` notebooks directly in VS Code (Jupyter extension required) and select the `venv` folder as the kernel/interpreter.

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
└── INSTALL.md               # This file
```

## Running the Pipeline

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

## Common Issues

### Memory issues with large climate datasets
Use Dask for out-of-core processing:
```python
import dask.array as da
import xarray as xr
```

## Questions?

Refer to the Jupyter notebooks for detailed analysis and examples.
