# Installation Guide

## Option 1: Using Conda (Recommended)

```bash
# Create and activate the environment
conda env create -f environment.yml
conda activate toad-amazon

# Start Jupyter Lab
jupyter lab
```

## Option 2: Using pip

```bash
# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Start Jupyter Lab
jupyter lab
```

## Option 3: Using both (if Conda for system dependencies, pip for specific packages)

```bash
# Create minimal conda environment with system libraries
conda create -n toad-amazon python=3.13 geopandas cartopy rasterio gdal netCDF4

conda activate toad-amazon

# Install remaining packages via pip
pip install -r requirements.txt
```

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

### NetCDF4 or GDAL import errors on macOS
```bash
# Use conda to manage these system dependencies properly
conda install -c conda-forge gdal netCDF4
```

### Cartopy data download issues
Cartopy may need to download map data on first use. This requires internet connection.

### Memory issues with large climate datasets
Use Dask for out-of-core processing:
```python
import dask.array as da
import xarray as xr
```

## Questions?

Refer to the Jupyter notebooks for detailed analysis and examples.
