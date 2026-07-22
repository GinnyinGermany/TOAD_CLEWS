import glob
import os
import json
import numpy as np
import xarray as xr
import geopandas as gpd
import warnings
import matplotlib.pyplot as plt
from scipy import ndimage

# Suppress specific xarray/dask warnings for cleaner output
warnings.filterwarnings("ignore")

# ==========================================
# Project root directory (handles relative paths correctly)
# ==========================================
PREPROCESS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(PREPROCESS_DIR)

# ==========================================
# Phase 0: JSON Config Loader
# ==========================================
# All per-variable behaviour (preprocessing type/aggregation/unit_conversion,
# plot styling, TOAD shift/optimize settings) is driven by this JSON file
# instead of being hardcoded in the pipeline classes below.
DEFAULT_CONFIG_PATH = os.path.join(PREPROCESS_DIR, "config.json")

# Frequency directory CMIP6 uses per raw variable. Extend this as new raw
# variables are added to the JSON's "inputs"/target list.
VARIABLE_FREQUENCY = {
    "tas": "Amon",
    "pr": "Amon",
    "evspsbl": "Amon",
}
DEFAULT_FREQUENCY = "Lmon"

# Some entries in "variable_settings" are *derived* names (e.g. "AnnualRainfall")
# that do not correspond 1:1 to a raw CMIP6 variable id ("pr"). If the JSON
# does not carry an explicit "source_variable" key inside "preprocessing",
# this table is used as a fallback. Prefer adding "source_variable" to the
# JSON going forward -- that keeps everything fully config-driven.
DERIVED_SOURCE_VARIABLES = {
    "AnnualRainfall": "pr",
}


def load_pipeline_json(config_path=DEFAULT_CONFIG_PATH):
    """Load the pipeline's JSON configuration file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Pipeline config JSON not found: {config_path}")
    with open(config_path, "r") as f:
        return json.load(f)


def get_variable_config(pipeline_json, target_variable):
    """Return the variable_settings block for a given variable name."""
    settings = pipeline_json.get("variable_settings", {})
    if target_variable not in settings:
        raise KeyError(
            f"No configuration found for variable '{target_variable}' in JSON config. "
            f"Available variables: {list(settings.keys())}"
        )
    return settings[target_variable]


def resolve_source_variable(target_variable, preprocessing_config):
    """Resolve the raw CMIP6 variable id backing a (possibly derived) target variable."""
    return (
        preprocessing_config.get("source_variable")
        or DERIVED_SOURCE_VARIABLES.get(target_variable, target_variable)
    )


def frequency_for(variable):
    """Resolve the CMIP6 table frequency (Amon/Lmon/...) for a raw variable id."""
    return VARIABLE_FREQUENCY.get(variable, DEFAULT_FREQUENCY)


def seconds_per_month(time_coord):
    """Seconds in each month of a monthly time coordinate (for flux -> depth/total conversions)."""
    return time_coord.dt.days_in_month.astype("float64") * 86400.0


def flux_units_to_annual_total_units(original_units):
    """Best-effort conversion of a per-second flux unit string to an annual-total unit string."""
    if not original_units:
        return original_units
    stripped = (
        original_units.replace("s-1", "")
        .replace("s^-1", "")
        .replace("s**-1", "")
        .strip()
    )
    return f"{stripped} yr-1" if stripped else original_units


def compute_buffered_bbox(shapefile_path, buffer_degrees=5.0, bioma_filter="Amazonía"):
    """Computes a rectangular (lat_min, lat_max, lon_min, lon_max) bounding box
    from a biome shapefile's own extent plus a fixed-degree buffer margin,
    with longitudes converted to the 0-360 convention used by the raw CMIP6
    files (DataFetcher does no cropping itself; see SimpleProcessor).
    """
    boundary = gpd.read_file(shapefile_path)
    if bioma_filter is not None and "bioma" in boundary.columns:
        boundary = boundary[boundary["bioma"] == bioma_filter]
    if boundary.crs is not None and boundary.crs.to_string() != "EPSG:4326":
        boundary = boundary.to_crs(epsg=4326)

    lon_min, lat_min, lon_max, lat_max = boundary.total_bounds
    lat_min -= buffer_degrees
    lat_max += buffer_degrees
    lon_min -= buffer_degrees
    lon_max += buffer_degrees

    # Convert from -180/180 (shapefile convention) to 0-360 (raw CMIP6 convention)
    lon_min_360 = lon_min + 360 if lon_min < 0 else lon_min
    lon_max_360 = lon_max + 360 if lon_max < 0 else lon_max

    return lat_min, lat_max, lon_min_360, lon_max_360


DEFAULT_BIOME_SHAPEFILE_PATH = os.path.join(PROJECT_ROOT, "data", "Biome", "Biomas.shp")
DEFAULT_BASIN_SHAPEFILE_PATH = os.path.join(
    PROJECT_ROOT, "data", "AmazonBasinLimits", "amazon_sensulatissimo_gmm_v1.shp"
)


def load_boundary_shapefiles(
    biome_path=DEFAULT_BIOME_SHAPEFILE_PATH,
    basin_path=DEFAULT_BASIN_SHAPEFILE_PATH,
    bioma_filter="Amazonía",
):
    """Loads the Amazon biome and basin boundary shapefiles (both EPSG:4326),
    for use as reference outlines on plots -- NOT for masking/clipping data.
    """
    biome = gpd.read_file(biome_path)
    if bioma_filter is not None and "bioma" in biome.columns:
        biome = biome[biome["bioma"] == bioma_filter]
    if biome.crs is not None and biome.crs.to_string() != "EPSG:4326":
        biome = biome.to_crs(epsg=4326)

    basin = gpd.read_file(basin_path)
    if basin.crs is not None and basin.crs.to_string() != "EPSG:4326":
        basin = basin.to_crs(epsg=4326)

    return biome, basin


def biome_cell_mask(lat_values, lon_values, biome_gdf=None):
    """Boolean (lat, lon) mask of grid cells whose center falls inside the
    Amazon biome polygon. Longitudes may be in 0-360 or -180/180 convention.

    Unlike draw_boundary_overlays(), this IS meant for masking/averaging data
    -- use it for summary statistics (e.g. area-averaged trend plots) that
    should represent the Amazon biome specifically, not the wider rectangular
    processing domain.
    """
    from shapely.geometry import Point

    if biome_gdf is None:
        biome_gdf, _ = load_boundary_shapefiles()
    biome_union = biome_gdf.union_all()

    mask = np.zeros((len(lat_values), len(lon_values)), dtype=bool)
    for i, lat in enumerate(lat_values):
        for j, lon in enumerate(lon_values):
            lon180 = lon - 360 if lon > 180 else lon
            mask[i, j] = biome_union.contains(Point(lon180, lat))
    return mask


def draw_boundary_overlays(ax, lon_is_0_360, biome_gdf=None, basin_gdf=None, transform=None):
    """Draws the Amazon biome (solid black) and basin (dashed gray) outlines
    on a matplotlib Axes as reference lines only (no data masking).

    Args:
        ax: matplotlib (or cartopy GeoAxes) Axes to draw on.
        lon_is_0_360: True if the plotted data uses 0-360 longitude
            convention (raw CMIP6 files do); shifts the -180/180 shapefile
            geometries to match before plotting.
        biome_gdf, basin_gdf: pre-loaded GeoDataFrames (EPSG:4326). If either
            is None, loads the default shapefiles via load_boundary_shapefiles().
        transform: optional cartopy CRS (e.g. ccrs.PlateCarree()) to pass
            through to GeoDataFrame.plot(). Required when ax is a cartopy
            GeoAxes (TOAD's own map()/overview() plots default to one),
            otherwise the overlay geometry is misaligned with the map.
    """
    if biome_gdf is None or basin_gdf is None:
        default_biome, default_basin = load_boundary_shapefiles()
        biome_gdf = biome_gdf if biome_gdf is not None else default_biome
        basin_gdf = basin_gdf if basin_gdf is not None else default_basin

    biome_plot = biome_gdf.translate(xoff=360) if lon_is_0_360 else biome_gdf
    basin_plot = basin_gdf.translate(xoff=360) if lon_is_0_360 else basin_gdf

    plot_kwargs = {} if transform is None else {"transform": transform}
    basin_plot.plot(ax=ax, facecolor="none", edgecolor="dimgray", linewidth=1.2, linestyle="--", **plot_kwargs)
    biome_plot.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=1.5, **plot_kwargs)


# ==========================================
# Phase 1: Configuration Loader
# ==========================================
class PipelineConfig:
    """Centralized configuration manager for the pipeline.

    Loads model-agnostic pipeline settings (paths, spatial bbox) as before,
    but now also loads the per-variable "preprocessing" / "plot" / "toad"
    settings dynamically from the JSON config file.
    """

    def __init__(self, model_name, target_variable, config_path=DEFAULT_CONFIG_PATH):
        self.model_name = model_name
        self.target_variable = target_variable
        self.config_path = config_path

        self.pipeline_json = load_pipeline_json(config_path)

        json_model = self.pipeline_json.get("model")
        if json_model and json_model != model_name:
            warnings.warn(
                f"Config model '{json_model}' does not match requested model '{model_name}'."
            )

        try:
            self.variable_config = get_variable_config(self.pipeline_json, target_variable)
        except KeyError:
            # Some variables (e.g. "tas", used only internally for GWL) are not
            # expected to have an entry in variable_settings. Fall back to an
            # empty config rather than hard-failing.
            warnings.warn(
                f"'{target_variable}' has no entry in variable_settings; "
                "proceeding with empty preprocessing/plot/toad settings."
            )
            self.variable_config = {}

        self.preprocessing_config = self.variable_config.get("preprocessing", {})
        self.plot_config = self.variable_config.get("plot", {})
        self.toad_config = self.variable_config.get("toad", {})

        # Hardcoding built-in paths
        self.base_path = os.path.join(PROJECT_ROOT, "CMIP6 data")
        self.shapefile_path = os.path.join(PROJECT_ROOT, "data", "Biome", "Biomas.shp")
        self.basin_shapefile_path = os.path.join(
            PROJECT_ROOT, "data", "AmazonBasinLimits", "amazon_sensulatissimo_gmm_v1.shp"
        )

        # Spatial bounding box: the Amazonía biome polygon's own bounding box,
        # plus a fixed-degree buffer margin, rather than a hand-picked box or
        # a hard mask to the irregular biome polygon. Clipping to the jagged
        # biome polygon (previously done in SpatialMasker) removes real
        # spatial neighbors right at the boundary, which artificially lowers
        # local point-density there and fragments density-based clustering
        # (SpaceTimeDBSCAN) near the edges. A rectangular box (with margin)
        # keeps every boundary grid cell's true neighbors in the data, while
        # the biome/basin shapes are still drawn as reference outlines on
        # plots (see plot_spatial_maps_from_nc / toad_runner overlays).
        (
            self.lat_min,
            self.lat_max,
            self.lon_min,
            self.lon_max,
        ) = compute_buffered_bbox(self.shapefile_path, buffer_degrees=5.0)

    def get_dict(self):
        """Returns the configuration as a dictionary."""
        return {
            "model_name": self.model_name,
            "target_variable": self.target_variable,
            "base_path": self.base_path,
            "shapefile_path": self.shapefile_path,
            "basin_shapefile_path": self.basin_shapefile_path,
            "lat_min": self.lat_min,
            "lat_max": self.lat_max,
            "lon_min": self.lon_min,
            "lon_max": self.lon_max,
            "preprocessing": self.preprocessing_config,
            "plot": self.plot_config,
            "toad": self.toad_config,
        }


# ==========================================
# Phase 2: Data Fetcher
# ==========================================
class DataFetcher:
    """Dynamically loads NetCDF files and handles safe time alignment.

    Which raw files get loaded is now driven entirely by the JSON config:
    - "simple" variables load a single raw file (target variable, or its
      "source_variable"/DERIVED_SOURCE_VARIABLES mapping).
    - "custom" variables load every variable listed in preprocessing.inputs,
      aligning each one against the first input's time coordinate.
    """

    def __init__(self, config):
        self.config = config
        self.model = config["model_name"]
        self.base_dir = config["base_path"]
        self.preprocessing = config.get("preprocessing", {})

    def _build_filepath(self, variable, experiment, frequency=None):
        """Constructs the file path pattern and verifies existence."""
        frequency = frequency or frequency_for(variable)
        file_pattern = f"{variable}_{frequency}_{self.model}_{experiment}_r1i1p1f1_*.nc"
        full_path = os.path.join(self.base_dir, self.model, experiment, file_pattern)

        if not glob.glob(full_path):
            raise FileNotFoundError(f"No files found: {full_path}")
        return full_path

    def _safe_time_align(self, ds_target, ds_reference, target_name):
        """Validate and align monthly time coordinates between datasets."""
        if len(ds_target.time) != len(ds_reference.time):
            raise ValueError(f"Time length mismatch in {target_name}.")

        target_ym = ds_target.time.dt.strftime("%Y-%m").values
        ref_ym = ds_reference.time.dt.strftime("%Y-%m").values

        if not np.array_equal(target_ym, ref_ym):
            raise ValueError(f"Year-Month values do not match in {target_name}.")

        ds_target["time"] = ds_reference["time"]
        return ds_target

    def fetch_data(self):
        """Returns a dictionary of required raw datasets."""
        datasets = {}
        target_var = self.config["target_variable"]

        tas_files = glob.glob(self._build_filepath("tas", "1pctCO2"))
        with xr.open_dataset(tas_files[0]) as ds_first:
            datasets["attrs_raw"] = ds_first.attrs.copy()

        # 1. Load temperature data for GWL
        datasets["tas_1pct"] = xr.open_mfdataset(self._build_filepath("tas", "1pctCO2"))
        datasets["tas_pi"] = xr.open_mfdataset(
            self._build_filepath("tas", "piControl"), combine="by_coords"
        )

        proc_type = self.preprocessing.get("type", "simple")

        if proc_type == "custom":
            # 2. Custom variable branch: fetch every variable listed in "inputs"
            input_vars = self.preprocessing.get("inputs", [])
            if not input_vars:
                raise ValueError(
                    f"Custom preprocessing for '{target_var}' requires a non-empty "
                    "'inputs' list in the JSON config."
                )

            reference_ds = None
            for i, var in enumerate(input_vars):
                ds = xr.open_mfdataset(
                    self._build_filepath(var, "1pctCO2"),
                    combine="by_coords",
                )
                if i == 0:
                    reference_ds = ds
                else:
                    ds = self._safe_time_align(ds, reference_ds, var)
                datasets[var] = ds

        elif proc_type == "simple":
            # 3. Simple variable branch: single raw file
            source_var = resolve_source_variable(target_var, self.preprocessing)
            datasets["target_raw"] = xr.open_mfdataset(
                self._build_filepath(source_var, "1pctCO2")
            )

        else:
            raise ValueError(
                f"Unknown preprocessing type '{proc_type}' for variable '{target_var}'."
            )

        return datasets


# ==========================================
# Phase 3: Core Processors
# ==========================================
class BaseProcessor:
    def __init__(self, config):
        self.config = config
        self.target_var = config["target_variable"]
        self.preprocessing = config.get("preprocessing", {})
        self.lat_slice = slice(config["lat_min"], config["lat_max"])
        self.lon_slice = slice(config["lon_min"], config["lon_max"])

    def process(self, datasets):
        raise NotImplementedError


class SimpleProcessor(BaseProcessor):
    """Generic processor for "simple"-type variables.

    Behaviour is entirely driven by config["preprocessing"]:
      - time_frequency: resample frequency (e.g. "1YS")
      - aggregation: "mean" or "sum"
      - unit_conversion: "none" or "flux_to_annual_total"

    This single class replaces the old hardcoded DirectVariableProcessor and
    AnnualRainfallProcessor -- both cVeg-style mean variables and fFire /
    AnnualRainfall-style flux-integrated-sum variables are handled the same
    way, just with different config values.
    """

    AGGREGATIONS = {"mean": "mean", "sum": "sum"}

    def process(self, datasets):
        source_var = resolve_source_variable(self.target_var, self.preprocessing)
        raw_da = datasets["target_raw"][source_var]
        original_attrs = raw_da.attrs.copy()

        sliced_da = raw_da.sel(lat=self.lat_slice, lon=self.lon_slice)

        unit_conversion = self.preprocessing.get("unit_conversion", "none")
        aggregation = self.preprocessing.get("aggregation", "mean")
        time_frequency = self.preprocessing.get("time_frequency", "1YS")

        if unit_conversion == "none":
            converted_da = sliced_da
        elif unit_conversion == "flux_to_annual_total":
            # Convert a per-second flux to a per-month total before summing over the year.
            converted_da = sliced_da * seconds_per_month(sliced_da["time"])
        else:
            warnings.warn(
                f"Unknown unit_conversion '{unit_conversion}' for '{self.target_var}'; "
                "no conversion applied."
            )
            converted_da = sliced_da

        if aggregation not in self.AGGREGATIONS:
            raise ValueError(
                f"Unsupported aggregation '{aggregation}' for '{self.target_var}'. "
                f"Supported: {list(self.AGGREGATIONS.keys())}"
            )

        resampled = converted_da.resample(time=time_frequency)
        annual_da = getattr(resampled, self.AGGREGATIONS[aggregation])().load()

        annual_da.name = self.target_var
        annual_da.attrs = original_attrs
        if unit_conversion == "flux_to_annual_total":
            annual_da.attrs["units"] = flux_units_to_annual_total_units(
                original_attrs.get("units", "")
            )
            annual_da.attrs["long_name"] = (
                f"Annual total of {original_attrs.get('long_name', self.target_var)}"
            )
        annual_da.attrs["unit_conversion"] = unit_conversion
        annual_da.attrs["aggregation"] = aggregation

        return annual_da


# ---- Custom formula registry -------------------------------------------
# Each formula function has the signature (datasets, config) -> xr.DataArray
# and is looked up dynamically by config["preprocessing"]["formula_id"].
# Adding a new custom variable to the JSON only requires registering a new
# formula function here -- no changes needed elsewhere in the pipeline.

def _validate_mcwd_inputs(pr, evspsbl):
    """Validate time coordinates and physical units for the MCWD formula."""
    pr_year_month = pr["time"].dt.strftime("%Y-%m").values
    et_year_month = evspsbl["time"].dt.strftime("%Y-%m").values

    if not np.array_equal(pr_year_month, et_year_month):
        raise ValueError("pr and evspsbl do not have matching monthly timestamps.")

    pr_units = pr.attrs.get("units", "")
    et_units = evspsbl.attrs.get("units", "")

    expected_units = {
        "kg m-2 s-1",
        "kg m-2 s^-1",
        "kg m**-2 s**-1",
    }

    if pr_units not in expected_units:
        warnings.warn(
            f"Unexpected pr units: {pr_units}. Expected a water flux in kg m-2 s-1."
        )

    if et_units not in expected_units:
        warnings.warn(
            f"Unexpected evspsbl units: {et_units}. Expected a water flux in kg m-2 s-1."
        )


def formula_mcwd_hydrological_year(datasets, config):
    """Calculate hydrological-year MCWD from precipitation and evapotranspiration.

    Registered under formula_id "MCWD_hydrological_year". Reads the two
    input variable names from config["preprocessing"]["inputs"] (order:
    [precipitation, evapotranspiration]) so the formula stays generic with
    respect to which raw variables back it.
    """
    preprocessing = config.get("preprocessing", {})
    inputs = preprocessing.get("inputs", ["pr", "evspsbl"])
    if len(inputs) < 2:
        raise ValueError(
            "formula_id 'MCWD_hydrological_year' requires two 'inputs': "
            "[precipitation_variable, evapotranspiration_variable]."
        )
    precip_var, evap_var = inputs[0], inputs[1]

    lat_slice = slice(config["lat_min"], config["lat_max"])
    lon_slice = slice(config["lon_min"], config["lon_max"])

    pr = datasets[precip_var][precip_var].sel(lat=lat_slice, lon=lon_slice)
    evspsbl = datasets[evap_var][evap_var].sel(lat=lat_slice, lon=lon_slice)

    _validate_mcwd_inputs(pr, evspsbl)

    # Ensure that both variables use exactly the same time coordinate
    evspsbl = evspsbl.assign_coords(time=pr["time"])

    # Convert monthly mean fluxes from kg m-2 s-1 to mm month-1
    monthly_precipitation = pr * seconds_per_month(pr["time"])
    monthly_evaporation = evspsbl * seconds_per_month(pr["time"])

    monthly_water_balance = (monthly_precipitation - monthly_evaporation).load()

    years = monthly_water_balance["time"].dt.year.values.astype(int)
    months = monthly_water_balance["time"].dt.month.values.astype(int)

    # Oct-Dec belong to the hydrological year ending in the next year
    hydrological_years = np.where(months >= 10, years + 1, years)

    time_frequency = preprocessing.get("time_frequency", "1YS")

    # Create annual timestamps using the same year-start convention as GWL
    annual_time = (
        pr.isel(lat=0, lon=0)
        .resample(time=time_frequency)
        .mean()
        ["time"]
    )

    time_lookup = {
        int(year): timestamp
        for year, timestamp in zip(
            annual_time.dt.year.values.astype(int),
            annual_time.values,
        )
    }

    annual_results = []

    for hydrological_year in np.unique(hydrological_years):
        indices = np.where(hydrological_years == hydrological_year)[0]

        # Skip incomplete Oct-Sep years
        if len(indices) != 12:
            continue

        water_balance_year = monthly_water_balance.isel(time=indices)
        template = water_balance_year.isel(time=0, drop=True)

        cumulative_water_deficit = np.zeros_like(template.values, dtype=np.float64)
        minimum_cwd = np.zeros_like(cumulative_water_deficit, dtype=np.float64)

        for time_index in range(12):
            monthly_balance = water_balance_year.isel(time=time_index).values

            cumulative_water_deficit = np.minimum(
                0.0, cumulative_water_deficit + monthly_balance
            )
            minimum_cwd = np.minimum(minimum_cwd, cumulative_water_deficit)

        # Keep negative MCWD values
        mcwd_values = minimum_cwd
        hydro_year = int(hydrological_year)

        if hydro_year not in time_lookup:
            continue

        mcwd_year = xr.DataArray(
            mcwd_values,
            coords=template.coords,
            dims=template.dims,
            name=config["target_variable"],
        ).expand_dims(time=[time_lookup[hydro_year]])

        annual_results.append(mcwd_year)

    mcwd_annual = xr.concat(annual_results, dim="time").sortby("time")

    mcwd_annual.name = config["target_variable"]
    mcwd_annual.attrs = {
        "long_name": (
            "Maximum cumulative water deficit "
            "for the October-September hydrological year"
        ),
        "units": "mm",
        "calculation": (
            "min(CWD); CWD_i = min("
            "0, CWD_i-1 + precipitation_i - evapotranspiration_i)"
        ),
        "hydrological_year": "October to September",
        "unit_conversion": preprocessing.get("unit_conversion", "none"),
    }

    return mcwd_annual


FORMULA_REGISTRY = {
    "MCWD_hydrological_year": formula_mcwd_hydrological_year,
}


class CustomProcessor(BaseProcessor):
    """Dispatches to a registered formula function based on config's "formula_id"."""

    def process(self, datasets):
        formula_id = self.preprocessing.get("formula_id")
        if not formula_id:
            raise ValueError(
                f"Custom preprocessing for '{self.target_var}' requires a 'formula_id' "
                "in the JSON config."
            )

        formula_fn = FORMULA_REGISTRY.get(formula_id)
        if formula_fn is None:
            raise KeyError(
                f"No formula registered for formula_id '{formula_id}'. "
                f"Available: {list(FORMULA_REGISTRY.keys())}"
            )

        return formula_fn(datasets, self.config)


class ProcessorFactory:
    """Routes data to the correct processor based on config["preprocessing"]["type"]."""

    @staticmethod
    def get_processor(config):
        proc_type = config.get("preprocessing", {}).get("type", "simple")

        if proc_type == "custom":
            return CustomProcessor(config)
        if proc_type == "simple":
            return SimpleProcessor(config)

        raise ValueError(
            f"Unknown preprocessing type '{proc_type}' for variable "
            f"'{config['target_variable']}'."
        )


class GWLCalculator:
    """Calculate GWL following the Terpstra et al. approach."""

    @staticmethod
    def calculate_global_annual_mean(ds):
        tas = ds["tas"]

        weights = np.cos(np.deg2rad(ds["lat"]))
        weights.name = "weights"

        global_monthly_tas = (
            tas
            .weighted(weights)
            .mean(dim=["lat", "lon"], skipna=True)
        )

        global_annual_tas = (
            global_monthly_tas
            .resample(time="1YS")
            .mean(dim="time")
            .load()
        )

        return global_annual_tas

    @staticmethod
    def compute_gwl(tas_1pct_ds, tas_pi_ds):
        # Annual global mean temperatures
        gmt_1pct = GWLCalculator.calculate_global_annual_mean(tas_1pct_ds)
        gmt_pi = GWLCalculator.calculate_global_annual_mean(tas_pi_ds)

        # Terpstra et al.: mean over the complete piControl simulation
        gmt_preindustrial = float(gmt_pi.mean(dim="time").values)

        # Terpstra et al.: Gaussian smoothing with sigma = 10 years
        gmt_1pct_smoothed_values = ndimage.gaussian_filter1d(
            gmt_1pct.values,
            sigma=10,
            mode="nearest",
        )

        gmt_1pct_smoothed = xr.DataArray(
            gmt_1pct_smoothed_values,
            dims=["time"],
            coords={"time": gmt_1pct["time"]},
            name="GMT_smoothed",
        )

        # GWL relative to the piControl climatological mean
        gwl = gmt_1pct_smoothed - gmt_preindustrial

        gwl.name = "GWL"
        gwl.attrs["units"] = "°C"
        gwl.attrs["baseline"] = (
            "Mean GMT over the complete piControl simulation"
        )
        gwl.attrs["smoothing"] = (
            "Gaussian filter with sigma=10 years, mode='nearest'"
        )

        print(f"piControl baseline: {gmt_preindustrial:.4f} K")
        print(
            "GWL range:",
            f"{float(gwl.min()):.3f} to "
            f"{float(gwl.max()):.3f} K",
        )

        return gwl



# ==========================================
# Phase 5: TOAD Formatter & Exporter
# ==========================================
class TOADExporter:
    """Formats and exports to a TOAD-compatible NetCDF file.

    Note: this deliberately does NOT read or persist config["toad"] (the
    shift/optimize settings) into the exported file. The JSON config is the
    single source of truth for those settings, and the TOAD-running stage
    (toad_runner.py) re-reads the same JSON directly instead of relying on
    anything embedded in the .nc file. This keeps preprocessing outputs
    decoupled from TOAD-run parameters, so re-tuning shift/optimize settings
    never requires regenerating the .nc file.
    """
    def __init__(self, config):
        self.model = config["model_name"]
        self.variable = config["target_variable"]

    def export(self, masked_da, gwl_da, raw_attrs, rolling_years=None):
        merged_ds = xr.merge([masked_da, gwl_da], join="inner")
        final_ds = merged_ds.set_coords("GWL").swap_dims({"time": "GWL"}).sortby("GWL")

        final_ds.attrs = raw_attrs
        final_ds.attrs["nominal_resolution"] = raw_attrs.get("nominal_resolution", "Not Available")
        final_ds.attrs["rolling_years"] = (
            "none" if rolling_years is None or rolling_years <= 1 else rolling_years
        )
        final_ds.attrs["processing"] = smoothing_title(rolling_years)

        vars_to_drop = [var for var in final_ds.data_vars if var != self.variable]
        final_ds = final_ds.drop_vars(vars_to_drop)

        output_dir = os.path.join(PREPROCESS_DIR, "processed_data", self.model)
        os.makedirs(output_dir, exist_ok=True)
        suffix = smoothing_suffix(rolling_years)
        output_filename = f"TOAD_{self.variable}_{self.model}{suffix}.nc"
        output_path = os.path.join(output_dir, output_filename)

        final_ds.to_netcdf(output_path)
        print(f"[Export Success] Created NetCDF at: {output_path}")
        return output_path


# ==========================================
# Function 1: Generate and Save .nc file
# ==========================================
def generate_nc(model_name, target_variable, rolling_years=None, config_path=DEFAULT_CONFIG_PATH):
    """Executes Phase 1 to 5 to process raw CMIP6 data and save it as a .nc file.

    rolling_years defaults to whatever is set in the JSON config's
    preprocessing.rolling_years for this variable; pass an explicit value to
    override it.
    """
    print(f"\n--- Generating .nc for {model_name} | {target_variable} ---")

    # Paths and per-variable settings are loaded dynamically from the JSON config.
    config_manager = PipelineConfig(model_name, target_variable, config_path=config_path)
    config = config_manager.get_dict()

    if rolling_years is None:
        rolling_years = config["preprocessing"].get("rolling_years")

    fetcher = DataFetcher(config)
    datasets = fetcher.fetch_data()

    processor = ProcessorFactory.get_processor(config)
    target_da = processor.process(datasets)

    target_da = apply_centered_rolling_mean(
        target_da,
        rolling_years,
    )

    gwl_da = GWLCalculator.compute_gwl(
        datasets["tas_1pct"],
        datasets["tas_pi"],
    )

    print("Target time size:", target_da.sizes.get("time", 0))
    print("GWL time size:", gwl_da.sizes.get("time", 0))

    print("Target first times:")
    print(target_da["time"].values[:3])

    print("GWL first times:")
    print(gwl_da["time"].values[:3])

    common_times = np.intersect1d(
        target_da["time"].values,
        gwl_da["time"].values,
    )

    print("Common time size:", len(common_times))

    exporter = TOADExporter(config)
    output_nc_path = exporter.export(
        target_da,
        gwl_da,
        datasets["attrs_raw"],
        rolling_years=rolling_years,
    )

    return output_nc_path

# ==========================================
# Function 2: Caculating centered rolling mean
# ==========================================

def apply_centered_rolling_mean(data_array, rolling_years):
    """Apply centered rolling mean along time dimension."""

    if rolling_years is None or rolling_years <= 1:
        return data_array

    smoothed = (
        data_array
        .rolling(
            time=rolling_years,
            center=True,
            min_periods=rolling_years,
        )
        .mean()
        .dropna(dim="time", how="all")
    )

    smoothed.name = data_array.name
    smoothed.attrs = data_array.attrs.copy()
    smoothed.attrs["processing"] = (
        f"{rolling_years}-year centered rolling mean"
    )
    smoothed.attrs["rolling_years"] = rolling_years

    return smoothed


def smoothing_suffix(rolling_years):
    """Return filename suffix for smoothing setting."""

    if rolling_years is None or rolling_years <= 1:
        return ""

    return f"_smooth{rolling_years}yr"


def smoothing_title(rolling_years):
    """Return plot title label for smoothing setting."""

    if rolling_years is None or rolling_years <= 1:
        return "Annual values"

    return f"{rolling_years}-year centered rolling mean"

# ==========================================
# Function 3: Load .nc file and Plot
# ==========================================
def plot_trend_from_nc(model_name, target_variable, rolling_years=None, config_path=DEFAULT_CONFIG_PATH):
    """Loads a previously processed .nc file and generates a spatial mean plot.

    Plot styling (color/ylabel) and the default rolling_years come from the
    JSON config's "plot" / "preprocessing" blocks for this variable.
    """
    pipeline_json = load_pipeline_json(config_path)
    try:
        var_cfg = get_variable_config(pipeline_json, target_variable)
    except KeyError:
        var_cfg = {}

    if rolling_years is None:
        rolling_years = var_cfg.get("preprocessing", {}).get("rolling_years")

    suffix = smoothing_suffix(rolling_years)
    file_path = os.path.join(PREPROCESS_DIR, "processed_data", model_name, f"TOAD_{target_variable}_{model_name}{suffix}.nc")
    if not os.path.exists(file_path):
        print(f"[Error] File not found: {file_path}")
        return

    print(f"\n--- Plotting data from {file_path} ---")
    ds = xr.open_dataset(file_path)

    # Calculate Spatial Mean over the Amazon biome only (the processing
    # domain is a wider rectangular buffer around the biome; averaging over
    # it directly would dilute the Amazon signal with non-Amazon cells).
    mask_2d = biome_cell_mask(ds["lat"].values, ds["lon"].values)
    mask_da = xr.DataArray(mask_2d, dims=["lat", "lon"], coords={"lat": ds["lat"], "lon": ds["lon"]})
    spatial_mean = ds[target_variable].where(mask_da).mean(dim=["lon", "lat"], skipna=True)

    meta = var_cfg.get("plot", {"color": "blue", "ylabel": f"Mean {target_variable}"})
    gwl_vals = spatial_mean["GWL"].values
    var_vals = spatial_mean.values

    plt.figure(figsize=(8, 5))
    plt.plot(
        gwl_vals,
        var_vals,
        color=meta.get("color", "blue"),
        linewidth=2.5,
        label=meta.get("label", f"{model_name} {target_variable}"),
    )

    plt.title(
        f"Amazon Area-Averaged {target_variable} Response to GWL "
        f"({model_name}, {smoothing_title(rolling_years)})",
        fontsize=13,
    )
    plt.xlabel("Global Warming Level (GWL) [°C]", fontsize=11)
    plt.ylabel(meta.get("ylabel", f"Mean {target_variable}"), fontsize=11)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    output_dir = os.path.join(PREPROCESS_DIR, "plots", model_name)
    os.makedirs(output_dir, exist_ok=True)
    output_filename = f"TrendPlot_{target_variable}_{model_name}{suffix}.png"
    output_path = os.path.join(output_dir, output_filename)

    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"[Plot Success] Saved graph at: {output_path}")

# ==========================================
# Function 4: Plot GWL graph for model
# ==========================================
def plot_gwl_by_model_year(model_name, config_path=DEFAULT_CONFIG_PATH):
    """Calculate and plot GWL against 1pctCO2 model year."""

    # "tas" has no entry in variable_settings (it's only used internally for
    # GWL) -- PipelineConfig tolerates that and falls back to empty settings.
    config = PipelineConfig(
        model_name=model_name,
        target_variable="tas",
        config_path=config_path,
    ).get_dict()

    fetcher = DataFetcher(config)

    tas_1pct = xr.open_mfdataset(
        fetcher._build_filepath("tas", "1pctCO2"),
        combine="by_coords",
    )

    tas_pi = xr.open_mfdataset(
        fetcher._build_filepath("tas", "piControl"),
        combine="by_coords",
    )

    gwl = GWLCalculator.compute_gwl(
        tas_1pct,
        tas_pi,
    )

    model_years = np.arange(1, gwl.sizes["time"] + 1)

    plt.figure(figsize=(8, 5))

    plt.plot(
        model_years,
        gwl.values,
        linewidth=2.5,
        label=model_name,
    )

    plt.axhline(
        0,
        color="black",
        linestyle="--",
        linewidth=1,
    )

    plt.title(
        f"Global Warming Level by Model Year ({model_name})",
        fontsize=13,
    )

    plt.xlabel("1pctCO2 Model Year", fontsize=11)
    plt.ylabel("Global Warming Level (GWL) [°C]", fontsize=11)

    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    output_dir = os.path.join(PREPROCESS_DIR, "plots", model_name)
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(
        output_dir,
        f"GWL_ModelYear_{model_name}.png",
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    tas_1pct.close()
    tas_pi.close()

    print(f"[Plot Success] Saved GWL graph at: {output_path}")


# ==========================================
# Function 5: Load .nc file and Plot Spatial Maps
# ==========================================
def plot_spatial_maps_from_nc(model_name, target_variable, rolling_years=None, config_path=DEFAULT_CONFIG_PATH):
    """Loads a processed .nc file and generates a 1x6 grid of spatial maps.

    cmap/label styling and the default rolling_years come from the JSON
    config's "plot" / "preprocessing" blocks for this variable.
    """
    pipeline_json = load_pipeline_json(config_path)
    try:
        var_cfg = get_variable_config(pipeline_json, target_variable)
    except KeyError:
        var_cfg = {}

    if rolling_years is None:
        rolling_years = var_cfg.get("preprocessing", {}).get("rolling_years")

    suffix = smoothing_suffix(rolling_years)
    file_path = os.path.join(PREPROCESS_DIR, "processed_data", model_name, f"TOAD_{target_variable}_{model_name}{suffix}.nc")
    if not os.path.exists(file_path):
        print(f"[Error] File not found: {file_path}")
        return

    print(f"\n--- Plotting spatial maps from {file_path} ---")
    ds = xr.open_dataset(file_path)
    data_array = ds[target_variable]

    meta = var_cfg.get("plot", {"cmap": "viridis", "label": f"{target_variable}"})

    # GFDL-ESM4: GWL 0.0 ~ 4.0
    target_gwls = [0.0, 1.0, 1.5, 2.0, 3.0, 4.0]

    # For other models, uncomment below and adjust GWL range accordingly
    # gwl_min = float(data_array["GWL"].min())
    # gwl_max = float(data_array["GWL"].max())
    # target_gwls = [gwl_min + (gwl_max - gwl_min) * i / 5 for i in range(6)]

    v_min = float(data_array.min(skipna=True))
    v_max = float(data_array.max(skipna=True))

    lon_is_0_360 = bool(data_array["lon"].max() > 180)
    biome_gdf, basin_gdf = load_boundary_shapefiles()

    fig, axes = plt.subplots(1, 6, figsize=(24, 4.5), sharex=True, sharey=True)

    for ax, target in zip(axes, target_gwls):
        da_slice = data_array.sel(GWL=target, method="nearest")
        actual_gwl = float(da_slice["GWL"])

        mesh = ax.pcolormesh(
            da_slice["lon"],
            da_slice["lat"],
            da_slice.values,
            cmap=meta.get("cmap", "viridis"),
            vmin=v_min,
            vmax=v_max,
            shading="auto",
        )

        draw_boundary_overlays(ax, lon_is_0_360, biome_gdf=biome_gdf, basin_gdf=basin_gdf)
        ax.set_title(f"GWL: {target}°C\n(Actual: {actual_gwl:.2f}°C)", fontsize=11)
        ax.set_xlabel("Longitude [°E]", fontsize=9)

    axes[0].set_ylabel("Latitude [°N]", fontsize=10)

    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(mesh, cax=cbar_ax)
    cbar.set_label(meta.get("label", target_variable), fontsize=11)

    plt.subplots_adjust(right=0.90)
    plt.suptitle(
        f"Amazon ${target_variable}$ Spatial Distribution by Global Warming Level "
        f"({model_name}, {smoothing_title(rolling_years)})",
        fontsize=14,
        y=1.05,
    )

    output_dir = os.path.join(PREPROCESS_DIR, "plots", model_name)
    os.makedirs(output_dir, exist_ok=True)
    output_filename = f"SpatialMap_{target_variable}_{model_name}{suffix}.png"
    output_path = os.path.join(output_dir, output_filename)

    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()

    print(f"[Plot Success] Saved spatial map at: {output_path}")