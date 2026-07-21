"""
toad_runner.py

Wraps the TOAD shift-detection / clustering / plotting stage into reusable
functions so it can be driven programmatically for any (model, variable)
pair produced by preprocess.py.

Design note (method B): the JSON config is the single source of truth for
TOAD's shift/optimize settings, and this module re-reads it directly via
preprocess.load_pipeline_json() / preprocess.get_variable_config() -- the
same JSON file used at preprocessing time. The processed .nc file produced
by preprocess.py does NOT carry any "toad" settings in its attrs (see
TOADExporter in preprocess.py); it is purely the masked, GWL-indexed data.
This means:
  - Re-tuning shift_direction/shift_selection/optimize_params only requires
    editing the JSON and re-running the TOAD stage -- no need to regenerate
    the .nc file.
  - There is exactly one place (the JSON) where TOAD-run parameters live,
    so there's no attrs<->JSON round-trip/type-conversion layer to maintain
    and no risk of the two silently drifting apart.

Schema (current): config.json's "toad" block per variable now names both
algorithms explicitly under "shift":
  - "shift_method"   (e.g. "ASDETECT()")       -> used by td.compute_shifts
  - "cluster_method" (e.g. "SpaceTimeDBSCAN")  -> used by td.compute_clusters
Both are resolved dynamically via SHIFT_METHOD_REGISTRY /
CLUSTER_METHOD_REGISTRY below instead of being hardcoded, so adding a new
algorithm only means registering it here once.

"toad.optimize" now carries "optimize_objective" directly, and
"optimize_n_trials" is read from config.json too (build_run_config_from_json
checks both toad.optimize.optimize_n_trials and a toad-level sibling key, to
be robust to either placement). optimize_direction ("maximize"/"minimize")
is still not part of the schema, so it stays a module-level default
(DEFAULT_OPTIMIZE_DIRECTION) -- add it to the JSON and thread it through
build_run_config_from_json() if per-variable control is ever needed.

The old "optimize.optimize" boolean on/off toggle no longer exists in the
schema; optimization is treated as always-on.
"""

import ast
import io
import json
import logging
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

import preprocess

# ------------------------------------------------------------------
# Script location and project root (for dynamic path resolution)
# ------------------------------------------------------------------
TOAD_RUNNER_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(TOAD_RUNNER_DIR)

# ------------------------------------------------------------------
# TOAD package import
# ------------------------------------------------------------------
def _ensure_toad_on_path():
    """Add the toad-main directory (sibling of this project's parent dir) to
    sys.path, matching the original notebook's path setup."""
    toad_main_path = os.path.join(PROJECT_ROOT, "toad-main")
    if toad_main_path not in sys.path:
        sys.path.append(toad_main_path)


_ensure_toad_on_path()

from toad import TOAD  # noqa: E402
from toad.shifts import ASDETECT  # noqa: E402
from toad.clustering.methods.space_time_dbscan import SpaceTimeDBSCAN  # noqa: E402

import custom_shift_method  # noqa: E402


# ------------------------------------------------------------------
# Method registries: config.json now names the shift-detection algorithm
# and the clustering algorithm explicitly (toad.shift.shift_method /
# toad.shift.cluster_method) instead of hardcoding ASDETECT in this module.
# Add an entry here whenever a new TOAD method needs to be reachable from
# the JSON.
#
# Note the asymmetry (carried over from the original notebook): TOAD wants
# an *instance* for compute_shifts(method=...) but the bare *class* for
# compute_clusters(method=...). resolve_shift_method() instantiates;
# resolve_cluster_method() does not.
# ------------------------------------------------------------------
SHIFT_METHOD_REGISTRY = {
    "ASDETECT": ASDETECT,
}
CLUSTER_METHOD_REGISTRY = {
    "SpaceTimeDBSCAN": SpaceTimeDBSCAN,
}


def _clean_method_name(raw_name):
    """Normalize a JSON method string like 'ASDETECT()' -> 'ASDETECT'."""
    if raw_name is None:
        return None
    return raw_name.replace("(", "").replace(")", "").strip()


def resolve_shift_method(name, params=None):
    """Look up a shift-detection method by name and return an *instance*
    (TOAD's compute_shifts expects method=SomeClass()).

    Args:
        name: Shift method name (e.g., "FilteredASDETECT", "ASDETECT")
        params: Dictionary of customize parameters. For FilteredASDETECT,
                default values (lmin, lmax, segmentation) are merged with
                customize params from JSON.
    """
    if params is None:
        params = {}

    key = _clean_method_name(name)

    if key == "FilteredASDETECT":
        # FilteredASDETECT default values (applied to all variables)
        filtered_asdetect_defaults = {
            "lmin": 5,
            "lmax": None,
            "segmentation": "two_sided"
        }
        # Merge defaults with customize params from JSON
        # (JSON params override defaults)
        merged_params = {**filtered_asdetect_defaults, **params}
        return custom_shift_method.FilteredASDETECT(**merged_params)
    elif key == "ASDETECT":
        return ASDETECT(**params)
    else:
        raise KeyError(
            f"No shift method registered for '{name}'. "
            f"Available: FilteredASDETECT, ASDETECT"
        )


def resolve_cluster_method(name):
    """Look up a clustering method by name and return the bare *class*
    (TOAD's compute_clusters expects method=SomeClass, not an instance)."""
    key = _clean_method_name(name)
    method_cls = CLUSTER_METHOD_REGISTRY.get(key)
    if method_cls is None:
        raise KeyError(
            f"No cluster method registered for '{name}'. "
            f"Available: {list(CLUSTER_METHOD_REGISTRY.keys())}"
        )
    return method_cls


# ------------------------------------------------------------------
# Defaults for keys config.json still does not carry.
# config.json's "toad.optimize" block now DOES carry optimize_objective,
# and "toad" now carries optimize_n_trials -- see build_run_config_from_json.
# optimize_direction ("maximize"/"minimize") is still not in the schema, so
# it stays a module-level default. Likewise, the old boolean toggle
# "optimize.optimize" (on/off) was dropped from the schema entirely; we
# treat optimization as always-on to match the current workflow.
# ------------------------------------------------------------------
DEFAULT_OPTIMIZE_OBJECTIVE = "mean_spatial_autocorrelation"
DEFAULT_OPTIMIZE_DIRECTION = "maximize"
DEFAULT_OPTIMIZE_N_TRIALS = 1000

# Matches toad.utils.DEFAULT_SHIFT_THRESHOLD -- the minimum shift magnitude
# required for a point to be used in clustering / counted as a transition in
# time_of_max_shift_map. Kept here so both are driven by the same JSON value
# ("toad.shift.shift_threshold") instead of silently using TOAD's library
# default independently in each call site.
DEFAULT_SHIFT_THRESHOLD = 0.5

RESULTS_DIR = Path(os.path.join(TOAD_RUNNER_DIR, "results_toad"))


# ------------------------------------------------------------------
# Log parsing helpers (optimizer prints its best trial to the logger;
# TOAD does not return it directly, so we capture + regex-parse it)
# ------------------------------------------------------------------
def remove_ansi_codes(text):
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


def parse_best_optimization_from_log(log_text):
    log_text = remove_ansi_codes(log_text)
    result = {}

    best_pattern = (
        r"Best \(#(?P<best_trial>\d+)\): "
        r"score (?P<best_score>[-+0-9.eE]+), "
        r"params (?P<params>\{.*?\})"
    )
    best_match = re.search(best_pattern, log_text)
    if best_match is not None:
        params = ast.literal_eval(best_match.group("params"))
        result.update({
            "best_trial": int(best_match.group("best_trial")),
            "best_score": float(best_match.group("best_score")),
            "best_spatial_eps": params.get("spatial_eps"),
            "best_temporal_eps": params.get("temporal_eps"),
            "best_min_samples": params.get("min_samples"),
        })

    cluster_pattern = (
        r"New cluster variable\s+"
        r"(?P<cluster_variable>\S+): "
        r"Identified\s+"
        r"(?P<n_clusters>\d+)\s+clusters?\s+"
        r"in\s+(?P<n_clustered_points>\d+)\s+pts; "
        r"Left\s+(?P<noise_percent>[0-9.]+)%\s+as noise "
        r"\((?P<n_noise_points>\d+)\s+pts\)"
    )
    cluster_match = re.search(cluster_pattern, log_text)
    if cluster_match is not None:
        result.update({
            "cluster_variable": cluster_match.group("cluster_variable"),
            "n_clusters_from_log": int(cluster_match.group("n_clusters")),
            "n_clustered_points": int(cluster_match.group("n_clustered_points")),
            "noise_percent": float(cluster_match.group("noise_percent")),
            "n_noise_points": int(cluster_match.group("n_noise_points")),
        })

    return result


# ------------------------------------------------------------------
# Run-directory / labeling helpers
# ------------------------------------------------------------------
def smoothing_label(rolling_years):
    if rolling_years is None or rolling_years <= 1:
        return "annual"
    return f"smooth{rolling_years}yr"


def safe_float_label(value):
    return str(value).replace(".", "p")


def optimize_range_label(optimize_params):
    sp = optimize_params["spatial_eps"]
    te = optimize_params["temporal_eps"]
    ms = optimize_params["min_samples"]
    return (
        f"sp{int(sp[0])}-{int(sp[1])}_"
        f"t{safe_float_label(te[0])}-{safe_float_label(te[1])}_"
        f"ms{int(ms[0])}-{int(ms[1])}"
    )


def make_run_dir(config):
    """Create shift/ and cluster/ directories with new structure.

    For optimization runs, creates timestamped subdirectories to preserve
    multiple optimization experiments with different params/objectives.

    Returns:
        (shift_dir, cluster_dir, shift_run_name, cluster_run_name)
    """
    from datetime import datetime

    model_var_path = RESULTS_DIR / config['model'] / config['variable']

    # Shift directory (same for all clustering experiments)
    shift_run_name = (
        f"{config['model']}/"
        f"{config['variable']}/"
        f"shift"
    )
    shift_dir = model_var_path / "shift"
    shift_dir.mkdir(parents=True, exist_ok=True)

    # Cluster directory with parameters or optimization
    if config['optimize']:
        # Add timestamp to distinguish multiple optimization runs
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        cluster_subdir = f"optimization_{timestamp}"
    else:
        # Build cluster folder name from clustering_params
        spatial_eps = config['clustering_params'].get('spatial_eps', '?')
        temporal_eps = config['clustering_params'].get('temporal_eps', '?')
        min_samples = config['clustering_params'].get('min_samples', '?')
        cluster_subdir = f"eps{spatial_eps}_temp{temporal_eps}_min{min_samples}"

    cluster_run_name = (
        f"{config['model']}/"
        f"{config['variable']}/"
        f"cluster/{cluster_subdir}"
    )
    cluster_dir = model_var_path / "cluster" / cluster_subdir
    cluster_dir.mkdir(parents=True, exist_ok=True)

    return shift_dir, cluster_dir, shift_run_name, cluster_run_name


def flatten_dict(d, prefix=""):
    out = {}
    if not isinstance(d, dict):
        return out
    for key, value in d.items():
        new_key = f"{prefix}_{key}" if prefix else key
        if isinstance(value, dict):
            out.update(flatten_dict(value, new_key))
        else:
            out[new_key] = value
    return out


# ------------------------------------------------------------------
# Saving cluster results
# ------------------------------------------------------------------
def save_optimized_toad_run_light(td, config, opt_result=None, best_optimization=None):
    run_dir, run_name = make_run_dir(config)

    metadata = config.copy()
    metadata["run_name"] = run_name

    if best_optimization is not None:
        metadata.update(best_optimization)
    if opt_result is not None:
        metadata["opt_result_repr"] = repr(opt_result)

    try:
        cluster_variable = metadata.get("cluster_variable", config["variable"])
        cluster_ids = [
            int(cid) for cid in td.get_cluster_ids(cluster_variable) if int(cid) != -1
        ]
    except Exception as e:
        cluster_ids = []
        metadata["cluster_id_error"] = str(e)

    metadata["n_clusters"] = len(cluster_ids)
    metadata["cluster_ids"] = cluster_ids

    with open(run_dir / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=4, default=str)

    rows = []
    if len(cluster_ids) == 0:
        rows.append({**metadata, "cluster_id": None})
    else:
        for cluster_id in cluster_ids:
            row = {**metadata, "cluster_id": cluster_id}
            try:
                time_stats = td.stats(cluster_variable).time.all_stats(cluster_id)
                row.update(flatten_dict(time_stats, prefix="time"))
            except Exception as e:
                row["time_stats_error"] = str(e)
            try:
                general_stats = td.stats(cluster_variable).general.all_stats(cluster_id)
                row.update(flatten_dict(general_stats, prefix="general"))
            except Exception as e:
                row["general_stats_error"] = str(e)
            rows.append(row)

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(run_dir / "cluster_stats.csv", index=False)

    print(f"[Saved] {run_dir}")
    print(f"[Saved] {run_dir / 'metadata.json'}")
    print(f"[Saved] {run_dir / 'cluster_stats.csv'}")

    return run_dir, stats_df


# ------------------------------------------------------------------
# Saving shift and cluster results (separated)
# ------------------------------------------------------------------
def save_shift_results(td, shift_dir, config):
    """Save shift detection results to shift/ directory."""
    shift_dir = Path(shift_dir)

    # Create metadata for shift
    shift_metadata = {
        "model": config["model"],
        "variable": config["variable"],
        "rolling_years": config["rolling_years"],
        "shift_method": config["shift_method_name"],
        "shift_params": config["shift_params"],
        "shift_direction": config["shift_direction"],
        "shift_selection": config["shift_selection"],
        "shift_threshold": config["shift_threshold"],
    }

    with open(shift_dir / "config_snapshot.json", "w") as f:
        json.dump(shift_metadata, f, indent=4, default=str)

    print(f"[Saved] {shift_dir / 'config_snapshot.json'}")

    # Save shift plots
    plot_dir = shift_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # 1. Maximum shift magnitude map
    td.plot.max_shift_map()
    plt.savefig(plot_dir / "01_max_shift_map.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 2. Time of maximum shift map
    td.plot.time_of_max_shift_map(shift_threshold=config["shift_threshold"])
    plt.savefig(plot_dir / "02_time_of_max_shift_map.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[Saved plots] {plot_dir}")


def save_cluster_results(td, cluster_dir, config, opt_result=None, best_optimization=None, cluster_variable=None, max_clusters=5):
    """Save clustering results to cluster/{params}/ directory."""
    from datetime import datetime

    cluster_dir = Path(cluster_dir)

    if cluster_variable is None:
        cluster_variable = config["variable"]

    # Create metadata for cluster
    cluster_metadata = {
        "model": config["model"],
        "variable": config["variable"],
        "rolling_years": config["rolling_years"],
        "cluster_method": config["cluster_method_name"],
        "optimize": config["optimize"],
        "timestamp": datetime.now().isoformat(),
    }

    if config["optimize"]:
        cluster_metadata.update({
            "optimize_objective": config["optimize_objective"],
            "optimize_params": config["optimize_params"],
            "optimize_n_trials": config["optimize_n_trials"],
        })
    else:
        cluster_metadata["clustering_params"] = config["clustering_params"]

    if best_optimization is not None:
        cluster_metadata.update(best_optimization)
    if opt_result is not None:
        cluster_metadata["opt_result_repr"] = repr(opt_result)

    with open(cluster_dir / "config_snapshot.json", "w") as f:
        json.dump(cluster_metadata, f, indent=4, default=str)

    # Get cluster IDs and stats
    try:
        cluster_ids = [
            int(cid) for cid in td.get_cluster_ids(cluster_variable) if int(cid) != -1
        ]
    except Exception as e:
        cluster_ids = []
        cluster_metadata["cluster_id_error"] = str(e)

    cluster_metadata["n_clusters"] = len(cluster_ids)
    cluster_metadata["cluster_ids"] = cluster_ids

    # Save cluster stats
    rows = []
    if len(cluster_ids) == 0:
        rows.append({**cluster_metadata, "cluster_id": None})
    else:
        for cluster_id in cluster_ids:
            row = {**cluster_metadata, "cluster_id": cluster_id}
            try:
                time_stats = td.stats(cluster_variable).time.all_stats(cluster_id)
                row.update(flatten_dict(time_stats, prefix="time"))
            except Exception as e:
                row["time_stats_error"] = str(e)
            try:
                general_stats = td.stats(cluster_variable).general.all_stats(cluster_id)
                row.update(flatten_dict(general_stats, prefix="general"))
            except Exception as e:
                row["general_stats_error"] = str(e)
            rows.append(row)

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(cluster_dir / "cluster_stats.csv", index=False)

    with open(cluster_dir / "config_snapshot.json", "w") as f:
        json.dump(cluster_metadata, f, indent=4, default=str)

    print(f"[Saved] {cluster_dir / 'config_snapshot.json'}")
    print(f"[Saved] {cluster_dir / 'cluster_stats.csv'}")

    # Save cluster plots
    plot_dir = cluster_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    selected_cluster_ids = cluster_ids[:max_clusters]

    # 3. Timeseries - one plot per cluster
    if len(selected_cluster_ids) > 0:
        for cluster_id in selected_cluster_ids:
            td.plot.timeseries(var=cluster_variable, cluster_ids=[cluster_id])
            plt.savefig(plot_dir / f"03_timeseries_cluster_{cluster_id}.png", dpi=300, bbox_inches="tight")
            plt.close()
    else:
        print("[Warning] No valid cluster IDs found for timeseries plot.")

    # 4. Overview - global
    td.plot.overview(var=cluster_variable, map_style={"projection": "global"})
    plt.savefig(plot_dir / "04_overview_global.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 5. Overview - aggregated
    td.plot.overview(var=cluster_variable, mode="aggregated")
    plt.savefig(plot_dir / "05_overview_aggregated.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[Saved plots] {plot_dir}")

    return stats_df


# ------------------------------------------------------------------
# Saving plots (legacy)
# ------------------------------------------------------------------
def save_toad_plots(td, run_dir, cluster_variable, max_clusters=5):
    plot_dir = Path(run_dir) / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    print("Using cluster variable:", cluster_variable)

    cluster_ids = [
        int(cid) for cid in td.get_cluster_ids(cluster_variable) if int(cid) != -1
    ]
    selected_cluster_ids = cluster_ids[:max_clusters]
    print("Cluster IDs:", selected_cluster_ids)

    # 1. Maximum shift magnitude map
    td.plot.max_shift_map()
    plt.savefig(plot_dir / "01_max_shift_map.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 2. Time of maximum shift map
    td.plot.time_of_max_shift_map()
    plt.savefig(plot_dir / "02_time_of_max_shift_map.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 3. Timeseries - one plot per cluster
    if len(selected_cluster_ids) > 0:
        for cluster_id in selected_cluster_ids:
            td.plot.timeseries(var=cluster_variable, cluster_ids=[cluster_id])
            plt.savefig(plot_dir / f"03_timeseries_cluster_{cluster_id}.png", dpi=300, bbox_inches="tight")
            plt.close()
    else:
        print("[Warning] No valid cluster IDs found for timeseries plot.")

    # 4. Overview - global
    td.plot.overview(var=cluster_variable, map_style={"projection": "global"})
    plt.savefig(plot_dir / "04_overview_global.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 5. Overview - aggregated
    td.plot.overview(var=cluster_variable, mode="aggregated")
    plt.savefig(plot_dir / "05_overview_aggregated.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[Saved plots] {plot_dir}")
    return plot_dir


# ------------------------------------------------------------------
# JSON config -> RUN_CONFIG
# ------------------------------------------------------------------
def build_run_config_from_json(
    config_path,
    model,
    variable,
    rolling_years,
    optimize_objective=None,
    optimize_direction=None,
    optimize_n_trials=None,
):
    """Build the RUN_CONFIG dict by reading the variable's "toad" block
    (shift + optimize) straight from the pipeline JSON -- the same file
    preprocess.py used to generate the .nc file. Nothing is read from the
    .nc file itself; the JSON is the single source of truth for these
    settings.

    optimize_objective / optimize_direction / optimize_n_trials passed here
    are treated as explicit overrides; when left as None (the default),
    the value is taken from the JSON, falling back to the module-level
    DEFAULT_* constants only if the JSON is also missing it.

    Schema notes:
      - toad.shift now carries "shift_method" (e.g. "ASDETECT()") and
        "cluster_method" (e.g. "SpaceTimeDBSCAN") as separate keys, resolved
        below via resolve_shift_method()/resolve_cluster_method().
      - toad.optimize now carries "optimize_objective" directly.
      - "optimize_n_trials" is read from toad.optimize first, falling back
        to a sibling key directly under "toad" -- config.json currently
        places it as a sibling of "shift"/"optimize" rather than nested
        inside "optimize", so both locations are checked to be robust to
        either layout.
      - The old "optimize.optimize" boolean on/off toggle no longer exists
        in the schema; optimization is treated as always-on.
    """
    pipeline_json = preprocess.load_pipeline_json(config_path)
    var_cfg = preprocess.get_variable_config(pipeline_json, variable)
    toad_cfg = var_cfg.get("toad", {})

    shift_cfg = toad_cfg.get("shift", {})
    cluster_cfg = toad_cfg.get("cluster", {})

    # Clustering configuration
    optimize_enabled = cluster_cfg.get("optimize", True)

    if optimize_enabled:
        # optimize=true: use optimization params (ranges)
        optimize_params_raw = cluster_cfg.get("optimize_params", {})
        optimize_params = {k: tuple(v) for k, v in optimize_params_raw.items()}
        resolved_objective = (
            optimize_objective
            or cluster_cfg.get("optimize_objective")
            or DEFAULT_OPTIMIZE_OBJECTIVE
        )
        resolved_n_trials = (
            optimize_n_trials
            or cluster_cfg.get("optimize_n_trials")
            or toad_cfg.get("optimize_n_trials")
            or DEFAULT_OPTIMIZE_N_TRIALS
        )
        clustering_params = {}
    else:
        # optimize=false: use fixed clustering params
        optimize_params = {}
        resolved_objective = None
        resolved_n_trials = None
        clustering_params = cluster_cfg.get("clustering_params", {})

    shift_method_name = shift_cfg.get("shift_method", "ASDETECT")
    shift_params = shift_cfg.get("shift_params", {})
    cluster_method_name = cluster_cfg.get("cluster_method", "SpaceTimeDBSCAN")
    shift_threshold = shift_cfg.get("shift_threshold", DEFAULT_SHIFT_THRESHOLD)

    resolved_direction = optimize_direction or DEFAULT_OPTIMIZE_DIRECTION

    return {
        "model": model,
        "variable": variable,
        "rolling_years": rolling_years,
        "shift_method_name": shift_method_name,
        "shift_params": shift_params,
        "cluster_method_name": cluster_method_name,
        "method": cluster_method_name,
        "shift_direction": shift_cfg.get("shift_direction", "negative"),
        "shift_selection": shift_cfg.get("shift_selection", "local"),
        "shift_threshold": shift_threshold,
        "optimize": optimize_enabled,
        "optimize_params": optimize_params,
        "optimize_objective": resolved_objective,
        "optimize_direction": resolved_direction,
        "optimize_n_trials": resolved_n_trials,
        "clustering_params": clustering_params,
    }


# ------------------------------------------------------------------
# Public entry point
# ------------------------------------------------------------------
def run_toad_for_variable(
    nc_path,
    model,
    variable,
    rolling_years,
    config_path,
    max_clusters=5,
    optimize_objective=None,
    optimize_direction=None,
    optimize_n_trials=None,
):
    """Runs the full TOAD stage for a single processed .nc file:
    shift detection -> (optimized) clustering -> save result tables -> save plots.

    config_path points at the same pipeline JSON used to preprocess nc_path;
    the "toad" (shift/optimize) settings, including which shift-detection
    and clustering algorithms to use, are read from it directly rather than
    from the .nc file's attrs.

    Returns (run_dir, stats_df).
    """
    ds = xr.open_dataset(nc_path)
    td = TOAD(ds, time_dim="GWL")

    run_config = build_run_config_from_json(
        config_path, model, variable, rolling_years,
        optimize_objective=optimize_objective,
        optimize_direction=optimize_direction,
        optimize_n_trials=optimize_n_trials,
    )

    shift_method_instance = resolve_shift_method(
        run_config["shift_method_name"],
        params=run_config["shift_params"]
    )
    cluster_method_cls = resolve_cluster_method(run_config["cluster_method_name"])

    # Step 1: shift detection -- method now comes from toad.shift.shift_method
    td.drop_shifts()
    td.compute_shifts(var=variable, method=shift_method_instance)

    # Step 2: clustering, with log capture so we can recover the optimizer's
    # best trial (TOAD only prints it, it isn't returned directly).
    log_stream = io.StringIO()
    handler = logging.StreamHandler(log_stream)
    handler.setLevel(logging.INFO)

    logger = td.logger
    old_level = logger.level
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    try:
        # Build clustering arguments
        cluster_kwargs = {
            "var": variable,
            "method": cluster_method_cls,
            "shift_direction": run_config["shift_direction"],
            "shift_selection": run_config["shift_selection"],
            "shift_threshold": run_config["shift_threshold"],
            "optimize": run_config["optimize"],
        }

        if run_config["optimize"]:
            # Optimization mode: use optimization params (ranges)
            cluster_kwargs.update({
                "optimize_params": run_config["optimize_params"],
                "optimize_objective": run_config["optimize_objective"],
                "optimize_direction": run_config["optimize_direction"],
                "optimize_n_trials": run_config["optimize_n_trials"],
            })
        else:
            # Single clustering mode: compute_clusters() has no top-level
            # spatial_eps/temporal_eps/min_samples kwargs -- when method is
            # passed as a bare class (not optimizing), TOAD instantiates it
            # with NO arguments internally (`method()`), which would raise
            # a TypeError for SpaceTimeDBSCAN (spatial_eps/temporal_eps are
            # required). Fixed params must instead be applied by
            # pre-instantiating the method here.
            cluster_kwargs["method"] = cluster_method_cls(**run_config["clustering_params"])

        opt_result = td.compute_clusters(**cluster_kwargs)
    finally:
        logger.removeHandler(handler)
        logger.setLevel(old_level)

    log_text = log_stream.getvalue()
    best_optimization = parse_best_optimization_from_log(log_text)
    print(best_optimization)

    # Step 3: Create shift/ and cluster/ directories
    shift_dir, cluster_dir, shift_run_name, cluster_run_name = make_run_dir(run_config)

    # Step 4: Save shift detection results
    save_shift_results(
        td=td,
        shift_dir=shift_dir,
        config=run_config,
    )

    # Step 5: Save clustering results and plots
    stats_df = save_cluster_results(
        td=td,
        cluster_dir=cluster_dir,
        config=run_config,
        opt_result=opt_result,
        best_optimization=best_optimization,
        cluster_variable=best_optimization.get("cluster_variable", variable),
        max_clusters=max_clusters,
    )

    return cluster_dir, stats_df