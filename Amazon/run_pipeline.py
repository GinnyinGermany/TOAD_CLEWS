"""
run_pipeline.py

Top-level orchestrator for the full CMIP6 -> preprocess -> TOAD workflow.

For every model config JSON found in configs/, and for every variable listed
in that JSON's "variable_settings":

  1. preprocess.py stage
     - generate_nc(): raw CMIP6 -> masked, GWL-indexed, TOAD-ready .nc file
     - plot_trend_from_nc(): area-averaged trend vs. GWL plot
     - plot_spatial_maps_from_nc(): 1x6 spatial map grid across GWL levels
     - plot_gwl_by_model_year(): once per model (variable-independent)

  2. toad_runner.py stage
     - opens the .nc file produced above
     - runs TOAD shift detection (ASDETECT) and (optimized) clustering
       (SpaceTimeDBSCAN), reading shift_direction/shift_selection/
       optimize_params directly from the same JSON config (config_path is
       passed straight through) -- the .nc file itself carries no TOAD
       parameters, so re-tuning them never requires regenerating the .nc
     - saves TOAD's cluster metadata/stats table and diagnostic plots under
       results_toad/{model}/{variable}/...

Run with:  python run_pipeline.py

Adding a new model or variable requires no changes to this file: drop a new
JSON into configs/, or add an entry to an existing JSON's variable_settings.
"""

import argparse
import glob
import json
import os
import traceback

import preprocess
import toad_runner

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_DIR = os.path.join(SCRIPT_DIR, "configs")
PROCESSED_DATA_DIR = os.path.join(SCRIPT_DIR, "processed_data")


def discover_configs(config_dir=CONFIG_DIR, filter_model=None):
    """Scan config_dir for model JSON files.

    Returns a list of (config_path, model_name, [variable_name, ...]) tuples,
    one per JSON file.

    Args:
        config_dir: Directory containing config JSON files
        filter_model: If specified, only return configs for this model
    """
    entries = []
    for config_path in sorted(glob.glob(os.path.join(config_dir, "*.json"))):
        with open(config_path, "r") as f:
            pipeline_json = json.load(f)
        model_name = pipeline_json["model"]
        variables = list(pipeline_json.get("variable_settings", {}).keys())

        # Filter by model if specified
        if filter_model and model_name != filter_model:
            continue

        entries.append((config_path, model_name, variables))
    return entries


def get_nc_path(model_name, variable, rolling_years=None):
    """Construct the expected NC file path for a given model/variable/rolling_years."""
    from preprocess import smoothing_suffix
    suffix = smoothing_suffix(rolling_years)
    nc_filename = f"TOAD_{variable}_{model_name}{suffix}.nc"
    return os.path.join(PROCESSED_DATA_DIR, model_name, nc_filename)


def run_preprocessing_stage(config_path, model_name, variable, skip_if_exists=False):
    """Runs preprocess.py's generate_nc + both plot functions for one variable.

    Returns (nc_path, rolling_years) so the TOAD stage can reuse the exact
    rolling_years value that was actually applied.

    Args:
        skip_if_exists: If True and NC file exists, skip preprocessing
    """
    pipeline_json = preprocess.load_pipeline_json(config_path)
    var_cfg = preprocess.get_variable_config(pipeline_json, variable)
    rolling_years = var_cfg.get("preprocessing", {}).get("rolling_years")

    nc_path = get_nc_path(model_name, variable, rolling_years)

    # Check if NC file already exists
    if skip_if_exists and os.path.exists(nc_path):
        print(f"[Cache] NC file already exists: {nc_path}")
        return nc_path, rolling_years

    nc_path = preprocess.generate_nc(
        model_name,
        variable,
        rolling_years=rolling_years,
        config_path=config_path,
    )

    preprocess.plot_trend_from_nc(
        model_name,
        variable,
        rolling_years=rolling_years,
        config_path=config_path,
    )

    preprocess.plot_spatial_maps_from_nc(
        model_name,
        variable,
        rolling_years=rolling_years,
        config_path=config_path,
    )

    return nc_path, rolling_years


def run_toad_stage(nc_path, model_name, variable, rolling_years, config_path):
    """Feeds a processed .nc file into TOAD and saves its cluster results + plots.

    config_path is passed through so toad_runner reads the same JSON's
    "toad" (shift/optimize) settings directly -- the .nc file itself carries
    no TOAD-run parameters.

    Returns (run_dir, stats_df) from toad_runner.
    """
    return toad_runner.run_toad_for_variable(
        nc_path=nc_path,
        model=model_name,
        variable=variable,
        rolling_years=rolling_years,
        config_path=config_path,
    )


def main():
    parser = argparse.ArgumentParser(
        description="TOAD Pipeline: Preprocessing + Shift Detection + Clustering Optimization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run everything for all models and variables
  python run_pipeline.py

  # Process only GFDL-ESM4 model
  python run_pipeline.py --model GFDL-ESM4

  # Process only cVeg variable of GFDL-ESM4
  python run_pipeline.py --model GFDL-ESM4 --variable cVeg

  # Re-run TOAD only (skip preprocessing if NC exists)
  python run_pipeline.py --model GFDL-ESM4 --variable cVeg --stage toad-only

  # Re-run preprocessing only
  python run_pipeline.py --model GFDL-ESM4 --variable cVeg --stage preprocess-only
        """
    )
    parser.add_argument(
        "--model",
        help="Process only this model (e.g., GFDL-ESM4). If not specified, processes all models."
    )
    parser.add_argument(
        "--variable",
        help="Process only this variable (e.g., cVeg). Requires --model. If not specified, processes all variables."
    )
    parser.add_argument(
        "--stage",
        choices=["all", "preprocess-only", "toad-only"],
        default="all",
        help="Which stage(s) to run. 'toad-only' skips preprocessing if NC file exists."
    )

    args = parser.parse_args()

    # Validation
    if args.variable and not args.model:
        parser.error("--variable requires --model to be specified")

    # Determine what to run
    run_preprocess = args.stage in ["all", "preprocess-only"]
    run_toad = args.stage in ["all", "toad-only"]
    skip_preprocess_if_exists = args.stage == "toad-only"

    # Print configuration
    print(f"\n{'=' * 60}\nPipeline Configuration\n{'=' * 60}")
    print(f"Stage: {args.stage}")
    if args.model:
        print(f"Model filter: {args.model}")
    if args.variable:
        print(f"Variable filter: {args.variable}")
    print(f"Will run preprocessing: {run_preprocess}")
    print(f"Will run TOAD: {run_toad}")
    print(f"Skip preprocessing if NC exists: {skip_preprocess_if_exists}")
    print()

    configs = discover_configs(filter_model=args.model)
    if not configs:
        print(f"[Error] No JSON config files found in {CONFIG_DIR}/")
        if args.model:
            print(f"        (filtered for model: {args.model})")
        return

    successes = []
    failures = []

    for config_path, model_name, all_variables in configs:
        print(f"\n{'=' * 60}\nModel: {model_name}  (config: {config_path})\n{'=' * 60}")

        # Filter variables if specified
        variables = [all_variables] if args.variable else all_variables
        if args.variable:
            if args.variable not in all_variables:
                print(f"[Error] Variable '{args.variable}' not found in model config")
                print(f"        Available: {all_variables}")
                continue
            variables = [args.variable]

        # Model year <-> GWL plot is variable-independent
        if run_preprocess:
            try:
                preprocess.plot_gwl_by_model_year(model_name, config_path=config_path)
            except Exception:
                print(f"[Warning] plot_gwl_by_model_year failed for {model_name}")
                traceback.print_exc()

        for variable in variables:
            label = f"{model_name}/{variable}"
            try:
                if run_preprocess:
                    print(f"\n--- [{label}] Preprocessing ---")
                    nc_path, rolling_years = run_preprocessing_stage(
                        config_path, model_name, variable,
                        skip_if_exists=skip_preprocess_if_exists
                    )
                else:
                    # Load rolling_years from config for TOAD-only runs
                    pipeline_json = preprocess.load_pipeline_json(config_path)
                    var_cfg = preprocess.get_variable_config(pipeline_json, variable)
                    rolling_years = var_cfg.get("preprocessing", {}).get("rolling_years")
                    nc_path = get_nc_path(model_name, variable, rolling_years)

                    if not os.path.exists(nc_path):
                        raise FileNotFoundError(
                            f"NC file not found: {nc_path}\n"
                            f"Please run preprocessing first or use --stage all"
                        )

                if run_toad:
                    print(f"\n--- [{label}] TOAD ---")
                    run_dir, stats_df = run_toad_stage(
                        nc_path, model_name, variable, rolling_years, config_path
                    )
                else:
                    run_dir = None

                successes.append(
                    {"model": model_name, "variable": variable, "run_dir": str(run_dir) if run_dir else "skipped"}
                )

            except Exception as e:
                print(f"[Failed] {label}: {e}")
                traceback.print_exc()
                failures.append(
                    {"model": model_name, "variable": variable, "error": str(e)}
                )

    print(f"\n{'=' * 60}\nPipeline summary\n{'=' * 60}")
    print(f"Succeeded: {len(successes)}")
    for s in successes:
        print(f"  [OK]   {s['model']}/{s['variable']} -> {s['run_dir']}")

    print(f"Failed: {len(failures)}")
    for f in failures:
        print(f"  [FAIL] {f['model']}/{f['variable']}: {f['error']}")


if __name__ == "__main__":
    main()