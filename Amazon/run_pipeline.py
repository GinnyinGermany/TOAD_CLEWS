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

import glob
import json
import os
import traceback

import preprocess
import toad_runner

CONFIG_DIR = "./configs"


def discover_configs(config_dir=CONFIG_DIR):
    """Scan config_dir for model JSON files.

    Returns a list of (config_path, model_name, [variable_name, ...]) tuples,
    one per JSON file.
    """
    entries = []
    for config_path in sorted(glob.glob(os.path.join(config_dir, "*.json"))):
        with open(config_path, "r") as f:
            pipeline_json = json.load(f)
        model_name = pipeline_json["model"]
        variables = list(pipeline_json.get("variable_settings", {}).keys())
        entries.append((config_path, model_name, variables))
    return entries


def run_preprocessing_stage(config_path, model_name, variable):
    """Runs preprocess.py's generate_nc + both plot functions for one variable.

    Returns (nc_path, rolling_years) so the TOAD stage can reuse the exact
    rolling_years value that was actually applied.
    """
    pipeline_json = preprocess.load_pipeline_json(config_path)
    var_cfg = preprocess.get_variable_config(pipeline_json, variable)
    rolling_years = var_cfg.get("preprocessing", {}).get("rolling_years")

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
    configs = discover_configs()
    if not configs:
        print(f"[Error] No JSON config files found in {CONFIG_DIR}/")
        return

    successes = []
    failures = []

    for config_path, model_name, variables in configs:
        print(f"\n{'=' * 60}\nModel: {model_name}  (config: {config_path})\n{'=' * 60}")

        # Model year <-> GWL plot is variable-independent, so it's generated
        # once per model rather than once per (model, variable) pair.
        try:
            preprocess.plot_gwl_by_model_year(model_name, config_path=config_path)
        except Exception:
            print(f"[Warning] plot_gwl_by_model_year failed for {model_name}")
            traceback.print_exc()

        for variable in variables:
            label = f"{model_name}/{variable}"
            try:
                print(f"\n--- [{label}] Preprocessing ---")
                nc_path, rolling_years = run_preprocessing_stage(
                    config_path, model_name, variable
                )

                print(f"\n--- [{label}] TOAD ---")
                run_dir, stats_df = run_toad_stage(
                    nc_path, model_name, variable, rolling_years, config_path
                )

                successes.append(
                    {"model": model_name, "variable": variable, "run_dir": str(run_dir)}
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