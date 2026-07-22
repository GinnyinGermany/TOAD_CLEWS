# CMIP6 Raw Data

This directory is where raw CMIP6 NetCDF files go if you want to run the
preprocessing pipeline from scratch (`python run_pipeline.py --stage all` or
`--stage preprocess-only`). Raw files are **not** included in this
repository -- download them yourself from the [ESGF Grid](https://esgf-metagrid.cloud.dkrz.de/search)
(or a mirror).

Most users don't need this: the repo already ships the cropped, GWL-indexed
`.nc` outputs of this step under `Amazon/processed_data/`, which is enough to
run the TOAD stage (`--stage toad-only`) and reproduce all downstream results.

## Expected structure

```
data/cmip6/{model}/{experiment}/{variable}_{frequency}_{model}_{experiment}_r1i1p1f1_*.nc
```

For example:

```
data/cmip6/GFDL-ESM4/1pctCO2/tas_Amon_GFDL-ESM4_1pctCO2_r1i1p1f1_gr1_000101-010012.nc
data/cmip6/GFDL-ESM4/piControl/tas_Amon_GFDL-ESM4_esm-piControl_r1i1p1f1_gr1_000101-010012.nc
```

## Variables and experiments needed (per `Amazon/configs/*.json`)

| Model | Variables (1pctCO2) | Variables (piControl) |
|---|---|---|
| GFDL-ESM4 | tas, pr, evspsbl, cVeg, fFire, grassFrac, treeFrac | tas |
| TaiESM1 | tas, pr, evspsbl, cVeg, fFire | tas |

`tas` under `piControl` is always required (used as the pre-industrial GWL
baseline); every other variable is only needed under `1pctCO2`.

## License & Citation

We acknowledge the World Climate Research Programme, which, through its
Working Group on Coupled Modelling, coordinated and promoted CMIP6. We thank
the climate modeling groups for producing and making available their model
output, the Earth System Grid Federation (ESGF) for archiving the data and
providing access, and the multiple funding agencies who support CMIP6 and
ESGF.

Both models below are distributed under a **Creative Commons
Attribution-ShareAlike 4.0 International License**
(<https://creativecommons.org/licenses/by-sa/4.0/>), per each file's own
`license` attribute:

- **GFDL-ESM4** -- NOAA-GFDL (National Oceanic and Atmospheric
  Administration, Geophysical Fluid Dynamics Laboratory, Princeton, NJ, USA)
  - `1pctCO2`: <https://furtherinfo.es-doc.org/CMIP6.NOAA-GFDL.GFDL-ESM4.1pctCO2.none.r1i1p1f1>
  - `piControl`: <https://furtherinfo.es-doc.org/CMIP6.NOAA-GFDL.GFDL-ESM4.piControl.none.r1i1p1f1>

- **TaiESM1** -- AS-RCEC (Research Center for Environmental Changes,
  Academia Sinica, Taipei, Taiwan)
  - `1pctCO2`: <https://furtherinfo.es-doc.org/CMIP6.AS-RCEC.TaiESM1.1pctCO2.none.r1i1p1f1>
  - `piControl`: <https://furtherinfo.es-doc.org/CMIP6.AS-RCEC.TaiESM1.piControl.none.r1i1p1f1>

Full CMIP6 terms of use (citation requirements, model-name conventions,
publication registration): <https://pcmdi.llnl.gov/CMIP6/TermsOfUse>.

The `.nc` files under `Amazon/processed_data/` are a spatially-cropped,
GWL-indexed derivative of this output, shared here under the same CC BY-SA
4.0 terms.
