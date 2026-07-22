# Amazon Basin Boundary

This directory is where the Amazon basin boundary shapefile goes. It is
**not included** in this repository -- the source repository has no
LICENSE file, so redistributing a copy here isn't clearly permitted.

## Where to get it

Download (or clone) the shapefile from:
<https://github.com/gamamo/AmazonBasinLimits>

## Expected files

Place the downloaded shapefile (and its sidecar files) here as:

```
data/basin/amazon_sensulatissimo_gmm_v1.shp
data/basin/amazon_sensulatissimo_gmm_v1.shx
data/basin/amazon_sensulatissimo_gmm_v1.dbf
data/basin/amazon_sensulatissimo_gmm_v1.prj
... (other .shp sidecar files as provided)
```

`preprocess.py` uses this only as a reference outline drawn on plots (see
`DEFAULT_BASIN_SHAPEFILE_PATH` / `load_boundary_shapefiles()` /
`draw_boundary_overlays()`) -- it is never used to mask or crop data.
