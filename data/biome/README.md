# Amazon Biome Boundary

This directory is where the Amazon biome boundary shapefile goes. It is
**not included** in this repository -- RAISG's terms of use restrict the
data to non-profit use, require it to be redistributed unaltered, and are
silent on redistribution through third-party repositories, so we don't
vendor a copy here.

## Where to get it

Download the biome shapefile from RAISG's map portal:
<https://www.raisg.org/en/maps/> (registration required). Look for the
"Biomas" / biome layer.

## Expected files

Place the downloaded shapefile (and its sidecar files) here as:

```
data/biome/Biomas.shp
data/biome/Biomas.shx
data/biome/Biomas.dbf
data/biome/Biomas.prj
... (other .shp sidecar files as provided)
```

`preprocess.py` filters this layer to the "Amazonía" biome via the
`bioma` attribute column (see `DEFAULT_BIOME_SHAPEFILE_PATH` /
`compute_buffered_bbox()` / `load_boundary_shapefiles()`).

## Citation

Cite this dataset using the citation format given in RAISG's metadata for
the biome layer (see the "citation" field in the dataset's metadata on the
RAISG portal). RAISG's terms of use also require attributing the original
national data sources listed in that metadata.
