# Interpolate freezing rain

```
mamba create -n map_digitization -c conda-forge python=3.10 gdal matplotlib pillow numpy scipy cartopy shapely xarray netcdf4
```

Create a PNG file of the map. For example:

Then add some control points -- x, y, lon, lat

```
gdal_translate -gcp 68 0 -122.775754 49.0 \
               -gcp 41 5 -124.677630 48.36 \
               -gcp 770 62 -66.982600 44.822 \
               -gcp 371 503 -97.203 25.969157 \
               -gcp 744 126 -70.222 42.062038 \
               -gcp 404 38 -95.139693 49.374 \
               -gcp 66 338 -117.102033 32.553 \
               -gcp 645 490 -81.059149 25.226 \
               -gcp 201 385 -108.206270 31.3331 \
               -gcp 685 243 -75.933826 37.1264 \
               -gcp 655 109 -76.346708 44.104 \
               -gcp 413 282 -94.618023 36.998 \
               clean_fzra_with_adds.png clean_fzra_with_adds.tif
```

now you have a georefenced tiff file called `clean_fzra_with_adds.tif`. You can check with this code:

```
python plot-gdal-pts.py clean_fzra_with_adds.tif 
```

if happy, regrid to a rectilinear grid

```
RESOLUTION=0.05
MINLON=-128.0
MAXLON=-65.0
MINLAT=24.0
MAXLAT=50.0
gdalwarp -t_srs 'EPSG:4326' -tr $RESOLUTION $RESOLUTION -te $MINLON $MINLAT $MAXLON $MAXLAT clean_fzra_with_adds.tif clean_fzra_with_adds_rect.tif
```

Now run the code:

```
python img_to_netcdf.py clean_fzra_with_adds_rect.tif
```
