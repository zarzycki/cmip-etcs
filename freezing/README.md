# Precipitation typing for climate data

Prerequisites:

- numpy
- xarray
- metpy
- tqdm (this is used for timing loops and is optional scientifically, although code will need to be cleaned of its calls)
- calpreciptype (**See below!**)

The first four can be installed with `conda install` or `pip`. The last one comes with this repo.

## Workflow

This code leans quite heavily on both xarray and numpy, with metpy used to do some thermodynamic calculations. I could probably generalize a bit more beyond xarray but eh.

1. Load **T**, **Q**, and **PS** data.
	- Q can be RH, etc., but must be converted to specific humidity at this point.
	- Assuming T is on a higher resolution vertical grid, Q can be interpolated to the same vertical coordinate.
2. Calculate **pmid** and **pint**.
	- If the data is on *constant pressure* surfaces, use the `constant_press_prep` function to process.
	- If the data is on *model* surfaces, pmid and pint need to be defined. This is probably via the use of a hybrid calculation (e.g., P(z,i,j) = a(k)\*P0 + b(k)\*PS(i,j)) but for other models a 4-D (ntime, nlev, nlat, nlon) pressure array may be provided.
3. Run `calc_ptype` wrapper with **T**, **Q**, **pmid**, **pint** as inputs (plus some integers describing array sizes).

Some important notes:

- One should be able to reproduce ptypes from dev by checking out the repo and setting `dataset = "ERA5"` (constant pressure surfaces) or `dataset = "LENS"` (model surfaces).
- All arrays should be ordered **top-to-bottom**! These can be flipped in step 1 by using the ::-1 convention.
- For *constant pressure* surfaces, it is unclear what resolution is needed for accurate typing. ERA5 has been tested with deltaP of 25mb in the lower troposphere, 50mb in the mid-troposphere, and 100mb in the upper troposphere to stratosphere with decent results (i.e., freezing rain shows up spatially where historical events occurred). Gut tells me we need >10 levels in the troposphere, but we may want to do some stride tests with reference data!
- For *constant pressure* surfaces, if we have high resolution **T** (vertically) but low resolution **Q**, it seems reasonable to interpolate Q vertically to T given that the algos are more sensitive to the former variable *and* **Q** can be interpolated using a log formulation given C-C scaling. 

---

### calc_ptype

This function uses three different algorithms from NCEP and chooses either A) the dominant ptype where 2 or 3 algos agree. If *no* algorithms agree, chooses Bourgouin for consistency with Zarzycki (2018).

The inputs to `calc_ptype` are below... one note, setting zint = -1 (negative integer) will cause the code to call `calc_zint` which uses the hypsometric relationship to build interface heights iteratively from the surface. Technically, one could use an actual interface height variable from a model, but for comparison purposes, it may make sense to just internally calculate for everything to be apples-to-appples.

Reminder, **top-to-bottom**!!

```
# T (ntim,nlev,nlat,nlon) in K
# Q (ntim,nlev,nlat,nlon) in kg/kg
# pmid (ntim,nlev,nlat,nlon) in Pa
# pint (ntim,nlev+1,nlat,nlon) in Pa
# zint (ntim,nlev+1,nlat,nlon) in m
# ntim (int)
# nlev (int)
# nlon (int)
# nlat (int)
```

### calpreciptype

`calpreciptype.f90` is Fortran, so we need to use numpy to create wrappers/compile for use with Python. **This only needs to be done once**, although it may be good practice to recompile this every once in a while to account for compiler changes, etc. This can be done by making sure `numpy` is in the current environment and running:

```
python -m numpy.f2py --opt='-O3' -c calpreciptype.f90 -m calpreciptype
```

This should create a file that has a name *such as* `calpreciptype.cpython-38-darwin.so` (the actual name will differ based on system).

Functions in the Fortran can now be called, such as...

```
ptype = calpreciptype.calwxt_revised(...)
```

Interfaces can be printed using:

```
print(calpreciptype.calwxt_ramer.__doc__)
print(calpreciptype.calwxt_bourg.__doc__)
print(calpreciptype.calwxt_revised.__doc__)
print(calpreciptype.calwxt.__doc__)
```

### Sample data

```
YYYY=1998
MM=01
DD=06
STRIDE=4
ncks -O -d time,"${YYYY}-${MM}-${DD} 00:00:00","${YYYY}-${MM}-${DD} 23:00:00",6 -d latitude,0,,${STRIDE} -d longitude,0,,${STRIDE} /glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/${YYYY}${MM}/e5.oper.an.sfc.128_134_sp.ll025sc.${YYYY}${MM}0100_${YYYY}${MM}3123.nc PS.nc
ncks -O -d time,"${YYYY}-${MM}-${DD} 00:00:00","${YYYY}-${MM}-${DD} 23:00:00",6 -d latitude,0,,${STRIDE} -d longitude,0,,${STRIDE} /glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.pl/${YYYY}${MM}/e5.oper.an.pl.128_133_q.ll025sc.${YYYY}${MM}${DD}00_${YYYY}${MM}${DD}23.nc Q.nc
ncks -O -d time,"${YYYY}-${MM}-${DD} 00:00:00","${YYYY}-${MM}-${DD} 23:00:00",6 -d latitude,0,,${STRIDE} -d longitude,0,,${STRIDE} /glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.pl/${YYYY}${MM}/e5.oper.an.pl.128_130_t.ll025sc.${YYYY}${MM}${DD}00_${YYYY}${MM}${DD}23.nc T.nc
echo "DONE"
```