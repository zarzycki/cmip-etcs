# Precipitation typing!

Prerequisites:

- numpy
- xarray
- metpy
- tqdm
- calpreciptype (**See below!**)


### Workflow

This code leans quite heavily on both xarray and numpy, with metpy used to do some thermodynamic calculations. I could probably generalize a bit more beyond xarray but eh.

1. Load **T**, **Q**, and **PS** data.
	- Q can be RH, etc., but must be converted to specific humidity at this point.
	- Assuming T is on a higher resolution vertical grid, Q can be interpolated to the same vertical coordinate.
2. Calculate **pmid** and **pint**.
	- If the data is on *constant pressure* surfaces, use the `constant_press_prep` function to process.
	- If the data is on *model* surfaces, pmid and pint need to be defined. This is probably via the use of a hybrid calculation (e.g., P(z,i,j) = a(k)\*P0 + b(k)\*PS(i,j)) but for other models a 4-D (ntime, nlev, nlat, nlon) pressure array may be provided.
3. Run `calc_ptype` wrapper.

Some important notes:

- All arrays should be ordered **top-to-bottom**! These can be flipped in step 1 by using the ::-1 convention.

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