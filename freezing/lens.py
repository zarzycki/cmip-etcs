import os
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import calpreciptype
from tqdm import tqdm
import numba

### This function calculates ptype based on bourgoin
# ptype = ntim x nlat   x nlon
# T     = ntim x nlev   x nlat x nlon
# Q     = ntim x nlat   x nlon
# pmid  = ntim x nlev   x nlat x nlon 
# pint  = ntim x nlev+1 x nlat x nlon
# zint2 = ntim x nlev+1 x nlat x nlon
# ntim (int)
# nlon (int)
# nlat (int)
#@numba.jit
def calc_ptype(ptype,T,Q,pmid,pint,zint2,ntim,nlon,nlat):

  for zz in range(ntim):
    for jj in tqdm(range(nlon)):
      for ii in range(nlat):
        ptype[zz,ii,jj] = calpreciptype.calwxt_bourg(np.random.rand(2),9.80665,T[zz,:,ii,jj],Q[zz,:,ii,jj],pmid[zz,:,ii,jj],pint[zz,:,ii,jj],zint2[zz,:,ii,jj])
  
  return ptype 
  



#python -m numpy.f2py -c calpreciptype.f90 -m calpreciptype
print(calpreciptype.calwxt_ramer.__doc__)
print(calpreciptype.calwxt_bourg.__doc__)
print(calpreciptype.calwxt_revised.__doc__)
print(calpreciptype.calwxt.__doc__)

grav=9.81
d608=0.608
rog=287.04/grav
h1=1.0
d00=0.0

FINIX=24

#--  data file name
fname  = "_T.nc"
ds = xr.open_dataset(fname)
T = ds.T[0:FINIX,:,:,:]

fname  = "_Q.nc"
ds = xr.open_dataset(fname)
Q = ds.Q[0:FINIX,:,:,:]

fname  = "_PS.nc"
ds = xr.open_dataset(fname)
PS=ds.PS[0:FINIX,:,:]

print(ds.hyam)

pmid = ds.hyam*ds.P0 + ds.hybm*PS
pmid.name="PMID"
pmid.attrs['units'] = 'Pa'
pmid = pmid.transpose('time', 'lev', 'lat', 'lon')

pint = ds.hyai*ds.P0 + ds.hybi*ds.PS
pint.name="PINT"
pint.attrs['units'] = 'Pa'
pint = pint.transpose('time', 'ilev', 'lat', 'lon')

rh = mpcalc.relative_humidity_from_specific_humidity(Q, T, pmid)
rharr = xr.DataArray(np.asarray(rh), dims=('time', 'lev', 'lat', 'lon'), coords=Q.coords)
rharr.name="RH"

td = mpcalc.dewpoint_from_specific_humidity(Q, T, pmid).to('K')
tdarr = xr.DataArray(np.asarray(td), dims=('time', 'lev', 'lat', 'lon'), coords=Q.coords)
tdarr.name="TD"
tdarr.attrs['units'] = 'K'


TW = (T-273.15)*np.arctan(0.151977*((rharr*100.)+8.313659)**(0.5))+np.arctan((T-273.15) + (rharr*100.)) - np.arctan((rharr*100.) - 1.676331) + 0.00391838*(rharr*100.)**(3/2)*np.arctan(0.023101*(rharr*100.)) - 4.686035
TW = TW + 273.15
TWET = xr.DataArray(np.asarray(TW), dims=('time', 'lev', 'lat', 'lon'), coords=PS.coords)
TWET.name="TWET"
TWET.attrs['units'] = 'K'

ptypea = xr.DataArray(0, dims=('time', 'lat', 'lon'), coords=PS.coords)
ptypea.name="PTYPEA"
ptypeb = xr.DataArray(0, dims=('time', 'lat', 'lon'), coords=PS.coords)
ptypeb.name="PTYPEB"
ptypec = xr.DataArray(0, dims=('time', 'lat', 'lon'), coords=PS.coords)
ptypec.name="PTYPEC"
ptyped = xr.DataArray(0, dims=('time', 'lat', 'lon'), coords=PS.coords)
ptyped.name="PTYPED"
print("LOOP")

print("Calculating virtual temperature")
tq = mpcalc.virtual_temperature(T, Q)
print(type(tq))
testarr = np.array(tq)
print(type(tq))

#zint2 = calc_zint(pint.values,PHIS.values,TKV.values,p1dims,LOOPIX,nlevp1,doFlip=True)
zint2 = mpcalc.pressure_to_height_std(pint).to('m')

calc_ptype(ptypeb.values,T.values,Q.values,pmid.values,pint.values,zint2,T.sizes['time'],T.sizes['lon'],T.sizes['lat'])

#rh.transpose('time', 'lev', 'lat', 'lon')

#ptype = calpreciptype.calwxt_ramer(t,pmid,rh,td)

sum = xr.merge([ptypeb,T,TWET,tdarr,rharr])

sum.to_netcdf('saved_on_disk.nc')

