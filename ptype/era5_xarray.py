import os
import sys
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import calpreciptype
from tqdm import tqdm
import argparse
import numba

def MyWho():
  for name in globals().copy().keys():
    print(name)

def np_ffill(arr, axis):
    '''https://stackoverflow.com/a/60941040'''
    idx_shape = tuple([slice(None)] + [np.newaxis] * (len(arr.shape) - axis - 1))
    idx = np.where(~np.isnan(arr), np.arange(arr.shape[axis])[idx_shape], 0)
    np.maximum.accumulate(idx, axis=axis, out=idx)
    slc = [np.arange(k)[tuple([slice(None) if dim==i else np.newaxis
        for dim in range(len(arr.shape))])]
        for i, k in enumerate(arr.shape)]
    slc[axis] = idx
    return arr[tuple(slc)]

@numba.jit
def calc_ptype(ptype,T,Q,pmid,pint,zint2,ntim,nlon,nlat):
  grav=9.80665

  for zz in range(ntim):
    for jj in tqdm(range(nlon)):
      for ii in range(nlat):
        ptype[zz,ii,jj] = calpreciptype.calwxt_bourg(np.random.rand(2),grav,T[zz,:,ii,jj],Q[zz,:,ii,jj],pmid[zz,:,ii,jj],pint[zz,:,ii,jj],zint2[zz,:,ii,jj])
  
  return ptype 

@numba.jit
def calc_zint(pint,PHIS,TKV,p1dims,ntim,nlevp1,doFlip=False):

  # doFlip flips ordering

  grav=9.80665
  rog=287.04/grav
  
  # Now calculate zint
  zint=np.empty(p1dims)

  # Note, we always want to build zint from the "bottom" up (0 -> nlev) as an integral.
  # If p is top -> bottom, we can built a numpy zint by "flipping" the p and TKV indices
  if doFlip:
    for zz in range(ntim):
      zint[zz,0,:,:] = PHIS/grav
      for kk in range(nlevp1-1):
        kkf = nlevp1-kk-1
        zint[zz,kk+1,:,:] = zint[zz,kk,:,:] + rog * TKV[zz,kkf-1,:,:] * np.log(pint[zz,kkf,:,:]/pint[zz,kkf-1,:,:])
    
    # If zint doesn't match orientation of other vars, we can just flip the numpy array before
    # building an xr DataArray
    zint = zint[:,::-1,:,:]

  else:
    for zz in range(ntim):
      zint[zz,0,:,:] = PHIS/grav
      for kk in range(nlevp1-1):
        zint[zz,kk+1,:,:] = zint[zz,kk,:,:] + rog * TKV[zz,kk,:,:] * np.log(pint[zz,kk,:,:]/pint[zz,kk+1,:,:])

  return zint




#python -m numpy.f2py --opt='-O3' -c calpreciptype.f90 -m calpreciptype
#print(calpreciptype.calwxt_ramer.__doc__)
#print(calpreciptype.calwxt_bourg.__doc__)
#print(calpreciptype.calwxt_revised.__doc__)
#print(calpreciptype.calwxt.__doc__)

#ncar_pylib 20200417

parser = argparse.ArgumentParser()
parser.add_argument('stix',type=int)
args = parser.parse_args()

## Constants
grav=9.80665
d608=0.608
rog=287.04/grav
h1=1.0
d00=0.0

# Settings
stride=1

# Indices
STAIX=args.stix
FINIX=STAIX+6
LOOPIX=FINIX-STAIX

#--  data file name
fname  = "/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.pl/201301/e5.oper.an.pl.128_130_t.ll025sc.2013010100_2013010123.nc"
ds = xr.open_mfdataset(fname, coords="minimal", data_vars="minimal", compat="override", combine='by_coords')

### Check if the maximum time length is beyond finix, if so, truncate
timeArr = ds.time

### Fix indices
MAXIX=len(timeArr.values)
FINIX=min([FINIX, MAXIX])
LOOPIX=FINIX-STAIX

### Set formatted date for output
formattedDate = timeArr.dt.strftime('%Y%m%d%H')

print("Getting T")
T = ds.T[STAIX:FINIX,:,::stride,::stride]
T = T.rename({'time':'time','level':'lev','latitude':'lat','longitude':'lon'})
T.load()
ds.close()

# Do some coordinate stuff with T
TLLLcoords = T.coords
TLLcoords = {'time': T.coords['time'], 'lat': T.coords['lat'], 'lon': T.coords['lon']}

print("Getting Q")
fname="/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.pl/201301/e5.oper.an.pl.128_133_q.ll025sc.2013010100_2013010123.nc"
ds = xr.open_mfdataset(fname, coords="minimal", data_vars="minimal", compat="override", combine='by_coords')
Q = ds.Q[STAIX:FINIX,:,::stride,::stride]
Q = Q.rename({'time':'time','level':'lev','latitude':'lat','longitude':'lon'})
Q.load()
ds.close()

print("Getting PS")
## This handles model level data by reading PS and "building" p at model levels
fname  = "/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.sfc/201301/e5.oper.an.sfc.128_134_sp.ll025sc.2013010100_2013013123.nc"
ds = xr.open_mfdataset(fname, combine='by_coords')
PS=ds.SP[STAIX:FINIX,::stride,::stride]
PS = PS.rename({'time':'time','latitude':'lat','longitude':'lon'})
PS.load()
ds.close()

print("Expanding levs...")
### This code reads a level file and expands the dimensions to match lat/lon
lev = T.coords['lev']  #.values
pmid = lev.expand_dims({'time': T.coords['time'],'lat': T.coords['lat'], 'lon': T.coords['lon']},(0,2,3))
pmid.name="P4D"
pmid = pmid * 100.
pmid.attrs['units'] = 'Pa'

print("Getting PHIS")
## This handles model level data by reading PS and "building" p at model levels
fname  = "/glade/u/home/zarzycki/rda/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"
ds = xr.open_dataset(fname)
PHIS=ds.Z[0,::stride,::stride]
PHIS = PHIS.rename({'latitude':'lat','longitude':'lon'})
ds.close()


### Handle vertical coordinate stuff

# Set up a pint numpy array
pdims = np.asarray(T.shape)
nlev = pdims[1]
nlevp1 = nlev+1
p1dims = pdims
p1dims[1] = nlevp1
pintnp = np.empty(p1dims)

## Expand PS to a 4-D time x lev x lat x lon array
psexpand = PS.expand_dims({'lev': T.coords['lev']},(1))
psexpand.name="PS4D"
## Expand PS to a 4-D time x ilev x lat x lon array
psexpandi = PS.expand_dims({'ilev': pintnp[0,:,0,0]},(1))
psexpandi.name="PSI4D"

# First, set all "below ground" points to -1 in the 4D pmid array
pmid  = pmid.where(pmid < psexpand, -1.)

# Use this information to fix other variables as well.
# Use np_ffill which was shamelessly stolen from StackOverflow
T = T.where(pmid > 0)
T.values = np_ffill(T.values,1)
Q = Q.where(pmid > 0)
Q.values = np_ffill(Q.values,1)

#for zz in range(LOOPIX):
#  for jj in tqdm(range(T.sizes['lon'])):
#    for ii in range(T.sizes['lat']):
#      T[zz,:,ii,jj] = fill_below_ground(T[zz,:,ii,jj])

# Find the maximum valid (positive) value of P at each pt
pmax  = pmid.max(axis=1).expand_dims({'lev': T.coords['lev']},(1))
pmaxi = pmid.max(axis=1).expand_dims({'ilev': pintnp[0,:,0,0]},(1))

## Wherever pmid < 0, we are "below" ground, so set to lowest pmid layer we just found
pmid = pmid.where(pmid > 0, pmax)

# Set lowest model level P interface to PS.
model_lid=10.  # lid in Pa to avoid 0 issues
pintnp[:,nlevp1-1,:,:] = PS
pintnp[:,0,:,:] = pmid[:,0,:,:].values - (pmid[:,1,:,:].values - pmid[:,0,:,:].values)
pintnp[:,0,:,:] = np.where(pintnp[:,0,:,:] > 0, pintnp[:,0,:,:], model_lid)
## Calculate interface levels by "splitting"
pintnp[:,1:nlevp1-1,:,:] = np.add( pmid[:,0:nlev-1,:,:].values , pmid[:,1:nlev,:,:].values ) / 2.
## Calculate interface levels using S+B (1981) eq. 3.18  ### NOT WORKING!
#delpmid = np.subtract( pmid[:,1:nlev,:,:].values, pmid[:,0:nlev-1,:,:].values )
#pintnp[:,1:nlevp1-1,:,:] = np.where(delpmid > 0.00001, np.exp( np.divide ( np.subtract( np.multiply(pmid[:,1:nlev,:,:].values,np.log(pmid[:,1:nlev,:,:].values)), np.multiply(pmid[:,0:nlev-1,:,:].values,np.log(pmid[:,0:nlev-1,:,:].values))  ) , delpmid ) - 1.0 ), pmid[:,0:nlev-1,:,:].values)
## NOTE: The error between S+B and straight average in the troposphere is <0.06% (0.0006) for L60.

# Turn pint into an xarray
pint = xr.DataArray(np.asarray(pintnp), dims=('time', 'ilev', 'lat', 'lon'), coords={'time': T.coords['time'], 'ilev' : pintnp[0,:,0,0] , 'lat': T.coords['lat'], 'lon': T.coords['lon']})

# Wherever pint2 >= lowest model level p (pmax) we MUST set to PS for interface
# the lower interface to pbot must be PS by definition
pint = pint.where(pint < pmaxi, psexpandi)
pint.name="PINT2"
pint.attrs['units'] = 'Pa'
        
# Do some cleanup here
del psexpand, psexpandi, pintnp

# Now do some derived variables...
print("Calculating virtual temperature")
tq = mpcalc.virtual_temperature(T, Q)
TKV = xr.DataArray(np.asarray(tq), dims=('time', 'lev', 'lat', 'lon'), coords=TLLLcoords)
TKV.name="TKV"
del tq



# for name, size in sorted(((name, sys.getsizeof(value)) for name, value in globals().items()),
#                          key= lambda x: -x[1])[:10]:
#     print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
MyWho()


#zint = np.empty(T.sizes)
#TVBAR = np.empty(T.sizes)

#zint = xr.DataArray(np.empty(   ), dims=('time', 'lev', 'lat', 'lon'), coords=TLLLcoords)
#TVBAR = xr.DataArray(np.empty(   ), dims=('time', 'lev', 'lat', 'lon'), coords=TLLLcoords)

# 
# c calculate geopotential height if initial z available [hydrostatic eq]
# 
#       ZH(1) = ZSFC
#       DO NL = 2,NLVL
# C same as [ln(p(nl)+ln(p(nl-1)]
#           TVBAR = (TKV(NL)*DLOG(P(NL))+TKV(NL-1)*DLOG(P(NL-1)))/
#      +            DLOG(P(NL)*P(NL-1))
#           ZH(NL) = ZH(NL-1) + RDAG*TVBAR*DLOG(P(NL-1)/P(NL))
#       END DO


#>> # Now calculate zint
#>> zint=np.empty(p1dims)
#>> 
#>> # Note, we always want to build zint from the "bottom" up (0 -> nlev) as an integral.
#>> # If p is top -> bottom, we can built a numpy zint by "flipping" the p and TKV indices
#>> for zz in tqdm(range(LOOPIX)):
#>>   zint[zz,0,:,:] = PHIS/grav
#>>   for kk in range(nlevp1-1):
#>>     
#>>     #print(kk,' ',pint[zz,kk+1,0,0].values)
#>>     #TVBAR[zz,kk+1,:,:] = ( TKV[zz,kk+1,:,:]*np.log(pint[zz,kk+1,:,:]) + TKV[zz,kk,:,:]*np.log(pint[zz,kk,:,:]) ) / np.log(pint[zz,kk,:,:] * pint[zz,kk+1,:,:])
#>>     #zint[zz,kk+1,:,:] = zint[zz,kk+1,:,:] + rog * TVBAR[zz,kk+1,:,:] * np.log(pint[zz,kk,:,:]/pint[zz,kk+1,:,:])
#>>     kkf = nlevp1-kk-1
#>>     zint[zz,kk+1,:,:] = zint[zz,kk,:,:] + rog * TKV[zz,kkf-1,:,:] * np.log(pint[zz,kkf,:,:]/pint[zz,kkf-1,:,:])
#>> 
#>> # If zint doesn't match orientation of other vars, we can just flip the numpy array before
#>> # building an xr DataArray
#>> zint = zint[:,::-1,:,:]

### ZINT function
zint = calc_zint(pint,PHIS,TKV,p1dims,LOOPIX,nlevp1,doFlip=True)
#>>>
# Built xr DataArray
zint2 = xr.DataArray(np.asarray(zint), dims=('time', 'ilev', 'lat', 'lon'), coords={'time': T.coords['time'], 'ilev' : pint[0,:,0,0] , 'lat': T.coords['lat'], 'lon': T.coords['lon']})
zint2.name="ZINT2"
del zint
  



#zint_std=mpcalc.pressure_to_height_std(pint2).to('m')
#zint2 = xr.DataArray(np.asarray(zint_std), dims=('time', 'ilev', 'lat', 'lon'), coords={'time': T.coords['time'], 'ilev' : pint[0,:,0,0], 'lat': T.coords['lat'], 'lon': T.coords['lon']})
#zint2.name="ZINT2"
#zint2.attrs['units'] = 'm'







#### This is all non-model specific stuff

##rh = mpcalc.relative_humidity_from_specific_humidity(Q, T, pmid)
##rharr = xr.DataArray(np.asarray(rh), dims=('time', 'lev', 'lat', 'lon'), coords=TLLLcoords)
##rharr.name="RH"

##td = mpcalc.dewpoint_from_specific_humidity(Q, T, pmid).to('K')
##tdarr = xr.DataArray(np.asarray(td), dims=('time', 'lev', 'lat', 'lon'), coords=TLLLcoords)
##tdarr.name="TD"
##tdarr.attrs['units'] = 'K'

##TW = (T-273.15)*np.arctan(0.151977*((rharr*100.)+8.313659)**(0.5))+np.arctan((T-273.15) + (rharr*100.)) - np.arctan((rharr*100.) - 1.676331) + 0.00391838*(rharr*100.)**(3/2)*np.arctan(0.023101*(rharr*100.)) - 4.686035
##TW = TW + 273.15
##TWET = xr.DataArray(np.asarray(TW), dims=('time', 'lev', 'lat', 'lon'), coords=TLLLcoords)
##TWET.name="TWET"
##TWET.attrs['units'] = 'K'

#ptypea = xr.DataArray(0, dims=('time', 'lat', 'lon'), coords=TLLcoords)
#ptypea.name="PTYPEA"
ptypeb = xr.DataArray(0, dims=('time', 'lat', 'lon'), coords=TLLcoords)
ptypeb.name="PTYPEB"
#ptypec = xr.DataArray(0, dims=('time', 'lat', 'lon'), coords=TLLcoords)
#ptypec.name="PTYPEC"

calc_ptype(ptypeb.values,T.values,Q.values,pmid.values,pint.values,zint2.values,LOOPIX,T.sizes['lon'],T.sizes['lat'])

#for zz in range(LOOPIX):
#  for jj in tqdm(range(T.sizes['lon'])):
#    for ii in range(T.sizes['lat']):
 
     #tt = T[zz,:,ii,jj]
      #qq = Q[zz,:,ii,jj]
      #pp = pmid[zz,:,ii,jj]
      #pi = pint[zz,:,ii,jj]
      #zi = mpcalc.pressure_to_height_std(pint[zz,:,ii,jj]).to('m')
      #rr = rharr[zz,:,ii,jj]
      #dd = tdarr[zz,:,ii,jj]
      #tw = TWET[zz,:,ii,jj]

      #if jj == 0:
        #print(grav)
        #print(tt)
        #print(qq)
        #print(pp)
        #print(pi)
        #print(zi)
      
      # everything needs to be arranged top to bottom
      #ptypea[zz,ii,jj] = calpreciptype.calwxt_ramer(tt,pp,rr,dd)
      #ptypeb[zz,ii,jj] = calpreciptype.calwxt_bourg(np.random.rand(2),grav,tt,qq,pp,pi,zi)
      #ptypec[zz,ii,jj] = calpreciptype.calwxt(tt,qq,pp,pi,d608,rog,1.0E-10,zi,tw)
      #ptyped[zz,ii,jj] = calpreciptype.calwxt_revised(tt,qq,pp,pi,d608,rog,1.0E-10,zi,tw)

      #ptypea[zz,ii,jj] = calpreciptype.calwxt_ramer(T[zz,::-1,ii,jj],pmid[zz,::-1,ii,jj],rharr[zz,::-1,ii,jj],tdarr[zz,::-1,ii,jj])
      

####ptypeb[zz,ii,jj] = calpreciptype.calwxt_bourg(np.random.rand(2),grav,T[zz,:,ii,jj],Q[zz,:,ii,jj],pmid[zz,:,ii,jj],pint[zz,:,ii,jj],zint2[zz,:,ii,jj])
      #ptypec[zz,ii,jj] = calpreciptype.calwxt_bourg(np.random.rand(2),grav,T[zz,::-1,ii,jj],Q[zz,::-1,ii,jj],pmid[zz,::-1,ii,jj],pint[zz,::-1,ii,jj],zint2[zz,::-1,ii,jj])

      # add logic for pmid dims

      #ptypea[zz,ii,jj] = calpreciptype.calwxt_ramer(T[zz,:,ii,jj],pmid[zz,:,ii,jj],rharr[zz,:,ii,jj],tdarr[zz,:,ii,jj])
      #ptypeb[zz,ii,jj] = calpreciptype.calwxt_bourg(np.random.rand(2),grav,T[zz,:,ii,jj],Q[zz,:,ii,jj],pmid[zz,:,ii,jj],pint[zz,:,ii,jj],mpcalc.pressure_to_height_std(pint[zz,:,ii,jj]).to('m'))
      #ptypec[zz,ii,jj] = calpreciptype.calwxt(T[zz,:,ii,jj],Q[zz,:,ii,jj],pmid[zz,:,ii,jj],pint[zz,:,ii,jj],d608,rog,1.0E-10,mpcalc.pressure_to_height_std(pint[zz,:,ii,jj]).to('m'),TWET[zz,:,ii,jj])
      #ptyped[zz,ii,jj] = calpreciptype.calwxt_revised(T[zz,:,ii,jj],Q[zz,:,ii,jj],pmid[zz,:,ii,jj],pint[zz,:,ii,jj],d608,rog,1.0E-10,mpcalc.pressure_to_height_std(pint[zz,:,ii,jj]).to('m'),TWET[zz,:,ii,jj])

      #print(ptype[zz,ii,jj])

#rh.transpose('time', 'lev', 'lat', 'lon')

#ptype = calpreciptype.calwxt_ramer(t,pmid,rh,td)

sum = xr.merge([ptypeb,T,pmid,zint2,pint,TKV,Q])

outFileName = os.path.join('/glade/u/home/zarzycki/scratch/ERA5-ptype/' + formattedDate[0].values + '_ptype.nc')

sum.to_netcdf(outFileName, unlimited_dims = "time", encoding={'time':{'units':'days since 1950-01-01'}} )

