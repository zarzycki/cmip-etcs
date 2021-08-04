import os
import sys
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
import calpreciptype
from tqdm import tqdm

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

# Inputs:
# T (ntim,nlev,nlat,nlon) in K
# Q (ntim,nlev,nlat,nlon) in kg/kg
# pmid (ntim,nlev,nlat,nlon) in Pa
# pint (ntim,nlev+1,nlat,nlon) in Pa
# zint (ntim,nlev+1,nlat,nlon) in m
# ntim (int)
# nlev (int)
# nlon (int)
# nlat (int)
def calc_ptype(T,Q,pmid,pint,zint,ntim,nlev,nlat,nlon):

    # Define constants
    grav=9.80665
    d608=0.608
    rog=287.04/grav
    h1=1.0
    d00=0.0
    epsq=1.e-10
    rheps = 1.e-4

    # Define any coord stuff...
    nlevp1 = nlev+1

    # If zint passed in -> -1, then we will calculate zint internally using T, Q, and pint
    if isinstance(zint, int):
        if (zint < 0):
            print("Calculating interface heights")
            zint = calc_zint(pint.values,mpcalc.virtual_temperature(T, Q).values,ntim,nlevp1)

    print("Calculating relative humidity")
    rh = mpcalc.relative_humidity_from_specific_humidity(pmid, T, Q)
    # Where RH > 100%, set to 100%
    rh = rh.where(rh < 1.0, 1.0)
    # Where RH < 0 (high in atmo), just set to tiny positive number
    rh = rh.where(rh > rheps, rheps)
    rh = xr.DataArray(np.asarray(rh), dims=('time', 'lev', 'lat', 'lon'), coords=T.coords)
    rh.name="RH"

    print("Calculating dew point temperature")
    #CMZ: calculating Td from q leads to nans in upper troposphere. Using corrected RH instead seems good!
    #td = mpcalc.dewpoint_from_specific_humidity(pmid, T, Q).metpy.convert_units('K')
    td = mpcalc.dewpoint_from_relative_humidity(T, rh).metpy.convert_units('K')
    td = xr.DataArray(np.asarray(td), dims=('time', 'lev', 'lat', 'lon'), coords=T.coords)
    td.name="TD"
    
    print("Calculating wet bulb temperature")
    TW = (T-273.15)*np.arctan(0.151977*((rh*100.)+8.313659)**(0.5))+np.arctan((T-273.15) + (rh*100.)) - np.arctan((rh*100.) - 1.676331) + 0.00391838*(rh*100.)**(3/2)*np.arctan(0.023101*(rh*100.)) - 4.686035
    TW = TW + 273.15
    TW = xr.DataArray(np.asarray(TW), dims=('time', 'lev', 'lat', 'lon'), coords=T.coords)
    TW.name="TWET"
    TW.attrs['units'] = 'K'
    
    # We can use this to dump diags re: the thermodynamic quantities from metpy -- set to true to dump.
    print_diags=False
    if print_diags:
        print("NOTE: Printing diags inside of calc_ptype!")
        diags = xr.merge([T,Q,pmid,pint,rh,td,TW])
        outFileName = os.path.join('./diag.nc')
        diags.to_netcdf(outFileName, unlimited_dims = "time")
    
    # These are all xarray, but lets check and convert to numpy to speed things up
    # xr arrays can get passed into the precip type calcs, but about 10x slower
    print("Converting xarray -> numpy")
    if hasattr(T, 'values'):
        T = T.values
    if hasattr(Q, 'values'):
        Q = Q.values
    if hasattr(pmid, 'values'):
        pmid = pmid.values
    if hasattr(pint, 'values'):
        pint = pint.values
    if hasattr(rh, 'values'):
        rh = rh.values
    if hasattr(td, 'values'):
        td = td.values
    if hasattr(TW, 'values'):
        TW = TW.values

    print("Calculating PTYPE")
    # Initialize array to hold ptypes for each algo
    a = np.array([0,0,0])
    # Initialize ptype array TLL
    ptype = np.empty((ntim, nlat, nlon))
    
    for zz in range(ntim):
        for jj in tqdm(range(nlon)):
            for ii in range(nlat):
                # Calculate ptype at this time + location with 3 different algos
                a[0] = calpreciptype.calwxt_bourg(np.random.rand(2),grav,T[zz,:,ii,jj],Q[zz,:,ii,jj],pmid[zz,:,ii,jj],pint[zz,:,ii,jj],zint[zz,:,ii,jj])
                a[1] = calpreciptype.calwxt_ramer(T[zz,:,ii,jj],pmid[zz,:,ii,jj],rh[zz,:,ii,jj],td[zz,:,ii,jj])
                a[2] = calpreciptype.calwxt_revised(T[zz,:,ii,jj],Q[zz,:,ii,jj],pmid[zz,:,ii,jj],pint[zz,:,ii,jj],d608,rog,epsq,zint[zz,:,ii,jj],TW[zz,:,ii,jj])

                # Count number of ptypes from each algo for this cell
                counts = np.bincount(a)

                if (counts < 2).all():
                    # This means we have 1, 1, 1 for ptype, so pick Bourgouin
                    ptype[zz,ii,jj] = a[0]
                else:
                    # This means we have at least 2/1 (or, better, 3/0) agreement, pick "winner"
                    ptype[zz,ii,jj] = np.argmax(counts)
                    
    return ptype

def calc_zint(pint,TKV,ntim,nlevp1,doFlip=False):

    grav=9.80665
    rog=287.04/grav
  
    p1dims = np.asarray(pint.shape)
    
    # Now calculate zint
    zint=np.empty(p1dims)

    # Note, we always want to build zint from the "bottom" up (0 -> nlev) as an integral.
    # If p is top -> bottom, we can built a numpy zint by "flipping" the p and TKV indices
    if doFlip:
        for zz in range(ntim):
            zint[zz,0,:,:] = 0.0
            for kk in tqdm(range(nlevp1-1)):
                kkf = nlevp1-kk-1
                zint[zz,kk+1,:,:] = zint[zz,kk,:,:] + rog * TKV[zz,kkf-1,:,:] * np.log(pint[zz,kkf,:,:]/pint[zz,kkf-1,:,:])
    
        # If zint doesn't match orientation of other vars, we can just flip the numpy array before
        # building an xr DataArray
        zint = zint[:,::-1,:,:]

    else:
        for zz in range(ntim):
            zint[zz,0,:,:] = 0.0
            for kk in tqdm(range(nlevp1-1)):
                zint[zz,kk+1,:,:] = zint[zz,kk,:,:] + rog * TKV[zz,kk,:,:] * np.log(pint[zz,kk,:,:]/pint[zz,kk+1,:,:])

    return zint

def constant_press_prep(T,Q,PS,lev,ntime,nlev,nlat,nlon):
    
    # Constants go here:
    model_lid=10.  # Pa
    
    print("Expanding levs...")
    ### This code takes a single-dimensioned lev xr.DataArray and expands it to 4dim ntim,nlev,nlat,nlon
    pmid = lev.expand_dims({'time': T.coords['time'],'lat': T.coords['lat'], 'lon': T.coords['lon']},(0,2,3))
    pmid.name="P4D"
    pmid = pmid * 100.
    pmid.attrs['units'] = 'Pa'
    
    # Set up a pint np array
    pdims = [ntime,nlev,nlat,nlon]
    nlevp1 = nlev+1
    p1dims = pdims
    p1dims[1] = nlevp1
    pintnp = np.empty(p1dims)

    ## Expand PS xr.DataArray to a 4-D time x lev x lat x lon array
    psexpand = PS.expand_dims({'lev': T.coords['lev']},(1))
    psexpand.name="PS4D"
    ## Expand PS xr.DataArray to a 4-D time x ilev x lat x lon array
    psexpandi = PS.expand_dims({'ilev': pintnp[0,:,0,0]},(1))
    psexpandi.name="PSI4D"

    print("Filling below-ground points in T, Q, and midpoint pressures...")
    # Set all "below ground" points to -1 in the 4D pmid array
    pmid  = pmid.where(pmid < psexpand, -1.)

    # Use this information to fix other variables as well.
    # This sets all "below ground" points to the same as the lowest valid value (lowest above-ground "layer").
    T = T.where(pmid > 0)
    T.values = np_ffill(T.values,1)
    Q = Q.where(pmid > 0)
    Q.values = np_ffill(Q.values,1)

    # Find the maximum valid (positive) value of P at each pt
    pmax  = pmid.max(axis=1).expand_dims({'lev': T.coords['lev']},(1))
    pmaxi = pmid.max(axis=1).expand_dims({'ilev': pintnp[0,:,0,0]},(1))

    # Wherever pmid < 0, we are "below" ground, so set to lowest pmid layer we just found
    # IOW, we are setting all below ground points to the lowest above-ground valid value
    pmid = pmid.where(pmid > 0, pmax)

    print("Building interface pressures...")
    # First fill in "edge" cases -> top and bottom...
    # Set lowest model level P interface to PS.
    pintnp[:,nlevp1-1,:,:] = PS
    # Set highest model level to an extrapolation
    pintnp[:,0,:,:] = pmid[:,0,:,:].values - (pmid[:,1,:,:].values - pmid[:,0,:,:].values) / 2.
    # Define a model lid to avoid p -> 0 errors.
    pintnp[:,0,:,:] = np.where(pintnp[:,0,:,:] > 0, pintnp[:,0,:,:], model_lid)
    
    # Calculate remainder of "inner" interface levels by splitting pmid layers in half
    pintnp[:,1:nlevp1-1,:,:] = np.add( pmid[:,0:nlev-1,:,:].values , pmid[:,1:nlev,:,:].values ) / 2.

    # Turn pint into an xr.DataArray to match pmid, T, and Q
    pint = xr.DataArray(np.asarray(pintnp), dims=('time', 'ilev', 'lat', 'lon'), coords={'time': T.coords['time'], 'ilev' : pintnp[0,:,0,0] , 'lat': T.coords['lat'], 'lon': T.coords['lon']})
    pint.name="PINT"
    pint.attrs['units'] = 'Pa'
    
    # Wherever pint >= lowest model level p (pmax) we are either A: below ground or B: at the sfc interface
    # Therefore, set all of these to PS since the lower interface to pbot must be PS by definition
    pint = pint.where(pint < pmaxi, psexpandi)

    # Do some cleanup here for funsies
    del psexpand, psexpandi, pintnp
    
    # Return xarray DataArrays for T, Q, pmid, and pint to main function
    return T, Q, pmid, pint

### MAIN PROGRAM BEGINS HERE!

dataset = "LENS"
print("dataset "+dataset)

if dataset == "LENS":
    
    #--  data file name
    fname  = "./sample-LENS/T.nc"
    ds = xr.open_dataset(fname)
    T = ds.T[:,:,:,:]
    ds.close()
    
    fname  = "./sample-LENS/Q.nc"
    ds = xr.open_dataset(fname)
    Q = ds.Q[:,:,:,:]
    ds.close()
    
    fname  = "./sample-LENS/PS.nc"
    ds = xr.open_dataset(fname)
    PS=ds.PS[:,:,:]
    
    pmid = ds.hyam*ds.P0 + ds.hybm*PS
    pmid.name="PMID"
    pmid.attrs['units'] = 'Pa'
    pmid = pmid.transpose('time', 'lev', 'lat', 'lon')

    pint = ds.hyai*ds.P0 + ds.hybi*ds.PS
    pint.name="PINT"
    pint.attrs['units'] = 'Pa'
    pint = pint.transpose('time', 'ilev', 'lat', 'lon')

    ds.close()
    
    # Get time-lev-lat-lon coords
    # ... and time-lat-lon coords
    TLLLcoords = T.coords
    TLLcoords = {'time': T.coords['time'], 'lat': T.coords['lat'], 'lon': T.coords['lon']}
    
    print ("Getting dim sizes")
    ntime = T.sizes['time']
    nlev  = T.sizes['lev']
    nlat  = T.sizes['lat']
    nlon  = T.sizes['lon']
    
else:
    
    # Settings
    stride=1
    numtimes=8
    # Indices
    STAIX=0
    FINIX=STAIX+numtimes

    #--  data file name
    fname = "./sample-ERA5/T.nc"
    ds    = xr.open_mfdataset(fname, coords="minimal", data_vars="minimal", compat="override", combine='by_coords')

    ### Check if the maximum time length is beyond finix, if so, truncate
    timeArr = ds.time
    MAXIX=len(timeArr.values)
    
    ### Fix indices
    FINIX=min([FINIX, MAXIX])
    
    # Calculate loop indices
    LOOPIX=FINIX-STAIX

    ### Set formatted date for output
    formattedDate = timeArr.dt.strftime('%Y%m%d%H')

    ### Get Temperature
    print("Getting T")
    T = ds.T[STAIX:FINIX,:,::stride,::stride]
    T = T.rename({'time':'time','level':'lev','latitude':'lat','longitude':'lon'})
    T.load()
    ds.close()

    # Get time-lev-lat-lon coords
    # ... and time-lat-lon coords
    TLLLcoords = T.coords
    TLLcoords = {'time': T.coords['time'], 'lat': T.coords['lat'], 'lon': T.coords['lon']}

    ### Get specific humidity
    print("Getting Q")
    fname="./sample-ERA5/Q.nc"
    ds = xr.open_mfdataset(fname, coords="minimal", data_vars="minimal", compat="override", combine='by_coords')
    Q = ds.Q[STAIX:FINIX,:,::stride,::stride]
    Q = Q.rename({'time':'time','level':'lev','latitude':'lat','longitude':'lon'})
    Q.load()
    ds.close()

    ## This handles model level data by reading PS and "building" p at model levels
    print("Getting PS")
    fname  = "./sample-ERA5/PS.nc"
    ds = xr.open_mfdataset(fname, combine='by_coords')
    PS=ds.SP[STAIX:FINIX,::stride,::stride]
    PS = PS.rename({'time':'time','latitude':'lat','longitude':'lon'})
    PS.load()
    ds.close()

    print ("Getting dim sizes")
    ntime = LOOPIX
    nlev  = T.sizes['lev']
    nlat  = T.sizes['lat']
    nlon  = T.sizes['lon']

    ## From constant pressure level T, Q, + PS, build pmid and pint arrays
    T, Q, pmid, pint = constant_press_prep(T,Q,PS,T.coords['lev'],ntime,nlev,nlat,nlon)

## Calculate ptype
ptype = calc_ptype(T,Q,pmid,pint,-1,ntime,nlev,nlat,nlon)

## Write stuff to disk
print("Writing to disk")
ptype = xr.DataArray(ptype, dims=('time', 'lat', 'lon'), coords=TLLcoords)
ptype.name="PTYPE"

output = xr.merge([ptype])

if dataset == "LENS":
    output.to_netcdf('LENS-ptype.nc')
else:
    outFileName = os.path.join('./ERA-ptype.nc')
    output.to_netcdf(outFileName, unlimited_dims = "time", encoding={'time':{'units':'days since 1950-01-01'}} )