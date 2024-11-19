import xarray as xr
import numpy as np
import glob
import matplotlib.pyplot as plt
import os

# Create img directory if it doesn't exist
os.makedirs('img', exist_ok=True)

def make_plot(precip, year, prefix, lons, lats):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    # Handle different time dimension names
    time_dim = 'time' if 'time' in precip.dims else 'initial_time0_hours'
    mean_precip = precip.mean(time_dim) * 86400
    if 'forecast_time1' in mean_precip.dims:
        mean_precip = mean_precip.mean('forecast_time1')

    # Set up the projection and plot
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Set map extent
    ax.set_extent([-90.285645, -62.325439, 36.958937, 47.517426], crs=ccrs.PlateCarree())

    # Add map features
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
    ax.gridlines(draw_labels=True)

    # Plot precipitation data
    cf = ax.contourf(lons, lats, mean_precip,
                     levels=np.linspace(0, 10, 21),
                     transform=ccrs.PlateCarree())

    plt.colorbar(cf, label='Precipitation Rate (mm/day)')
    plt.title(f'{prefix} Mean Precipitation Rate {year}')
    plt.savefig(f'img/{prefix}_{year}.png', dpi=300, bbox_inches='tight')
    plt.close()

def analyze_CFSR(min_year=None, max_year=None, file_pattern='./CFSR/regrid_CFSR.PREC.*.nc'):
    print(file_pattern)
    for file in sorted(glob.glob(file_pattern)):
        year = int(file.split('.')[-2])
        if (min_year and year < min_year) or (max_year and year > max_year):
            continue
        ds = xr.open_dataset(file)
        mask = (ds.lon >= LON_MIN) & (ds.lon <= LON_MAX) & (ds.lat >= LAT_MIN) & (ds.lat <= LAT_MAX)
        precip = ds.PRATE.where(mask)

        # Count only valid (non-masked) points
        valid_points = ~np.isnan(precip)
        exceeds_threshold = np.logical_and(precip > THRESHOLD_KG_M2_S, valid_points)
        precip_fraction = np.sum(exceeds_threshold) / np.sum(valid_points)

        # Open corresponding PTYPE file
        ptype_file = f'../ptypes/CFSR/CFSR_PTYPE_{year}.nc'
        ds_ptype = xr.open_dataset(ptype_file)
        mask2 = (ds_ptype.lon >= LON_MIN) & (ds_ptype.lon <= LON_MAX) & (ds_ptype.lat >= LAT_MIN) & (ds_ptype.lat <= LAT_MAX)
        ptype = ds_ptype.PTYPE.where(mask2)

        # Count PTYPE=4 points
        valid_ptype = ~np.isnan(ptype)
        is_fzra = np.logical_and(ptype == 4, valid_ptype)
        ptype_fraction = np.sum(is_fzra) / np.sum(valid_ptype)

        print(f"Year {year}: {precip_fraction:.3%} of pixels exceed {THRESHOLD_MM_DAY} mm/day, {ptype_fraction:.3%} of pixels are FZRA")
        make_plot(precip, year, 'CFSR', ds.lon, ds.lat)

        ds.close()
        ds_ptype.close()


def analyze_20CR(min_year=None, max_year=None, file_pattern='./CR20/anl_mean_*_PR.nc'):
   print(file_pattern)
   for file in sorted(glob.glob(file_pattern)):
       year = int(file.split('_')[2])
       if (min_year and year < min_year) or (max_year and year > max_year):
           continue

       # Precipitation analysis
       ds = xr.open_dataset(file)
       mask = (ds.longitude >= LON_MIN) & (ds.longitude <= LON_MAX) & (ds.latitude >= LAT_MIN) & (ds.latitude <= LAT_MAX)
       precip = ds.prate.where(mask)
       valid_points = ~np.isnan(precip)
       exceeds_threshold = np.logical_and(precip > THRESHOLD_KG_M2_S, valid_points)
       precip_fraction = np.sum(exceeds_threshold) / np.sum(valid_points)

       # PTYPE analysis
       ptype_file = f'../ptypes/CR20/CR20_PTYPE_{year}.nc'
       ds_ptype = xr.open_dataset(ptype_file)
       mask2 = (ds_ptype.lon >= LON_MIN) & (ds_ptype.lon <= LON_MAX) & (ds_ptype.lat >= LAT_MIN) & (ds_ptype.lat <= LAT_MAX)
       ptype = ds_ptype.PTYPE.where(mask2)
       valid_ptype = ~np.isnan(ptype)
       is_fzra = np.logical_and(ptype == 4, valid_ptype)
       ptype_fraction = np.sum(is_fzra) / np.sum(valid_ptype)

       print(f"Year {year}: {precip_fraction:.3%} of pixels exceed {THRESHOLD_MM_DAY} mm/day, {ptype_fraction:.3%} of pixels are FZRA")
       make_plot(precip, year, '20CR', ds.longitude, ds.latitude)

       ds.close()
       ds_ptype.close()


def analyze_ERA5(min_year=None, max_year=None, file_pattern='./ERA5/NorthAm_ERA5_PREC_*_small.nc'):
   print(file_pattern)
   for file in sorted(glob.glob(file_pattern)):
       year = int(file.split('_')[3])
       if (min_year and year < min_year) or (max_year and year > max_year):
           continue

       # Precipitation analysis
       ds = xr.open_dataset(file)
       mask = (ds.longitude >= LON_MIN) & (ds.longitude <= LON_MAX) & (ds.latitude >= LAT_MIN) & (ds.latitude <= LAT_MAX)
       precip = ds.MTPR.isel(time=slice(None, None, 6)).where(mask)
       valid_points = ~np.isnan(precip)
       exceeds_threshold = np.logical_and(precip > THRESHOLD_KG_M2_S, valid_points)
       precip_fraction = np.sum(exceeds_threshold) / np.sum(valid_points)

       # PTYPE analysis - need to process multiple files per year
       ptype_files = sorted(glob.glob(f'../ptypes/ERA5/monthly/{year}/ERAML_PTYPE_*.nc'))
       all_valid_ptype = 0
       all_fzra_points = 0

       for pfile in ptype_files:
           ds_ptype = xr.open_dataset(pfile)
           mask2 = (ds_ptype.lon >= LON_MIN) & (ds_ptype.lon <= LON_MAX) & (ds_ptype.lat >= LAT_MIN) & (ds_ptype.lat <= LAT_MAX)
           ptype = ds_ptype.PTYPE.where(mask2)
           valid_ptype = ~np.isnan(ptype)
           is_fzra = np.logical_and(ptype == 4, valid_ptype)

           all_valid_ptype += np.sum(valid_ptype)
           all_fzra_points += np.sum(is_fzra)
           ds_ptype.close()

       ptype_fraction = all_fzra_points / all_valid_ptype

       print(f"Year {year}: {precip_fraction:.3%} of pixels exceed {THRESHOLD_MM_DAY} mm/day, {ptype_fraction:.3%} of pixels are FZRA")
       make_plot(precip, year, 'ERA5', ds.longitude, ds.latitude)
       ds.close()


def analyze_JRA(min_year=None, max_year=None, file_pattern='./JRA/fcst_phy2m125_tprat_*.nc'):
   print(file_pattern)
   for file in sorted(glob.glob(file_pattern)):
       year = int(file.split('_')[-1].split('.')[0])
       if (min_year and year < min_year) or (max_year and year > max_year):
           continue

       # Precipitation analysis
       ds = xr.open_dataset(file)
       mask = (ds.g0_lon_3 >= LON_MIN) & (ds.g0_lon_3 <= LON_MAX) & (ds.g0_lat_2 >= LAT_MIN) & (ds.g0_lat_2 <= LAT_MAX)
       precip = ds.TPRAT_GDS0_SFC_ave3h * (1/86400)  # Convert mm/day to kg/m²/s
       precip = precip.where(mask)
       valid_points = ~np.isnan(precip)
       exceeds_threshold = np.logical_and(precip > THRESHOLD_KG_M2_S, valid_points)
       precip_fraction = np.sum(exceeds_threshold) / np.sum(valid_points)

       # PTYPE analysis
       ptype_file = f'../ptypes/JRA/JRA_PTYPE_{year}.nc'
       ds_ptype = xr.open_dataset(ptype_file)
       mask2 = (ds_ptype.lon >= LON_MIN) & (ds_ptype.lon <= LON_MAX) & (ds_ptype.lat >= LAT_MIN) & (ds_ptype.lat <= LAT_MAX)
       ptype = ds_ptype.PTYPE.where(mask2)
       valid_ptype = ~np.isnan(ptype)
       is_fzra = np.logical_and(ptype == 4, valid_ptype)
       ptype_fraction = np.sum(is_fzra) / np.sum(valid_ptype)

       print(f"Year {year}: {precip_fraction:.3%} of pixels exceed {THRESHOLD_MM_DAY} mm/day, {ptype_fraction:.3%} of pixels are FZRA")
       make_plot(precip, year, 'JRA', ds.g0_lon_3, ds.g0_lat_2)

       ds.close()
       ds_ptype.close()


# Analysis parameters
MIN_YEAR = 1980
MAX_YEAR = 2019
LON_MIN = 279
LON_MAX = 285
LAT_MIN = 40
LAT_MAX = 44
THRESHOLD_MM_DAY = 1
THRESHOLD_KG_M2_S = THRESHOLD_MM_DAY * (1/86400) * 0.001
print(THRESHOLD_KG_M2_S)

# Run analysis
analyze_CFSR(min_year=MIN_YEAR, max_year=MAX_YEAR)
analyze_20CR(min_year=MIN_YEAR, max_year=MAX_YEAR)
#analyze_ERA5(min_year=MIN_YEAR, max_year=MAX_YEAR)
analyze_JRA(min_year=MIN_YEAR, max_year=MAX_YEAR)