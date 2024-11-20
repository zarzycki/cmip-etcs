import xarray as xr
import numpy as np
import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap

def ensure_output_dirs():
    """Create output directories if they don't exist"""
    os.makedirs('nc_files', exist_ok=True)
    os.makedirs('img', exist_ok=True)

def calculate_freezing_rain_hours(ptype_file, prec_file, year, dataset):
    """
    Calculate hours of freezing rain for a given year and individual mask statistics.
    """
    # Load the datasets
    ds_ptype = xr.open_dataset(ptype_file)
    ds_prec = xr.open_dataset(prec_file)

    # Get the raw arrays
    freezing_rain_mask = (ds_ptype.PTYPE == 4).values
    precipitation_mask = (ds_prec.PRECT > 10).values

    # Check shapes match
    if freezing_rain_mask.shape != precipitation_mask.shape:
        raise ValueError(f"Shape mismatch for {dataset} year {year}: "
                        f"PTYPE shape {freezing_rain_mask.shape} != "
                        f"PREC shape {precipitation_mask.shape}")

    # Calculate individual statistics
    ptype_hours = (freezing_rain_mask.sum(axis=0) * 6).astype(float)
    precip_hours = (precipitation_mask.sum(axis=0) * 6).astype(float)

    # Calculate combined statistics
    freezing_rain_events = freezing_rain_mask & precipitation_mask
    annual_hours = (freezing_rain_events.sum(axis=0) * 6).astype(float)

    # Convert all to xarray DataArrays
    annual_hours = xr.DataArray(
        annual_hours,
        coords={'lat': ds_ptype.lat, 'lon': ds_ptype.lon},
        dims=['lat', 'lon']
    )

    ptype_hours = xr.DataArray(
        ptype_hours,
        coords={'lat': ds_ptype.lat, 'lon': ds_ptype.lon},
        dims=['lat', 'lon']
    )

    precip_hours = xr.DataArray(
        precip_hours,
        coords={'lat': ds_ptype.lat, 'lon': ds_ptype.lon},
        dims=['lat', 'lon']
    )

    # Add metadata
    annual_hours.attrs['units'] = 'hours'
    annual_hours.attrs['long_name'] = f'Hours of freezing rain in {year} ({dataset})'
    annual_hours.attrs['dataset'] = dataset

    ptype_hours.attrs['units'] = 'hours'
    ptype_hours.attrs['long_name'] = f'Hours of PTYPE=4 in {year} ({dataset})'
    ptype_hours.attrs['dataset'] = dataset

    precip_hours.attrs['units'] = 'hours'
    precip_hours.attrs['long_name'] = f'Hours of precipitation > 1 mm/day in {year} ({dataset})'
    precip_hours.attrs['dataset'] = dataset

    # Print summary statistics
    print(f"\nSummary Statistics for {dataset} {year}:")
    print(f"Mean hours with PTYPE=4: {ptype_hours.mean().values:.2f}")
    print(f"Mean hours with PREC>1: {precip_hours.mean().values:.2f}")
    print(f"Mean hours of freezing rain: {annual_hours.mean().values:.2f}")
    print(f"Max hours of PTYPE=4: {ptype_hours.max().values:.2f}")
    print(f"Max hours of PREC>1: {precip_hours.max().values:.2f}")
    print(f"Max hours of freezing rain: {annual_hours.max().values:.2f}")

    return {
        'freezing_rain': annual_hours,
        'ptype_only': ptype_hours,
        'precip_only': precip_hours
    }

def create_custom_colormap():
    """Create a custom colormap matching the example"""
    colors = [
        '#2A0359',  # Dark purple
        '#284187',  # Navy blue
        '#4F90CD',  # Light blue
        '#50C878',  # Green
        '#FFEB7F',  # Light yellow
        '#FF9B42',  # Orange
        '#8B2500'   # Dark red
    ]
    return LinearSegmentedColormap.from_list('custom_cmap', colors)

def create_fzra_colormap():
    """Create a colormap matching the FZRA flags plot"""
    colors = [
        '#FFFFFF',  # White
        '#FFEB7F',  # Light yellow
        '#FF9B42',  # Orange
        '#FF69B4',  # Pink
        '#8B3A62',  # Purple
        '#000080',  # Navy blue
        '#000000'   # Black
    ]
    return LinearSegmentedColormap.from_list('fzra_cmap', colors)

def plot_climatology(filename, dataset):
    """
    Create map plots of climatological mean hours for all three metrics
    """
    # Load the data
    ds = xr.open_dataset(filename, decode_times=False)

    # Extract arrays and coordinates
    lons = ds['lon'].values
    lats = ds['lat'].values

    # Dictionary of variables to plot with their specific settings
    plot_vars = {
        'climatology_freezing_rain': {
            'title': 'FZRA Hours',
            'levels': np.arange(0, 60, 10),
            'cmap': create_custom_colormap(),
            'cbar_label': 'Hours/year'
        },
        'climatology_ptype': {
            'title': 'FZRA Flags',
            'levels': np.arange(0, 501, 100),  # 0 to 500 in steps of 100
            'cmap': create_fzra_colormap(),
            'cbar_label': 'Freezing rain flags per year'
        },
        'climatology_precip': {
            'title': 'Precipitation >1mm/day Hours',
            'levels': np.arange(0, 2800, 200),
            'cmap': create_custom_colormap(),
            'cbar_label': 'Hours/year'
        }
    }

    # Create mesh grid for plotting
    lon_mesh, lat_mesh = np.meshgrid(lons, lats)

    for var_name, settings in plot_vars.items():
        # Get data for this variable
        data = ds[var_name].values

        # Print diagnostic information
        print(f"\nPlotting {dataset} {var_name}:")
        print("Data shape:", data.shape)
        print("Valid data range:", np.nanmin(data), "to", np.nanmax(data), "hours")
        print("Number of NaNs:", np.sum(np.isnan(data)))

        # Create figure and axis
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())

        # Create the filled contour plot
        plot = plt.contourf(lon_mesh, lat_mesh, data,
                           levels=settings['levels'],
                           cmap=settings['cmap'],
                           extend='max',
                           transform=ccrs.PlateCarree())

        # Add coastlines and state boundaries
        ax.coastlines()
        ax.add_feature(cfeature.STATES, linestyle='-', edgecolor='black', linewidth=0.5)

        # Set extent
        ax.set_extent([-100, -60, 24, 54], crs=ccrs.PlateCarree())

        # Add gridlines
        gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
        gl.top_labels = False
        gl.right_labels = False

        # Format lat/lon labels
        ax.set_xticks(np.arange(-100, -55, 10))
        ax.set_yticks(np.arange(25, 55, 5))
        lon_formatter = plt.matplotlib.ticker.FormatStrFormatter('%.1f°')
        lat_formatter = plt.matplotlib.ticker.FormatStrFormatter('%.1f°')
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)

        # Add colorbar at the bottom
        cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.03])
        cbar = plt.colorbar(plot, cax=cbar_ax, orientation='horizontal')
        cbar.set_label(settings['cbar_label'])

        # Add title
        plt.suptitle(f'{dataset} (1980-2019): {settings["title"]} (Mean)',
                    y=0.95, fontsize=12)

        # Adjust layout and save
        plt.subplots_adjust(bottom=0.15)  # Make room for colorbar
        output_file = os.path.join('img', f'climatology_{var_name}_{dataset}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Plot saved as {output_file}")

def process_dataset(dataset, start_year, end_year, datadir="./"):
    """
    Process a single dataset and create its climatology
    """
    print(f"\nProcessing {dataset} dataset...")

    # Initialize empty dictionaries to store yearly results
    yearly_fzra = {}
    yearly_ptype = {}
    yearly_precip = {}

    for year in range(start_year, end_year + 1):
        ptype_file = f"{datadir}/ptypes/{dataset}-{year}.nc"
        prec_file = f"{datadir}/precip/pre-{dataset}-{year}.nc"

        if not (glob.glob(ptype_file) and glob.glob(prec_file)):
            print(f"Missing files for {dataset} year {year}, skipping...")
            continue

        print(f"Processing {dataset} year {year}...")
        try:
            results = calculate_freezing_rain_hours(ptype_file, prec_file, year, dataset)
            yearly_fzra[year] = results['freezing_rain']
            yearly_ptype[year] = results['ptype_only']
            yearly_precip[year] = results['precip_only']
        except Exception as e:
            print(f"Error processing {dataset} year {year}: {str(e)}")
            continue

    if not yearly_fzra:
        print(f"No valid data found for {dataset}")
        return

    # Combine all years into single datasets for each variable
    concat_fzra = xr.concat(list(yearly_fzra.values()),
                           dim=pd.Index(list(yearly_fzra.keys()), name='year'))
    concat_ptype = xr.concat(list(yearly_ptype.values()),
                            dim=pd.Index(list(yearly_ptype.keys()), name='year'))
    concat_precip = xr.concat(list(yearly_precip.values()),
                             dim=pd.Index(list(yearly_precip.keys()), name='year'))

    # Calculate climatological means
    climo_fzra = concat_fzra.mean(dim='year')
    climo_ptype = concat_ptype.mean(dim='year')
    climo_precip = concat_precip.mean(dim='year')

    # Create output dataset
    ds = xr.Dataset({
        'yearly_freezing_rain': concat_fzra.astype('int64'),
        'yearly_ptype': concat_ptype.astype('int64'),
        'yearly_precip': concat_precip.astype('int64'),
        'climatology_freezing_rain': climo_fzra.astype(float),
        'climatology_ptype': climo_ptype.astype(float),
        'climatology_precip': climo_precip.astype(float)
    })

    # Add global attributes
    ds.attrs['description'] = f'Annual and climatological statistics for {dataset}'
    ds.attrs['creation_date'] = str(pd.Timestamp.now())
    ds.attrs['dataset'] = dataset
    ds.attrs['period'] = f'{start_year}-{end_year}'

    # Save the results
    output_file = os.path.join('nc_files', f'freezing_rain_hours_{dataset}.nc')
    ds.to_netcdf(output_file)
    print(f"Data saved to {output_file}")

    # Print climatological statistics
    print(f"\n{dataset} Climatological Statistics:")
    print(f"Mean hours/year with PTYPE=4: {float(climo_ptype.mean().values):.2f}")
    print(f"Mean hours/year with PREC>1: {float(climo_precip.mean().values):.2f}")
    print(f"Mean hours/year of freezing rain: {float(climo_fzra.mean().values):.2f}")

    # Create plots
    plot_climatology(output_file, dataset)

if __name__ == "__main__":
    FZRAPATH = '/glade/derecho/scratch/zarzycki/FZRA/'
    datasets = ['ERA5', 'CFSR', 'JRA', 'CR20']
    start_year = 1980
    end_year = 2016

    # Make sure relevant dirs exist locally in this folder
    ensure_output_dirs()

    # Process each dataset
    for dataset in datasets:
        process_dataset(dataset, start_year, end_year, datadir=FZRAPATH)