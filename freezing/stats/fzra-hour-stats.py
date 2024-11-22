import xarray as xr
import numpy as np
import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

def ensure_output_dirs():
    """Create output directories if they don't exist"""
    os.makedirs('nc_files', exist_ok=True)
    os.makedirs('img', exist_ok=True)

def calculate_freezing_rain_hours(ptype_file, prec_file, year, dataset):
    """
    Calculate hours of freezing rain and rates for a given year and individual mask statistics.
    """
    ds_ptype = xr.open_dataset(ptype_file)
    ds_prec = xr.open_dataset(prec_file)

    freezing_rain_mask = (ds_ptype.PTYPE == 4).values
    precipitation = ds_prec.PRECT.values

    if freezing_rain_mask.shape != precipitation.shape:
        raise ValueError(f"Shape mismatch for {dataset} year {year}")

    # Calculate existing metrics
    ptype_hours = (freezing_rain_mask.sum(axis=0) * 6).astype(float)
    precip_mask = (precipitation > 1)
    precip_hours = (precip_mask.sum(axis=0) * 6).astype(float)

    freezing_rain_events = freezing_rain_mask & precip_mask
    annual_hours = (freezing_rain_events.sum(axis=0) * 6).astype(float)

    # Calculate new freezing rain rate metric
    freezing_rain_rate = np.where(freezing_rain_mask, precipitation, 0).mean(axis=0)

    # Convert to DataArrays
    results = {}
    for name, data in [
        ('freezing_rain', annual_hours),
        ('ptype_only', ptype_hours),
        ('precip_only', precip_hours),
        ('fzra_rate', freezing_rain_rate)
    ]:
        results[name] = xr.DataArray(
            data,
            coords={'lat': ds_ptype.lat, 'lon': ds_ptype.lon},
            dims=['lat', 'lon']
        )

    # Add metadata
    results['fzra_rate'].attrs['units'] = 'mm/day'
    results['fzra_rate'].attrs['long_name'] = f'Mean freezing rain rate in {year} ({dataset})'
    results['fzra_rate'].attrs['dataset'] = dataset

    print(f"\nSummary Statistics for {dataset} {year}:")
    print(f"Mean hours with PTYPE=4: {results['ptype_only'].mean().values:.2f}")
    print(f"Mean hours with PREC>1: {results['precip_only'].mean().values:.2f}")
    print(f"Mean hours of freezing rain: {results['freezing_rain'].mean().values:.2f}")
    print(f"Mean freezing rain rate: {results['fzra_rate'].mean().values:.2f} mm/day")
    print(f"Max hours of PTYPE=4: {results['ptype_only'].max().values:.2f}")
    print(f"Max hours of PREC>1: {results['precip_only'].max().values:.2f}")
    print(f"Max hours of freezing rain: {results['freezing_rain'].max().values:.2f}")
    print(f"Max freezing rain rate: {results['fzra_rate'].max().values:.2f} mm/day")

    return results


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

def manual_fzra_climo_colormap():
    """Create a colormap with exactly discrete colors for FZRA flags"""
    colors = [
        '#ffffff',
        '#ffee18',
        '#ffaf40',
        '#ff7376',
        '#f43eb0',
        '#920fee',
        '#3900fe',
        '#0000df',
        '#000068',
        '#000000',
    ]
    return ListedColormap(colors)

def manual_fzra_hours_colormap():
    """Create a colormap with exactly 7 discrete colors for FZRA flags"""
    colors = [
        '#2A1334',  # Dark navy
        '#3690FD',  # Bright blue
        '#3DF779',  # Lime green
        '#DDD833',  # Golden yellow
        '#EE5115',  # Red-orange
        '#6F0604'
    ]
    return ListedColormap(colors)


def plot_climatology(filename, dataset):
    """Create map plots of climatological mean hours for all metrics"""
    ds = xr.open_dataset(filename, decode_times=False)
    lons, lats = ds['lon'].values, ds['lat'].values



    plot_vars = {
        'climatology_freezing_rain': {
            'title': 'FZRA Hours',
            'levels': np.arange(0, 61, 10),
            'cmap': manual_fzra_hours_colormap(),
            'cbar_label': 'Hours/year'
        },
        'climatology_ptype': {
            'title': 'FZRA Flags',
            'levels': np.arange(0, 550, 50),
            'cmap': manual_fzra_climo_colormap(),
            'cbar_label': 'Freezing rain flags per year'
        },
        'climatology_precip': {
            'title': 'Precipitation Hours',
            'levels': np.arange(1000, 6100, 200),
            'cmap': plt.cm.nipy_spectral,
            'cbar_label': 'Hours/year'
        },
        'climatology_fzra_rate': {
            'title': 'FZRA Rate',
            'levels': np.arange(0.00, 0.11, 0.01),
            'cmap': manual_fzra_climo_colormap(),
            'cbar_label': 'mm/day'
        }
    }

    # Create mesh grid for plotting
    lon_mesh, lat_mesh = np.meshgrid(lons, lats)

    for var_name, settings in plot_vars.items():
        # Get data for this variable
        data = ds[var_name].values

        # Create figure and axis
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())

        # Create the filled contour plot
        print(settings['levels'])
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
        #gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
        #gl.top_labels = False
        #gl.right_labels = False

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
        plt.subplots_adjust(bottom=0.15)
        output_file = os.path.join('img', f'climatology_{var_name}_{dataset}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Plot saved as {output_file}")

def process_dataset(dataset, start_year, end_year, datadir="./"):
    """
    Process a single dataset and create its climatology
    """
    yearly_results = {
        'fzra': {}, 'ptype': {}, 'precip': {}, 'rate': {}
    }

    for year in range(start_year, end_year + 1):
        ptype_file = f"{datadir}/ptypes/{dataset}-{year}.nc"
        prec_file = f"{datadir}/precip/pre-{dataset}-{year}.nc"

        if not (glob.glob(ptype_file) and glob.glob(prec_file)):
            print(f"Missing files for {dataset} year {year}, skipping...")
            continue

        try:
            results = calculate_freezing_rain_hours(ptype_file, prec_file, year, dataset)
            yearly_results['fzra'][year] = results['freezing_rain']
            yearly_results['ptype'][year] = results['ptype_only']
            yearly_results['precip'][year] = results['precip_only']
            yearly_results['rate'][year] = results['fzra_rate']
        except Exception as e:
            print(f"Error processing {dataset} year {year}: {str(e)}")
            continue

    if not yearly_results['fzra']:
        print(f"No valid data found for {dataset}")
        return

    # Combine years and calculate climatological means
    concat_data = {}
    climo_data = {}

    for metric in yearly_results:
        concat_data[metric] = xr.concat(
            list(yearly_results[metric].values()),
            dim=pd.Index(list(yearly_results[metric].keys()), name='year')
        )
        climo_data[metric] = concat_data[metric].mean(dim='year')

    # Create output dataset
    ds = xr.Dataset({
        'yearly_freezing_rain': concat_data['fzra'].astype('int64'),
        'yearly_ptype': concat_data['ptype'].astype('int64'),
        'yearly_precip': concat_data['precip'].astype('int64'),
        'yearly_fzra_rate': concat_data['rate'].astype(float),
        'climatology_freezing_rain': climo_data['fzra'].astype(float),
        'climatology_ptype': climo_data['ptype'].astype(float),
        'climatology_precip': climo_data['precip'].astype(float),
        'climatology_fzra_rate': climo_data['rate'].astype(float)
    })

    # Add metadata and save
    ds.attrs['description'] = f'Annual and climatological statistics for {dataset}'
    ds.attrs['creation_date'] = str(pd.Timestamp.now())
    ds.attrs['dataset'] = dataset
    ds.attrs['period'] = f'{start_year}-{end_year}'

    output_file = os.path.join('nc_files', f'freezing_rain_hours_{dataset}.nc')
    ds.to_netcdf(output_file)
    print(f"Data saved to {output_file}")

    # Print statistics including new rate metric
    print(f"\n{dataset} Climatological Statistics:")
    print(f"Mean freezing rain rate: {float(climo_data['rate'].mean().values):.2f} mm/day")

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