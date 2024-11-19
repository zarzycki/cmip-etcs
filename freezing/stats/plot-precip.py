import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def get_time_index(timestamp):
    """
    Convert timestamp to array index assuming 6-hourly data (4x daily)
    Format: YYYY-MM-DD-HH
    Returns index in the yearly array
    """
    dt = pd.to_datetime(timestamp, format='%Y-%m-%d-%H')
    day_of_year = dt.dayofyear - 1
    hour_index = dt.hour // 6
    return (day_of_year * 4) + hour_index

def plot_precipitation_snapshot(timestamp='2001-01-01-06', fzra_directory='./',
                              datasets=['JRA', 'CFSR', 'CR20', 'ERA5']):
    """
    Create panel plots of precipitation for all datasets at a specific timestamp
    with dynamic colorbar scaling and state boundaries
    """
    time_idx = get_time_index(timestamp)
    print(f"Using time index {time_idx} for timestamp {timestamp}")

    year = pd.to_datetime(timestamp, format='%Y-%m-%d-%H').year

    # First pass to get global maximum
    global_max = 0
    all_data = {}

    for dataset in datasets:
        try:
            filename = f'{fzra_directory}/pre-{dataset}-{year}.nc'
            ds = xr.open_dataset(filename)

            if 'PRECT' in ds:
                precip_var = ds['PRECT']
            elif 'pr' in ds:
                precip_var = ds['pr']
            else:
                print(f"Available variables in {dataset}: {list(ds.variables)}")
                continue

            plot_data = precip_var.isel(time=time_idx)
            all_data[dataset] = plot_data
            current_max = np.nanmax(plot_data.values)
            global_max = max(global_max, current_max)
            print(global_max)

        except Exception as e:
            print(f"Error processing {dataset} for max calculation: {str(e)}")

    # Set up the plot grid
    n_datasets = len(datasets)
    n_cols = 3
    n_rows = (n_datasets + n_cols - 1) // n_cols

    fig = plt.figure(figsize=(18, 4*n_rows))

    # Create axes with cartopy projection
    axs = []
    projection = ccrs.PlateCarree()
    for idx in range(n_rows * n_cols):
        if idx < n_datasets:
            ax = fig.add_subplot(n_rows, n_cols, idx + 1, projection=projection)
            axs.append(ax)

    # Set colormap settings using global maximum
    vmin, vmax = 0, global_max
    print(f"Using global maximum precipitation value: {global_max:.2f} mm/day")

    im = None

    # Process each dataset
    for idx, dataset in enumerate(datasets):
        print(f"Processing {dataset}")
        try:
            if dataset not in all_data:
                raise ValueError(f"No data available for {dataset}")

            plot_data = all_data[dataset]

            ax = axs[idx]

            # Add state boundaries and coastlines
            ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor='black')
            ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, edgecolor='black')

            # Create the precipitation plot
            im = ax.pcolormesh(plot_data.lon, plot_data.lat, plot_data,
                             transform=ccrs.PlateCarree(),
                             shading='auto', cmap='gist_ncar',
                             vmin=vmin, vmax=vmax)

            ax.set_title(f'{dataset}\ntime index: {time_idx}')
            ax.set_extent([234, 294, 23, 50], crs=ccrs.PlateCarree())

            # Add gridlines
            ax.gridlines(draw_labels=True, linestyle=':', alpha=0.5)

        except Exception as e:
            print(f"Error processing {dataset}: {str(e)}")
            ax = axs[idx]
            ax.text(0.5, 0.5, f'Error processing {dataset}\n{str(e)}',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(dataset)

    # Only add colorbar if we successfully plotted at least one dataset
    if im is not None:
        plt.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax, label='Precipitation Rate (mm/day)')
        plt.tight_layout(rect=[0, 0, 0.85, 0.95])
    else:
        plt.tight_layout()
        print("Warning: No datasets were successfully processed")

    plt.suptitle(f'Precipitation Comparison for {timestamp} (index {time_idx})', y=1.02)

    # Save the plot
    plt.savefig(f'precip_comparison_{timestamp}.png', bbox_inches='tight', dpi=300)
    print(f"Plot saved as precip_comparison_{timestamp}.png")
    plt.close()

FZRAPATH = '/glade/derecho/scratch/zarzycki/FZRA/precip/'
plot_precipitation_snapshot(timestamp='2005-08-29-12', fzra_directory=FZRAPATH)