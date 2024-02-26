import os
import numpy as np
import netCDF4

def calculate_rmse(data1, data2):
    return np.sqrt(np.mean((data1 - data2) ** 2))

# Specify the directory containing the NetCDF files
directory = '/glade/u/home/zarzycki/scratch/ETC-track'

# Specify the reference file
reference_file = 'histogram-06-RANGE_9.01-MINLEN_4-MAXGAP_1-MINEND_12.0-MINTIME_24h.nc'

# Load the reference data
ref_ds = netCDF4.Dataset(os.path.join(directory, reference_file))
ref_data = ref_ds.variables['density'][:]

# List to store RMSE values and filenames
rmse_values = []

# Loop over each file in the directory
for filename in os.listdir(directory):
    if filename.startswith('histogram-01-') and filename.endswith('.nc'):
        # Load the current dataset
        ds = netCDF4.Dataset(os.path.join(directory, filename))
        data = ds.variables['density'][:]

        data = data/6.

        # Calculate RMSE
        rmse = calculate_rmse(ref_data, data)

        # Append RMSE and filename to list
        rmse_values.append((rmse, filename))

        print(rmse,filename)

        # Close the dataset
        ds.close()

# Sort list by RMSE
rmse_values.sort()

# Print sorted list of RMSE values and filenames
for rmse, filename in rmse_values:
    print(f"{rmse:.4f} - {filename}")

# Don't forget to close the reference dataset
ref_ds.close()

