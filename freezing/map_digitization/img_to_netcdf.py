from osgeo import gdal
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader
from shapely.geometry import shape
from shapely.ops import unary_union
from shapely import vectorized
import xarray as xr
import os
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("image_file", type=str, help="Path to the georeferenced GeoTIFF file")
args = parser.parse_args()

# Path to the georeferenced GeoTIFF file
geotiff_path = args.image_file  # Get the image filename from command-line argument

# Use PIL to open and convert the image to RGB
with Image.open(geotiff_path) as img:
    # Check if the image has a palette and convert to RGB if necessary
    if img.mode != 'RGB':
        img = img.convert('RGB')

    # Convert the image to a numpy array
    img_array = np.array(img)

# Check for unexpected dimensions (e.g., presence of alpha channel)
print("Image shape:", img_array.shape)  # Should be (height, width, 3)

# Define a dictionary to map RGB colors to freezing rain hours
color_to_value = {
    (124, 189, 57): 0,
    (250, 244, 137): 3,
    (231, 205, 71): 6,
    (252, 166, 42): 9,
    (228, 112, 33): 12,
    (238, 76, 65): 15,
    (241, 151, 192): 18,
    (249, 226, 239): 21
}

# Initialize an array to hold freezing rain values
freezing_rain_hours = np.full(img_array.shape[:2], np.nan)

def color_match(rgb, target_rgb, tolerance=5):
    """Check if the color matches the target color within a given tolerance."""
    return np.all(np.abs(np.array(rgb) - np.array(target_rgb)) <= tolerance)

# was 20
tolerance = 18
# Iterate over the image and map colors to values with tolerance
for color, value in color_to_value.items():
    print(color)
    mask = np.apply_along_axis(lambda rgb: color_match(rgb, color, tolerance=tolerance), -1, img_array)
    freezing_rain_hours[mask] = value

# Check the number of pixels assigned to each category (for debugging)
unique, counts = np.unique(freezing_rain_hours[~np.isnan(freezing_rain_hours)], return_counts=True)
print("Value counts in extracted data:", dict(zip(unique, counts)))

output_directory="./figs/"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Plot each color match to visualize
for color, value in color_to_value.items():
    mask = np.apply_along_axis(lambda rgb: color_match(rgb, color, tolerance), -1, img_array)
    # Plot the result and save as PNG
    plt.figure()
    plt.imshow(mask, cmap='viridis')
    plt.colorbar(label='Freezing Rain Hours')
    plt.title('Extracted Freezing Rain Hours')
    plt.savefig(f"{output_directory}/freezing_rain_hours_{value}.png")
    plt.close()

# Plot the result
plt.figure()
plt.imshow(freezing_rain_hours, cmap='viridis')
plt.colorbar(label='Freezing Rain Hours')
plt.title('Extracted Freezing Rain Hours')
plt.savefig(f"{output_directory}/final_result.png")
plt.close()








# Open the georeferenced image to get geo-transform and projection info
dataset = gdal.Open(geotiff_path, gdal.GA_ReadOnly)
if not dataset:
    raise RuntimeError(f"Unable to open the GeoTIFF file: {geotiff_path}")

# Get the geo-transform and projection information
geotransform = dataset.GetGeoTransform()
projection = dataset.GetProjection()

# Print the geotransform and projection for debugging
print("Geotransform:", geotransform)
print("Projection:", projection)

# Define the corners of the image based on the geotransform
min_x = geotransform[0]
max_y = geotransform[3]
pixel_width = geotransform[1]
pixel_height = geotransform[5]
max_x = min_x + (pixel_width * dataset.RasterXSize)
min_y = max_y + (pixel_height * dataset.RasterYSize)

# Generate arrays of x (longitude) and y (latitude) coordinates
x_coords = np.arange(min_x, max_x, pixel_width)
y_coords = np.arange(max_y, min_y, pixel_height)

# Create a meshgrid of coordinates
lon_grid, lat_grid = np.meshgrid(x_coords, y_coords)

# Interpolate the extracted data to a regular grid using griddata
# Prepare known points and values
known_points = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))
known_values = freezing_rain_hours.ravel()

# Filter out NaN values for interpolation
mask = ~np.isnan(known_values)
known_points = known_points[mask]
known_values = known_values[mask]

# Define the grid for interpolation (same as input grid in this case)
interp_lon_grid, interp_lat_grid = lon_grid, lat_grid

# Interpolate the values
interpolated_data = griddata(
    known_points,
    known_values,
    (interp_lon_grid, interp_lat_grid),
    method='linear'
)

# Create a CONUS mask using cartopy's state boundaries
# Get US state geometries
shapefile_path = shapereader.natural_earth(resolution='10m', category='cultural', name='admin_1_states_provinces')
reader = shapereader.Reader(shapefile_path)
us_states = [shape(state.geometry) for state in reader.records() if state.attributes['admin'] == 'United States of America']

# Combine geometries to create a unified CONUS shape
conus_shape = unary_union(us_states)

# Use vectorized operations to create the CONUS mask
# Vectorized contains checks for each point in the grid
conus_mask = vectorized.contains(conus_shape, lon_grid, lat_grid)

# Path to the Great Lakes shapefile
great_lakes_shapefile = shapereader.natural_earth(resolution='10m', category='physical', name='lakes')

# Read the Great Lakes geometries
reader = shapereader.Reader(great_lakes_shapefile)
great_lakes = [shape(record.geometry) for record in reader.records() if record.attributes['name'] in ['Lake Superior', 'Lake Michigan', 'Lake Huron', 'Lake Erie', 'Lake Ontario']]

# Combine geometries to create a unified Great Lakes shape
great_lakes_shape = unary_union(great_lakes)

# Create a mask for the Great Lakes
great_lakes_mask = vectorized.contains(great_lakes_shape, lon_grid, lat_grid)

# Combine the CONUS and Great Lakes masks (keep only CONUS without Great Lakes)
final_mask = conus_mask & ~great_lakes_mask

# Apply the final mask to the interpolated data
interpolated_data[~final_mask] = np.nan

### Round to discretized data

# Extract the list of target freezing rain hour values
target_values = np.array(list(color_to_value.values()))

# Define a function to round to the nearest target value
def round_to_nearest(value):
    if np.isnan(value):
        return value  # Preserve NaN values
    # Find the closest value in target_values
    nearest_value = target_values[np.argmin(np.abs(target_values - value))]
    return nearest_value

# Apply the rounding function to each cell in the interpolated data
# Use vectorized operation for efficiency
rounded_interpolated_data = np.vectorize(round_to_nearest)(interpolated_data)

interpolated_data = rounded_interpolated_data





# Assuming `interpolated_data` contains the grid data, `lon_grid` and `lat_grid` contain the coordinates

# Define the mask for NaN values
nan_mask = np.isnan(interpolated_data)

# Get coordinates of NaN cells
nan_coords = np.column_stack((lon_grid[nan_mask], lat_grid[nan_mask]))

# Check if there are any NaN cells to process
if nan_coords.size > 0:
    # Build a k-d tree for efficient spatial queries
    tree = cKDTree(nan_coords)

    # Define the maximum distance in degrees (approximation)
    # 1 degree latitude ~ 111 km, so 100 km corresponds to ~0.9 degree latitude
    # For longitude, it depends on the latitude, so this is an approximation
    max_distance_km = 25
    lat_km_per_degree = 111.0
    max_distance_deg = max_distance_km / lat_km_per_degree

    # Flatten the grid coordinates
    grid_coords = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))

    # Query the k-d tree to find distances to the nearest NaN cell
    distances, _ = tree.query(grid_coords, distance_upper_bound=max_distance_deg)

    # Reshape the distances back to the grid shape
    distance_grid = distances.reshape(lon_grid.shape)

    # Apply the mask to cells within the distance threshold to set them to -999
    # Use np.logical_not(nan_mask) to avoid changing NaN cells themselves
    close_to_nan_mask = (distance_grid < max_distance_deg) & np.logical_not(nan_mask)
    interpolated_data[close_to_nan_mask] = -999

# Identify the mask of negative values (e.g., -99999)
negative_mask = (interpolated_data < 0)

# Identify the mask of valid positive values
valid_mask = (interpolated_data >= 0)

# Coordinates of valid positive values
valid_coords = np.column_stack((lon_grid[valid_mask], lat_grid[valid_mask]))
valid_values = interpolated_data[valid_mask]

# Check if there are valid values to use for replacement
if valid_coords.size > 0:
    # Build a k-d tree for valid positive values
    tree = cKDTree(valid_coords)

    # Coordinates of negative value cells
    negative_coords = np.column_stack((lon_grid[negative_mask], lat_grid[negative_mask]))

    # Query the k-d tree to find the nearest valid positive value for each negative cell
    _, nearest_idx = tree.query(negative_coords)

    # Replace each negative cell with the nearest valid positive value
    interpolated_data[negative_mask] = valid_values[nearest_idx]





# Plot the masked interpolated data with state outlines
plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.PlateCarree())  # Use Plate Carree projection for lat/lon data

# Plot the interpolated data
plt.imshow(interpolated_data, extent=[min_x, max_x, min_y, max_y], origin='upper', cmap='viridis', transform=ccrs.PlateCarree())
plt.colorbar(label='Freezing Rain Hours')
plt.title('Interpolated Freezing Rain Hours (Masked to CONUS)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Add state outlines
ax.add_feature(cfeature.STATES, edgecolor='white', linewidth=1)

# Set the extent to match the data bounds
ax.set_extent([min_x, max_x, min_y, max_y], crs=ccrs.PlateCarree())

plt.savefig(f"{output_directory}/freezing_rain_hours_netcdf.png")
plt.close()


# Create xarray Dataset to save as NetCDF
ds = xr.Dataset(
    {
        "freezing_rain_hours": (("latitude", "longitude"), interpolated_data)
    },
    coords={
        "longitude": x_coords,
        "latitude": y_coords,
    }
)

# Add metadata
ds.attrs['title'] = 'Interpolated Freezing Rain Hours (Masked to CONUS)'
ds.attrs['projection'] = projection
ds.attrs['geotransform'] = geotransform

# Specify compression and chunking options
encoding = {
    "freezing_rain_hours": {
        "zlib": True,
        "complevel": 1
    }
}

# Save as NetCDF file
netcdf_output_path='changnon_2004.nc'
ds.to_netcdf(netcdf_output_path, format='NETCDF4', encoding=encoding)
print(f"NetCDF file saved as {netcdf_output_path}")

# Clean up
dataset = None