import argparse
import matplotlib.pyplot as plt
from osgeo import gdal
import numpy as np

parser = argparse.ArgumentParser(description="Process an image file with GDAL and plot using matplotlib.")
parser.add_argument("image_file", type=str, help="Path to the image file")
args = parser.parse_args()

# Open the image with GDAL
image_file = args.image_file  # Get the image filename from command-line argument
dataset = gdal.Open(image_file)

# Read the image data
band = dataset.GetRasterBand(1)
img = band.ReadAsArray()

# Get GCPs
gcps = dataset.GetGCPs()

# Extract pixel and geographic coordinates of GCPs
pixel_x = [gcp.GCPPixel for gcp in gcps]
pixel_y = [gcp.GCPLine for gcp in gcps]
geo_x = [gcp.GCPX for gcp in gcps]
geo_y = [gcp.GCPY for gcp in gcps]

# Plot the image
plt.figure(figsize=(10, 10))
plt.imshow(img, cmap='gray')

# Overlay GCPs as red points
plt.scatter(pixel_x, pixel_y, color='red', marker='x', s=100)
for i, (px, py, gx, gy) in enumerate(zip(pixel_x, pixel_y, geo_x, geo_y)):
    plt.text(px, py, f'GCP {i+1}\n({gx:.2f}, {gy:.2f})', color='yellow', fontsize=12, ha='right')

plt.title("GCPs on Image")
plt.xlabel("Pixel X")
plt.ylabel("Pixel Y")
plt.grid()
plt.show()
