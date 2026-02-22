import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

file_path = "./changnon_2004.nc"
data = xr.open_dataset(file_path)

freezing_rain_hours = data["freezing_rain_hours"]

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_extent([-125, -67, 25, 50], crs=ccrs.PlateCarree())

levels = np.arange(0, 24, 1)
cmap = plt.cm.plasma_r

contour = ax.contourf(freezing_rain_hours.longitude,
                      freezing_rain_hours.latitude,
                      freezing_rain_hours,
                      levels=levels,
                      cmap=cmap,
                      transform=ccrs.PlateCarree())

cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Hours/year')

ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.set_title('Changnon (1932-2001): FZRA Hours')

plt.show()

