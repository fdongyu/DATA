#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 09:39:17 2021

@author: feng779
"""

"""
https://floodmap.modaps.eosdis.nasa.gov/
https://floodmap.modaps.eosdis.nasa.gov/readme.php
https://floodmap.modaps.eosdis.nasa.gov/Products/080W050N/2021/MWP_2021012_080W050N_3D3OT.tif

MWP_2012009_020E000S_2D2O.tif
MWP: MODIS Water Product (combines both MFW (flood water) and MSW (surface water), raster only)

MWP (MODIS Water Product):
Introduced March 2012. Currently delivered only in geotiff raster format, with the following pixel values:
0 : Insufficient data to make water determination (cloudy, missing images, swath gaps swaths, or bad data values).
1 : No water detected.
2 : Water detected AND coinciding with reference water (e.g., not flood).
3 : Water detected, beyond reference water, so is likely flood.

To display all surface water (eg., MSW), use all pixels >= 2.

Projection: latitude/longitude geographic
Datum: WGS-84
Pixel size: 0.002197 degrees square (approximately equal to 250 m at the equator)

The core MWP raster product is being generated, 
from which flood water (MFW) and surface water (MSW) are easily extracted (MFW: MWP = 3; MSW: MWP >= 2)

"""

"""
read tif:
    https://stackoverflow.com/questions/7569553/working-with-tiffs-import-export-in-python-using-numpy
    https://stackoverflow.com/questions/10376573/reading-tiff-in-python-and-matplotlib-using-gdal

Instruction for processing the data:
    https://geonetcast.wordpress.com/2020/05/15/python-scripts-flood-mapping-products/
    https://gist.github.com/jgomezdans/2636285
"""

import numpy as np
import matplotlib.pyplot as plt                              # Plotting library
import matplotlib.colors                                     # Matplotlib colors
import numpy as np                                           # Scientific computing with Python
import cartopy, cartopy.crs as ccrs                          # Plot maps
import cartopy.io.shapereader as shpreader                   # Import shapefiles
import time as t                                             # Time access and conversion
import os                                                    # Miscellaneous operating system interfaces
import sys                                                   # Import the "system specific parameters and functions" module
from mpl_toolkits.axes_grid1.inset_locator import inset_axes # Add a child inset axes to this existing axes.
from matplotlib.colors import LinearSegmentedColormap        # Linear interpolation for color maps
from datetime import datetime, timedelta                     # Library to convert julian day to dd-mm-yyyy
from osgeo import gdal, osr, ogr
#from mpl_toolkits.basemap import Basemap   

import pdb



#tiff_filepath = '/Users/feng779/OneDrive - PNNL/Documents/DATA/MODIS/2014/080W040N/MFW_2014001_080W040N_P14x3D3OT.tif'
#tiff_filepath = '/Users/feng779/OneDrive - PNNL/Documents/DATA/MODIS/080W050N/MWP_2012305_080W050N_2D2OT.tif'
#tiff_filepath = '/Users/feng779/OneDrive - PNNL/Documents/DATA/MODIS/080W040N/MWP_2012305_080W040N_2D2OT.tif'
tiff_filepath = '/Users/feng779/OneDrive - PNNL/Documents/DATA/MODIS/080W050N/MWP_2012311_080W050N_2D2OT.tif'

"""
extent
80W - 70W
40N - 30N
"""

# Define the UK OSNG, see <http://spatialreference.org/ref/epsg/27700/>
epsg_to=4326
wkt_from="+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"
wgs84= osr.SpatialReference ()
wgs84.ImportFromEPSG ( epsg_to )
modis = osr.SpatialReference ()
modis.ImportFromProj4( wkt_from )
tx = osr.CoordinateTransformation ( modis, wgs84 )

# Read the GeoTIFF as array
tf = gdal.Open(tiff_filepath)
data = tf.GetRasterBand(1).ReadAsArray()
print("Size:", data.shape)

data[data<3] = 0
pdb.set_trace()
# Getting the GeoTIFF extent
geoTransform = tf.GetGeoTransform()
x_size = tf.RasterXSize # Raster xsize
y_size = tf.RasterYSize # Raster ysize

# Define a cylindrical projection
projection_opts={'projection':'cyl','resolution':'l'}
# These are the extents in the native raster coordinates
min_lon = geoTransform[0]
max_lon = min_lon + geoTransform[1] * x_size
max_lat = geoTransform[3]
min_lat = max_lat + geoTransform[5] * y_size
# The full GeoTIFF extent
extent = [min_lon, min_lat, max_lon, max_lat]
print("Min Lon:", min_lon, "Min Lat:", min_lat, "Max Lon:", max_lon, "Max Lat:", max_lat)

# The extent I would like to plot
#plot_area = [-90.0, 20, -60, 50]
plot_area = extent

plt.rcParams.update({'font.size': 18})
dpi = 15
#fig = plt.figure(figsize=(data.shape[1]/dpi, data.shape[0]/dpi), dpi=dpi) 
fig = plt.figure(figsize=(12, 12.5)) 
# Use the Geostationary projection in cartopy
ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection=ccrs.PlateCarree())
# The full extent of the GeoTIFF
img_extent = [extent[0], extent[2], extent[1], extent[3]]
# The extent you want to plot (the limit is the full extent above)
ax.set_extent([plot_area[0], plot_area[2], plot_area[1], plot_area[3]], ccrs.PlateCarree())

# COLOR SETUP 1 colormap
#c_nodata = '#000000' # 14       (No Data)
# c_nodata = '#ffffff'
# c_land   = '#c0fcf1' # 17 | 43  (Land)
# c_clouds = '#a8f2ff' # 30 | 99  (Clouds)
# c_supra  = '#a8e2ff' # 56       (Supra Snow/Ice Water)
# c_snow   = '#a8ceff' # 70       (Snow)
# c_ice    = '#96abff' # 84       (Ice)
# c_water  = '#4f61ff' # 99 | 126 (Normal Open Water)
# c_water2 = '#2020fa'
# #Floodwater          # 100~200 | 140~254
# c_floodwater = ["#3dfaab", "#2df863", "#00f700", "#c6f800", "#faf690", "#f7f701", "#f5c400", "#f79029", "#f75e00", "#ff0000"]
# cmap = matplotlib.colors.ListedColormap([c_nodata, c_land, c_clouds, \
#                                          c_supra, c_snow, c_ice, c_water, c_water2])
#cmap_flood = matplotlib.colors.ListedColormap(c_floodwater)
#cmap_flood = 'Blues'
#boundaries = [0, 2, 15, 18, 31, 57, 71, 85, 100]

c_nodata = '#ffffff'
c2 = '#4679fa'
c20 = '#46befa'
c40 = '#61fa46'
c60 = '#e2fa46'
c80 = '#faa346'
c100 = '#fa4c46'
cmap = matplotlib.colors.ListedColormap([c_nodata, c2, c20, c40, c60, c80, c100])
boundaries = [0, 2, 20, 40, 60, 80, 100]

c_nodata = '#ffffff'
c1 = 'b'
c2 = '#ffed66'
c3 = '#ff5959'
cmap = matplotlib.colors.ListedColormap([c_nodata, c1, c3])
boundaries = [0, 1, 3]


print("Plotting the Map Elements...")



# Get only the flood from the full data
data_flood = data.astype(np.float64)

## find flood water, or mask surface water
#data_flood = data.astype(np.float64)
#data_flood[data_flood<=2] = 0

vmin = 0
vmax = 100
norm = matplotlib.colors.BoundaryNorm(boundaries, cmap.N, clip=True)
img1 = ax.imshow(data, origin='upper', extent=img_extent, cmap=cmap, norm=norm, zorder=1)

#print("Plotting the Floodwater fraction...")
#img2 = ax.imshow(data_flood, vmin=vmin, vmax=vmax, origin='upper', extent=img_extent, cmap=cmap_flood, zorder=2)

print("Plotting other elements...")

## colorbar
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="3%", pad=0.05)
cax = inset_axes(ax, width="3%", height="100%", loc='right', borderpad=-1.5)
cb = fig.colorbar(img1, cax=cax, orientation='vertical')
cb.ax.tick_params(labelsize=12)
cb.ax.yaxis.offsetText.set_fontsize(12)
cb.set_label('flood fraction (%)', fontsize=14)

# Add states and provinces
shapefile = list(shpreader.Reader('/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/boundary_lines/ne_10m_admin_1_states_provinces_lines/ne_10m_admin_1_states_provinces_lines.shp').geometries())
#ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='k',facecolor='none', linewidth=1.00, zorder=3)
# Add countries
shapefile = list(shpreader.Reader('/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/boundary_lines/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp').geometries())
#ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black',facecolor='none', linewidth=1.25, zorder=4)
# Add continents
shapefile = list(shpreader.Reader('/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/boundary_lines/ne_10m_coastline/ne_10m_coastline.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black',facecolor='none', linewidth=0.5, zorder=5)
# Add coastlines, borders and gridlines
ax.gridlines(color='k', alpha=0.5, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 1), ylocs=np.arange(-180, 180, 1), draw_labels=False, zorder=6)

## set x and y ticks
#xticks = [ r'{:2d}$\degree$'.format(x) for x in np.arange(-80, -70+1, 1, dtype=np.int64) ] 
xticks = np.arange(-80, -70+1, 1, dtype=np.int64)
yticks = np.arange(30,40+1,1, dtype=np.int64)
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
#plt.savefig('example.png')
#plt.close()

plt.show()




