#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 21:52:41 2021

@author: feng779
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors   
import matplotlib.tri as tri
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import gdal
import time
import os

import tilesMODIS

import sys
sys.path.append('../../utils') 
from myearth import findNearset1D, readShpPointLine, readShpPoly

import pdb

class MODIS_tools(object):
    """
    general class for processing MODIS data
    """

    def __init__(self, tiff_filepath, **kwargs):
        self.__dict__.update(kwargs)
        
        self.tiff_filepath = tiff_filepath
        
        if len(self.tiff_filepath) == 1:
            self.lonc, self.latc, self.data = self.readTile(self.tiff_filepath[0])
        elif len(self.tiff_filepath) == 2:
            self.lonc, self.latc, self.data = self.readMultiTile(self.tiff_filepath)
        else:
            print ('No option available, reduce the number of tiles!')
            

    def readTile(self, tilename, subset=True):
        """
        tilename : a single tile    
        """
        
        tf = gdal.Open(tilename)
        data = tf.GetRasterBand(1).ReadAsArray()
        print("Size:", data.shape)
        
        # get flood water only 
        # 3 : Water detected, beyond reference water, so is likely flood.
        data[data<3] = 0
        
        # Getting the GeoTIFF extent
        geoTransform = tf.GetGeoTransform()
        min_lon = int( round(geoTransform[0]) )
        max_lon = int( round(min_lon + geoTransform[1] * tf.RasterXSize) )
        max_lat = int( round(geoTransform[3]) )
        min_lat = int( round(max_lat + geoTransform[5] * tf.RasterYSize) )
        
        lonc = np.zeros([tf.RasterXSize])
        latc = np.zeros([tf.RasterYSize])
        
        lonv = np.zeros([tf.RasterXSize+1])
        latv = np.zeros([tf.RasterYSize+1])
        
        for i in range(tf.RasterXSize):
            lonc[i] = min_lon + (0.5+i) * (max_lon-min_lon)/tf.RasterXSize
        for i in range(tf.RasterYSize):
            latc[i] = min_lat + (0.5+i) * (max_lat-min_lat)/tf.RasterYSize
        for i in range(tf.RasterXSize+1):
            lonv[i] = min_lon + i * (max_lon-min_lon)/tf.RasterXSize
        for i in range(tf.RasterYSize+1):
            latv[i] = min_lat + i * (max_lat-min_lat)/tf.RasterYSize
        
        # need to flip the raster upside down
        data = np.flipud( data )
        
        if subset:
            print ('subsetting the orignal Tile')
            lonc, latc, data = self.subsetTile(lonc, latc, data)
        
        
        lonc, latc = np.meshgrid(lonc, latc)   
        #pdb.set_trace()
        return lonc.ravel(), latc.ravel(), data.ravel()
        
    def readMultiTile(self, tilelist, subset=True):
        
        lonc1, latc1, data1 = self.readTile(tilelist[0], subset=subset)
        lonc2, latc2, data2 = self.readTile(tilelist[1], subset=subset)
        
        lonc = np.concatenate([lonc1.ravel(),lonc2.ravel()])
        latc = np.concatenate([latc1.ravel(),latc2.ravel()])
        data = np.concatenate([data1.ravel(),data2.ravel()])
        
        return lonc, latc, data
        
        
    def subsetTile(self, lonc, latc, data):
        """
        subset a tile for fast processing
        """
        ## bbox of the subset region
        bbox = [-80.58, -74.3, 36.65, 43]
        ind_lon0 = findNearset1D(bbox[0], lonc)
        ind_lon1 = findNearset1D(bbox[1], lonc)
        ind_lat0 = findNearset1D(bbox[2], latc)
        ind_lat1 = findNearset1D(bbox[3], latc)
        
        return lonc[ind_lon0:ind_lon1+1], latc[ind_lat0:ind_lat1+1], \
            data[ind_lat0:ind_lat1+1, ind_lon0:ind_lon1+1]
        
        
        
    def tricontour(self, figname):
        """
        generate the contour plot
        the trigrid configuration is for two tiles
        """
        
        
        c_nodata = '#ffffff'
        c1 = 'b'
        c2 = '#ffed66'
        c3 = '#ff5959'
        #cmap = matplotlib.colors.ListedColormap([c_nodata, c1, c2, c3])
        #bounds = [0, 1, 2, 3]
        cmap = matplotlib.colors.ListedColormap([c_nodata, c1, c3])
        bounds = [0, 1, 3]
        
        
        plt.rcParams.update({'font.size': 18}) 
        #fig = plt.figure(figsize=(18,10))
        fig = plt.figure(figsize=(12,10))
        ax = fig.add_subplot(111)
        #cs = ax.contourf(self.lonc, self.latc, data, bounds, cmap=cmap, vmin=bounds[0], vmax=bounds[-1])
         
        ## tricontour
        triang = tri.Triangulation(self.lonc, self.latc)     
        cs = ax.tricontourf(triang, self.data, bounds, cmap=cmap, vmin=bounds[0], vmax=bounds[-1], extend='max')
        
        self.topo_lines(ax)
        
        #clip_shp = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/srb_reprojected/srb_reprojected.shp'
        clip_shp = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/Chesapeake_Bay_Watershed_Boundary/Chesapeake_Bay_Watershed_Boundary.shp'
        clip_shp2 = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/drbbnd_reprojected/drb_bnd_polygon_reproj.shp'
        #clip = self.shp_clip(ax, clip_shp)
        clip = self.shp_clip2(ax, clip_shp, clip_shp2)
        for contour in cs.collections:
            contour.set_clip_path(clip)
            
        self.multi_river_basin(ax)
            
        ax.set_xlim([-82, -72])
        ax.set_ylim([36, 44])
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_aspect('equal')
        fig.tight_layout()
        plt.savefig(figname)
        plt.close()
        
        # from mpl_toolkits.axes_grid1 import make_axes_locatable
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size="3%", pad=0.05)
        # cb = fig.colorbar(cs, cax=cax, orientation='vertical')
        # cb.ax.tick_params(labelsize=12)
        # cb.ax.yaxis.offsetText.set_fontsize(12)
        # cb.set_label('Flood water', fontsize=14)
               
        
    def topo_lines(self, ax):
        """
        function that reads topology shapefiles and overlay on the figure
        """

        shpfile = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/boundary_lines/ne_10m_coastline/ne_10m_coastline.shp'
        XY, field = readShpPointLine(shpfile)
        
        for line in XY:
            X = line[:,0]
            Y = line[:,1]
            ax.plot(X, Y, '-k', linewidth=0.1)
            
    def shp_clip(self, ax, clip_shp):
        """
        create clip based on a shapefile 
        """
        #shp = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/srb/srb.shp'
        ## reprojected in QGIS, set coordinate system to WGS84 EPSG(4326)
        ## use reproject layer tool
        ## https://gis.stackexchange.com/questions/35590/reprojecting-vector-layer-in-qgis
        
        XY, field = readShpPoly(clip_shp)
        
        for line in XY:
            X = line[:,0]
            Y = line[:,1]
            ax.plot(X, Y, '-k', linewidth=0.5)
        
        vertices = np.asarray(XY[0])
        codes = [Path.MOVETO] + [Path.LINETO]*(vertices.shape[0]-2) + [Path.CLOSEPOLY]
        
        clip = Path(vertices, codes)
        clip = PathPatch(clip, transform=ax.transData)
        #clip = PathPatch(clip, facecolor = 'white')
        #pdb.set_trace()
        return clip
    
    def shp_clip2(self, ax, clip_shp, clip_shp2):
        """
        create clip based on multi shapefile 
        """
        XY1, field1  = readShpPoly(clip_shp)
        XY2, field2  = readShpPoly(clip_shp2)
        
        XY = [XY1[0], XY2[0]]
        
        X = XY1[0][:,0]
        Y = XY1[0][:,1]
        ax.plot(X, Y, '-k', linewidth=1.0)
            
        X = XY2[0][:,0]
        Y = XY2[0][:,1]
        ax.plot(X, Y, '-k', linewidth=1.0)
        
        vertices = []
        codes = []
        for i in range(2):
            for j in range(len(XY[i])):
                vertices.append((XY[i][j][0], XY[i][j][1]))
                #pdb.set_trace()
            codes += [Path.MOVETO]
            codes += [Path.LINETO] * (len(XY[i]) -2)
            codes += [Path.CLOSEPOLY]
        
        #vertices = np.asarray(XY[0])
        #codes = [Path.MOVETO] + [Path.LINETO]*(vertices.shape[0]-2) + [Path.CLOSEPOLY]
        
        clip = Path(vertices, codes)
        clip = PathPatch(clip, transform=ax.transData)
        #clip = PathPatch(clip, facecolor = 'white')
        #pdb.set_trace()
        return clip
    
    def multi_river_basin(self, ax):
        """
        multiple river basin boundaries
        """
        
        shp_delaware = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/drbbnd_reprojected/drb_bnd_polygon_reproj.shp'
        shp_susquehanna = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/srb_reprojected/srb_reprojected.shp'
        shp_potomac = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Potomac_river_basin/Potomac_river_basin.shp'
        shp_patuxent = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Patuxent_river_basin/Patuxent_river_basin.shp'
        shp_choptank = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Choptank_river_basin/Choptank_river_basin.shp'
        shp_rappahannock = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Rappahannock_river_basin/Rappahannock_river_basin.shp'
        shp_mattaponi_pamunkey = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Mattaponi_Pamunkey_river_basin/Mattaponi_Pamunkey_river_basin.shp'
        shp_james_appomattox = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/James_river_basin/James_river_basin.shp'
        
        shp_susquehanna_bay = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/clipped_bay_boundary/Chesapeake_bay_poly.shp'
        shp_delaware_bay = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/clipped_bay_boundary/Delaware_bay_poly.shp'
    
        self.basin_bound(ax, shp_susquehanna)
        self.basin_bound(ax, shp_potomac)
        self.basin_bound(ax, shp_patuxent)
        self.basin_bound(ax, shp_choptank)
        self.basin_bound(ax, shp_rappahannock)
        self.basin_bound(ax, shp_mattaponi_pamunkey)
        self.basin_bound(ax, shp_james_appomattox)
        
        self.basin_bound(ax, shp_susquehanna_bay)
        self.basin_bound(ax, shp_delaware_bay)
        
    
    def basin_bound(self, ax, shpfile):
        
        XY, field = readShpPoly(shpfile)
        
        for line in XY:
            X = line[:,0]
            Y = line[:,1]
            ax.plot(X, Y, '-k', linewidth=0.3)
        


def plotTile(lonc, latc, lonv, latv, data):
    """
    imshow
    """
    import cartopy, cartopy.crs as ccrs      
    import matplotlib.colors 
    
    extent = [lonv.min(), latv.min(), lonv.max(), latv.max()]
    #pdb.set_trace()
    fig = plt.figure(figsize=(12, 12.5)) 
    ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection=ccrs.PlateCarree())
    img_extent = [extent[0], extent[2], extent[1], extent[3]]
    ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())
    
    c_nodata = '#ffffff'
    c2 = '#4679fa'
    c20 = '#46befa'
    c40 = '#61fa46'
    c60 = '#e2fa46'
    c80 = '#faa346'
    c100 = '#fa4c46'
    cmap = matplotlib.colors.ListedColormap([c_nodata, c2, c20, c40, c60, c80, c100])
    boundaries = [0, 2, 20, 40, 60, 80, 100]

    print("Plotting the Map Elements...")
    #vmin = 0
    #vmax = 100
    norm = matplotlib.colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    img1 = ax.imshow(data, origin='upper', extent=img_extent, cmap=cmap, norm=norm, zorder=1)
    
    return img1
    


def plotTile_patches(lonc, latc, lonv, latv, data):
    """
    plot the tile using matplotlib patcheCollection
    """
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    
    # need to flip the raster upside down
    data = np.flipud( data )
    
    ## delta for calculating cell vertices from cell center
    dlon = lonc[0] - lonc[1] / 2.
    dlat = latc[0] - latc[1] / 2.
    
    #lonp = np.zeros([len(lon), 4])
    #latp = np.zeros([len(lat), 4])
    patches = []
    for i, llon in enumerate(lonc):
        for j, llat in enumerate(latc):
            lon_w = llon - dlon
            lon_e = llon + dlon
            lat_n = llat + dlat
            lat_s = llat - dlat          
            ## counter-clockwise vertices for each cell
            xp = np.asarray([lon_w, lon_e, lon_e, lon_w])
            yp = np.asarray([lat_s, lat_n, lat_n, lat_s])
            patches.append(Polygon(np.vstack([xp, yp]).T))
        
    pdb.set_trace()
    cmap = 'viridis'
    pc = PatchCollection(patches, cmap=cmap)
    pc.set_array(data.flatten())
    pc.set_lw(0.1)
    
    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(111)
    ax.add_collection(pc)
    ax.set_xlim([lonv.min(), lonv.max()])
    ax.set_ylim([latv.min(), latv.max()])
    ax.set_aspect('equal')
    plt.show()
    
    
    
    #pdb.set_trace()
def plotTile_pcolor(lonc, latc, lonv, latv, data):
    """
    make plots using the pcolor function
    """
    
    # need to flip the raster upside down
    data = np.flipud( data )
    
    xv, yv = np.meshgrid(lonv, latv)    
    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(111)
    ax.pcolor(xv, yv, data, cmap='RdBu')
    ax.set_aspect('equal')
    plt.show()
    
    
    
    
if __name__ == "__main__":
    

    from datetime import datetime, timedelta
    #starttime = '2012-10-22'
    #endtime = '2012-10-28'
    starttime = '2012-10-22'
    endtime = '2012-11-13'
    starttime = datetime.strptime(starttime, '%Y-%m-%d')
    endtime = datetime.strptime(endtime, '%Y-%m-%d')
    times_data = [starttime + timedelta(days=x) for x in range(0, (endtime-starttime).days+1)]
    ddir = '/Users/feng779/OneDrive - PNNL/Documents/DATA/MODIS'
    for timei in times_data:
        print (timei)
        yday = timei.timetuple().tm_yday
        filepath1 = os.path.join(ddir, '080W050N/MWP_{}{:03d}_080W050N_2D2OT.tif'.format(timei.year, yday) )
        filepath2 = os.path.join(ddir, '080W040N/MWP_{}{:03d}_080W040N_2D2OT.tif'.format(timei.year, yday) )
        if os.path.exists(filepath1) and os.path.exists(filepath2):
            tiff_filepath = [filepath1, filepath2]
            MT = MODIS_tools(tiff_filepath)
            figname = '/Users/feng779/OneDrive - PNNL/Documents/CODE/DATA/MODIS/figures/{}.png'.format(datetime.strftime(timei, '%Y%m%d'))
            MT.tricontour(figname)
        #pdb.set_trace()
    
    