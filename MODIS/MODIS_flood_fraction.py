#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 13:35:57 2021

@author: feng779
"""


import numpy as np
import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.path import Path
import time
import gdal

#from MODIS_tools import MODIS_tools

import sys
sys.path.append('../../utils') 
from myearth import findNearset1D, readShpPointLine, readShpPoly

import pdb

class MODIS_flood_fraction(object):
    """
    class calculating MODIS flood fraction in the defined river basins
    """
    
    def __init__(self, timei, tiff_filepath, **kwargs):
        self.__dict__.update(kwargs)
        
        self.timei = timei
        self.tiff_filepath = tiff_filepath
        
        # if len(self.tiff_filepath) == 1:
        #     self.lonc, self.latc, self.data = self.readTile(self.tiff_filepath[0], subset=False)
        # elif len(self.tiff_filepath) == 2:
        #     self.lonc, self.latc, self.data = self.readMultiTile(self.tiff_filepath, subset=False)
        # else:
        #     print ('No option available, reduce the number of tiles!')
        
    def readMODIS(self):
        """
        read MODIS data
        """
        
        if len(self.tiff_filepath) == 1:
            lonc, latc, data = self.readTile(self.tiff_filepath[0], subset=True)
        elif len(self.tiff_filepath) == 2:
            lonc, latc, data = self.readMultiTile(self.tiff_filepath, subset=True)
        else:
            print ('No option available, reduce the number of tiles!')
        
        return lonc, latc, data
        
        
    def flood_fraction(self):
        """
        calculate flood fraction
        """
        clip_shp = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/srb_reprojected/srb_reprojected.shp'
        clip_shp1 = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/Chesapeake_Bay_Watershed_Boundary/Chesapeake_Bay_Watershed_Boundary.shp'
        clip_shp2 = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/drbbnd_reprojected/drb_bnd_polygon_reproj.shp'
        
        #ff = self.shp_clip2(clip_shp, clip_shp1, clip_shp2)
        ff = self.shp_clip_multi()
        
        return ff
    
    
    def shp_clip_multi(self):
        """
        clip based on multiple river basins
        """
        
        shp_chesapeake = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/Chesapeake_Bay_Watershed_Boundary/Chesapeake_Bay_Watershed_Boundary.shp'
        
        shp_delaware = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/drbbnd_reprojected/drb_bnd_polygon_reproj.shp'
        shp_susquehanna = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/srb_reprojected/srb_reprojected.shp'
        shp_potomac = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Potomac_river_basin/Potomac_river_basin.shp'
        shp_patuxent = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Patuxent_river_basin/Patuxent_river_basin.shp'
        shp_choptank = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Choptank_river_basin/Choptank_river_basin.shp'
        shp_rappahannock = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Rappahannock_river_basin/Rappahannock_river_basin.shp'
        shp_mattaponi_pamunkey = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Mattaponi_Pamunkey_river_basin/Mattaponi_Pamunkey_river_basin.shp'
        shp_james_appomattox = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/James_river_basin/James_river_basin.shp'
        
        
        shp_files = [shp_chesapeake, \
                     shp_delaware, shp_susquehanna, shp_potomac, shp_patuxent, \
                     shp_choptank, shp_rappahannock, shp_mattaponi_pamunkey, shp_james_appomattox]
            
        ff = []
        for shp_file in shp_files:
            print (shp_file)
            ff.append( self.shp_clip(shp_file) )
            
        return ff
    
    def shp_clip(self, clip_shp):
        """
        flood fraction in one shapefile
        """
        
        XY0, field0  = readShpPoly(clip_shp) 
        XY = XY0[0]
        
        self.bbox = [ XY[:,0].min(), XY[:,0].max(), XY[:,1].min(), XY[:,1].max() ]  ## update bbox
        lonc_basin, latc_basin, data_basin = self.readMODIS()  ## read MODIS subset data for every basin 
        
        vertices = []
        codes = []
        for j in range(len(XY)):
            vertices.append((XY[j][0], XY[j][1]))
        codes += [Path.MOVETO]
        codes += [Path.LINETO] * (len(XY) -2)
        codes += [Path.CLOSEPOLY]
        
        clip = Path(vertices, codes)
        starttimer = time.time()
        flags = clip.contains_points(np.vstack((lonc_basin, latc_basin)).T)
        endtimer = time.time()
        print ('Time elapsed: ', endtimer-starttimer)
        
        data = data_basin[flags]
        ff = len(data[data>=3]) / len(data)
            
        return ff
        
    
    def shp_clip2(self, clip_shp, clip_shp1, clip_shp2):
        """
        create clip based on multi shapefile 
        """
        
        XY0, field0  = readShpPoly(clip_shp)    # Susquehanna river basin
        XY1, field1  = readShpPoly(clip_shp1)   # Chesapeake Bay basin
        XY2, field2  = readShpPoly(clip_shp2)   # Delaware river basin
        
        ## test plot
        #plt.rcParams.update({'font.size': 18}) 
        #fig = plt.figure(figsize=(12,10))
        #ax = fig.add_subplot(111)
        
        XY = [XY0[0], XY1[0], XY2[0]]
        
        X1 = XY1[0][:,0]
        Y1 = XY1[0][:,1]
        
        #ax.plot(X1, Y1, '-k', linewidth=1.0)
        
        X2 = XY2[0][:,0]
        Y2 = XY2[0][:,1]
        
        #ax.plot(X2, Y2, '-k', linewidth=1.0)
        
        
        ff = []
        for i in range(len(XY)):
            vertices = []
            codes = []
            for j in range(len(XY[i])):
                vertices.append((XY[i][j][0], XY[i][j][1]))
            codes += [Path.MOVETO]
            codes += [Path.LINETO] * (len(XY[i]) -2)
            codes += [Path.CLOSEPOLY]
        
        
            clip = Path(vertices, codes)
            starttimer = time.time()
            flags = clip.contains_points(np.vstack((self.lonc, self.latc)).T)
            endtimer = time.time()
            print ('Time elapsed: ', endtimer-starttimer)
            data = self.data[flags]
            ff_tem = len(data[data>=3]) / len(data) 
            ff.append(ff_tem)
        
        
        print (self.timei)
        print (ff)
        
        return ff
            #ax.plot(self.lonc[flags], self.latc[flags], 'ok')
        #plt.show()
        
    
    def readMultiTile(self, tilelist, subset=True):
        
        lonc1, latc1, data1 = self.readTile(tilelist[0], subset=subset)
        lonc2, latc2, data2 = self.readTile(tilelist[1], subset=subset)
        
        lonc = np.concatenate([lonc1.ravel(),lonc2.ravel()])
        latc = np.concatenate([latc1.ravel(),latc2.ravel()])
        data = np.concatenate([data1.ravel(),data2.ravel()])
        
        return lonc, latc, data
        
        
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
            lonc, latc, data = self.subsetData(lonc, latc, data)
        
        
        lonc, latc = np.meshgrid(lonc, latc)   
        #pdb.set_trace()
        return lonc.ravel(), latc.ravel(), data.ravel()
    
    
    def subsetData(self, lonc, latc, data):
        """
        subset MODIS data based on each basin
        this increases the speed for flood fraction calculation
        """
        ## bbox of the subset region
        #bbox = [-79, -74.3, 38.6, 43]
        ind_lon0 = findNearset1D(self.bbox[0], lonc)
        ind_lon1 = findNearset1D(self.bbox[1], lonc)
        ind_lat0 = findNearset1D(self.bbox[2], latc)
        ind_lat1 = findNearset1D(self.bbox[3], latc)
        
        return lonc[ind_lon0:ind_lon1+1], latc[ind_lat0:ind_lat1+1], \
            data[ind_lat0:ind_lat1+1, ind_lon0:ind_lon1+1]
        
    
        
    def savetxt(self, filepath, tin, ff):
        """
        save flood fraction to a text file for plotting
        """
        fid = open(filepath, 'a')
        #line = 'Date         Chesapeake  Delaware    Susquehanna Potomac     Patuxent    Choptank    Rappahannock  Mattaponi   James\n'
        #fid.write(line)
        tin = datetime.strftime(tin, '%Y-%m-%d')
        line = '{:<}'.format(tin)
        for f in ff:
            f = round(f, 5)
            line += '     {:<6.5f}'.format(f)
        line += '\n'
        fid.write(line)
        fid.close()
        
    
    
if __name__ == "__main__":
    
    ## hurricane sandy
    starttime = '2012-10-12'
    #endtime = '2012-11-15'
    endtime = '2012-11-17'
    #endtime = '2012-12-31'
    ## year 2012
    #starttime = '2012-07-02'
    #endtime = '2012-12-31'
    
    starttime = datetime.strptime(starttime, '%Y-%m-%d')
    endtime = datetime.strptime(endtime, '%Y-%m-%d')
    times_data = [starttime + timedelta(days=x) for x in range(0, (endtime-starttime).days+1)]
    ddir = '/Users/feng779/OneDrive - PNNL/Documents/DATA/MODIS'
    txtpath = 'flood_fraction.txt'
    for timei in times_data:
        print (timei)
        yday = timei.timetuple().tm_yday
        filepath1 = os.path.join(ddir, '080W050N/MWP_{}{:03d}_080W050N_2D2OT.tif'.format(timei.year, yday) )
        filepath2 = os.path.join(ddir, '080W040N/MWP_{}{:03d}_080W040N_2D2OT.tif'.format(timei.year, yday) )
        if os.path.exists(filepath1) and os.path.exists(filepath2):
            tiff_filepath = [filepath1, filepath2]
            
            Mff = MODIS_flood_fraction(timei, tiff_filepath)
            
            ff = Mff.flood_fraction()
            Mff.savetxt(txtpath, timei, ff)
            #pdb.set_trace()
    