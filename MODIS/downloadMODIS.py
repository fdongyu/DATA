#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 15:56:44 2021

@author: feng779
"""

import os
from datetime import datetime, timedelta
import urllib

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import tilesMODIS

import pdb

def checkDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True


def getMODISData(tiles, starttime, endtime, ddir):
    """
    get MODIS data from a given time 
    """
    ## first check duplication in a list
    if checkDuplicates(tiles):
        raise IOError('Double check tile names, there are duplicates!')

    
    
    starttime = datetime.strptime(starttime, '%Y-%m-%d')
    endtime = datetime.strptime(endtime, '%Y-%m-%d')
    times_data = [starttime + timedelta(days=x) for x in range(0, (endtime-starttime).days+1)]
    
    # log_filepath = os.path.join(ddir, 'log.txt')
    # if os.path.exists(log_filepath):
    #     os.remove(log_filepath)

    # example 
    # 3-day composite
    # https://floodmap.modaps.eosdis.nasa.gov/Products/080W040N/2014/MFW_2014001_080W040N_P14x3D3OT.tif
    # 2-day composite
    # https://floodmap.modaps.eosdis.nasa.gov/Products/080W050N/2012/MWP_2012305_080W050N_2D2OT.tif
    # 1-day composite
    # https://floodmap.modaps.eosdis.nasa.gov/Products/080W050N/2012/MWP_2012305_080W050N_1D1OT.tif
    
    for i, tile in enumerate(tiles):
        ddir = os.path.join(ddir, tile)
        if not os.path.isdir(ddir):
            os.makedirs(ddir)
            
        for timei in times_data:
            yday = timei.timetuple().tm_yday
            ## 3-day 
            #url = 'https://floodmap.modaps.eosdis.nasa.gov/Products/{}/{}/MFW_{}{:03d}_{}_P14x3D3OT.tif'.format(tile, timei.year, timei.year, yday, tile)
            ## 2-day 
            url = 'https://floodmap.modaps.eosdis.nasa.gov/Products/{}/{}/MWP_{}{:03d}_{}_2D2OT.tif'.format(tile, timei.year, timei.year, yday, tile)
            ## 1-day
            #url = 'https://floodmap.modaps.eosdis.nasa.gov/Products/{}/{}/MWP_{}{:03d}_{}_1D1OT.tif'.format(tile, timei.year, timei.year, yday, tile)
            
            
            filepath = os.path.join(ddir, 'MWP_{}{:03d}_{}_2D2OT.tif'.format(timei.year, yday, tile) )
            if os.path.exists(filepath):
                os.remove(filepath) #this deletes the file
            
            try:
                print ("Downloading data: ", url)
                #urllib.urlretrieve(url, filepath)
                with open(filepath, 'wb') as f:
                    f.write(urllib.request.urlopen(url).read())
                
            except:
                print ("Data not available: ", url)
                #print ("Data not available: {}".format(url), file=open(log_filepath, "a"))
            
            
        #pdb.set_trace()
    
if __name__ == "__main__":
    ## hurricane sandy 1-day 
    ## https://floodmap.modaps.eosdis.nasa.gov/getTile.php?location=080W050N&day=303&year=2012&product=1
    
    
    ## Hurricane sandy
    #starttime = '2012-10-29'
    #endtime = '2012-11-13'
    ## year 2012
    starttime = '2012-07-02'
    endtime = '2012-12-31'

    ddir = '/Users/feng779/OneDrive - PNNL/Documents/DATA/MODIS'
    tiles = ['080W050N']
    #tiles = ['080W040N']
    getMODISData(tiles, starttime, endtime, ddir)


