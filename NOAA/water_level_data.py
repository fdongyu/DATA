#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 15:08:00 2021

@author: feng779
"""

import numpy as np
import os
import calendar

from NOAA_water_level import downloadTide

import sys
sys.path.append('/Users/feng779/OneDrive - PNNL/Documents/CODE/Utils')

import netcdfio

import pdb

# years = [1990,1991,1992,1993,1994,1995,1996,1997,1998,1999, \
#  	2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, \
#  	2010,2011,2012,2013,2014,2015,2016,2017,2018]
    
years = list(range(1980, 2021) )

months = np.asarray(range(12))+1
starttimes = []
endtimes = []
for year in years:
    for month in months:
        start_day = '01' 
        end_day = calendar.monthrange(year, month)[1]
        starttimes.append('%s-%s-%s'%(year, month, start_day))
        endtimes.append('%s-%s-%s'%(year, month, end_day))


#stationids = ['8775870']  ## Bob Hall Pier
#stationids = ['8574680'] ## Baltimore

## Baltimore Tolchester Beach, Annapolis, Cambridge, Solomons Island, Washington D.C., Lewisetta
#stationids = ['8574680', '8573364', '8575512', '8571892', '8577330', '8594900', '8635750', '8638610']
stationids = ['8632200']
wdir = os.getcwd() ## work directory
ddir = '/Users/feng779/OneDrive - PNNL/Documents/DATA/NOAA'

for staid in stationids:
    
    if not os.path.isdir(os.path.join(wdir, 'log')):
        os.makedirs(os.path.join(wdir, 'log'))
    #ff = open('log/water_level_station_%s.txt'%staid, 'w+')
    ff = 'log/water_level_station_%s.txt'%staid
    print ('Reading water elevation data at station %s ...'%staid)
    ## write to txt file for future reference
    print ('Reading water elevation data at station %s ...'%staid, file=open(ff, "a"))
    data = downloadTide(starttimes, endtimes, staid, ff)    
    #pdb.set_trace()
    
    if not os.path.isdir(os.path.join(ddir, 'water_level')):
        os.makedirs(os.path.join(ddir, 'water_level'))
    ncfile = os.path.join(ddir, 'water_level/NOAA_water_level_%s.nc'%staid)
    
    if os.path.exists(ncfile):
        os.remove(ncfile) #this deletes the file
    else:
        print("creating file:", ncfile)#add this to prevent errors
    
    globalatts = {'title':'NOAA oceanographic observation data'}
    netcdfio.writePointData2Netcdf_new(ncfile,data,globalatts)
    print ('end writing tide data...')
    
    