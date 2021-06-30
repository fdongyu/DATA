#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:42:51 2021

@author: feng779
"""

import os

from getUSGSnwis import getUSGSnwis, getUSGSgh

""" 
USGS site map:
    https://maps.waterdata.usgs.gov/mapper/index.html
"""

## Gauge height 
#stationids = ['01482170']
#stationids = ['01463500']
#stationids = ['01462500']

stationids = ['01463500',  # Delaware
              '01578310',  # Susquehanna
              '01647600',  # Potomac
              '01646500',  # Potomac
              '01594440',  # Patuxent
              '01491000',  # Choptank
              '01668000',  # Rappahannock
              '01673000',  # Pamunkey
              '02037500']  # James





vname = 'gauge_height'

starttime = '1980-01-01'
endtime = '2021-01-01'
      
wdir = os.getcwd() ## work directory
ddir = '/Users/feng779/OneDrive - PNNL/Documents/DATA/USGS/gage_height' ## data directory

if not os.path.isdir(ddir):
    os.makedirs(ddir)

ncfile = os.path.join(ddir, 'USGS_gage_height_hourly.nc')
if os.path.exists(ncfile):
     os.remove(ncfile) #this deletes the file
else:
     print("creating file:", ncfile)

getUSGSgh(stationids,starttime,endtime,ncfile, freq='H') 