# -*- coding: utf-8 -*-
"""
Created on Sun Apr  1 22:20:50 2018

@author: fengdongyu
"""

import numpy as np
from datetime import datetime

import netcdfio
import othertime

import pdb



def getTWDBriver(stationids,name_info, lat_info, lon_info, ncfile):
    """
    Main function for grabbing the data and saving it to a netcdf file
    """
    meta={}
    meta.update({'Station ID':[]})
    meta.update({'StationName':[]})
    meta.update({'Latitude':[]})
    meta.update({'Longitude':[]})
    
    time=[]
    discharge=[]
    
    for sid in stationids:
        meta = readTWDBmeta(sid,name_info, lat_info, lon_info, meta=meta)
        tt,dd = readTWDBtxt(sid,name_info)
        
        ## make up some bad data
        if sid == '08211000':
            starttime = datetime.strptime('2000-01-01-00', '%Y-%m-%d-%H')
            endtime = datetime.strptime('2002-10-01-00', '%Y-%m-%d-%H')
            ind0 = othertime.findNearest(starttime, tt)
            ind1 = othertime.findNearest(endtime, tt)
            dd[ind0:ind1] += 50.
        
        time.append(tt)
        discharge.append(dd)
        #pdb.set_trace()
        
    USGS2netcdf(ncfile,meta,time,discharge)
    

def readTWDBtxt(stationid,name_info):
    """
    Read the station data from TWDB text files
    """  
    
    scale_fac = 0.0283168  # Cubic feet to cubic metres conversion
    scale_fac2 = 1/1.98347 # acre-foot/day to cubic foot/second
    
    filename = 'river_data/%s_95_15.txt'%name_info[stationid]
    
    time = []
    discharge = []    
    f = open(filename, 'r')
    for s in f:
        line = s.split()
        if line[0] != 'year-month-day,inflow_afd':
            #pdb.set_trace()
            time.append(datetime.strptime(line[0][:10], '%Y-%m-%d'))
            discharge.append(float(line[0][11:])*scale_fac*scale_fac2)

    f.close()
    return np.asarray(time), np.asarray(discharge)
    
    

def readTWDBmeta(stationid, name_info, lat_info, lon_info, meta=None):
    
    lat = lat_info[stationid]
    lon = lon_info[stationid]
    name = name_info[stationid]
    
    if meta==None:
        meta={}
        meta.update({'Station ID':[]})
        meta.update({'StationName':[]})
        meta.update({'Latitude':[]})
        meta.update({'Longitude':[]})
    

    meta['Station ID'].append(stationid)
    meta['StationName'].append(name)
    meta['Latitude'].append(lat)
    meta['Longitude'].append(lon)    
    
    
    return meta
    
def USGS2netcdf(ncfile,meta,time,discharge):
    """
    Convert the USGS files to netcdf4 format
    """

    
    shpfile = ncfile[:-2]+'shp'
    
    varname = 'discharge'
    longname = 'Stream Discharge Rate'
    units = 'm3 s-1'
    
    ncdict=[]
    ii=-1
    for tt,dd in zip(time,discharge):
        ii+=1
        timeout = othertime.MinutesSince(tt,basetime=datetime(1970,1,1))
        ncdict=netcdfio.createObsDict(varname,longname,units,[dd],[timeout],\
                      [meta['Latitude'][ii]],[meta['Longitude'][ii]],[0.0],[meta['Station ID'][ii]],[meta['StationName'][ii]],ncdict=ncdict )
                      
    ## Global atts
    globalatts={'Title':'USGS stream gage discharge data'}
    # Write to a netcdf file
    netcdfio.writePointData2Netcdf(ncfile,ncdict,globalatts)
    # Write to a shape file
    netcdfio.pointNC2shp(ncfile,shpfile)


  
    
    
