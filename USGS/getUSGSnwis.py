# -*- coding: utf-8 -*-
"""
Download and convert USGS stream gauge data

Created on Thu Nov 15 19:14:03 2012

@author: mrayson
"""

import numpy as np
from datetime import datetime
import urllib.request, urllib.error, urllib.parse

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import sys
sys.path.append('/Users/feng779/OneDrive - PNNL/Documents/CODE/Utils')
import netcdfio
import othertime

import pdb


def getUSGSnwis(stationids,starttime,endtime,ncfile, freq='D'):
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
        meta = readUSGSmeta(sid,meta=meta)
        if freq == 'D':
            tt,dd = readUSGStxt_daily(sid,starttime,endtime)
        elif freq == 'H':
            tt,dd = readUSGStxt(sid,starttime,endtime)
        time.append(tt)
        discharge.append(dd)

    USGS2netcdf(ncfile,meta,time,discharge)


def readUSGStxt(stationid,starttime,endtime):
    """
    Read the daily station data from a web-service

    See here:
        http://waterdata.usgs.gov/nwis/news/?automated_retrieval_info#Examples
    """

    scale_fac = 0.0283168 # Cubic feet to cubic metres conversion

    #target_url = 'http://waterservices.usgs.gov/nwis/dv/?format=rdb&sites=%s&startDT=%s&endDT=%s&parameterCd=00060'%(stationid,starttime,endtime)
    target_url = 'http://nwis.waterdata.usgs.gov/usa/nwis/uv?cb_00060=on&format=rdb&site_no=%s&period=&begin_date=%s&end_date=%s'%(stationid,starttime,endtime)

    try:
        print('Opening: %s'%target_url)
        f = urllib.request.urlopen(target_url)
    except:
        raise Exception('cannot open url:\n%s'%target_url)

    #pdb.set_trace()
    StationID=[]
    time=[]
    discharge=[]
    # Read and convert the raw data
    for s in f:
        line = s.split()
        if line[0].decode()=='USGS' and line[5].decode() != 'Ice':
            StationID.append(line[1].decode())
            #time.append(datetime.strptime(line[2],'%Y-%m-%d'))
            #discharge.append(float(line[3])*scale_fac)
            time.append(datetime.strptime(line[2].decode()+' '+line[3].decode(),'%Y-%m-%d %H:%M'))
            if line[5].decode() == '' or line[5].decode() == 'Eqp':
                discharge.append(0*scale_fac)
            else:
                discharge.append(float(line[5].decode())*scale_fac)
                
		      
    f.close()
    
    return np.asarray(time), np.asarray(discharge)


def readUSGStxt_daily(stationid,starttime,endtime):
    """
    Read the daily station data from a web-service

    See here:
        http://waterdata.usgs.gov/nwis/news/?automated_retrieval_info#Examples
    """

    scale_fac = 0.0283168 # Cubic feet to cubic metres conversion

    target_url = 'http://waterservices.usgs.gov/nwis/dv/?format=rdb&sites=%s&startDT=%s&endDT=%s&parameterCd=00060'%(stationid,starttime,endtime)
    #target_url = 'http://nwis.waterdata.usgs.gov/usa/nwis/uv?cb_00060=on&format=rdb&site_no=%s&period=&begin_date=%s&end_date=%s'%(stationid,starttime,endtime)

    try:
        print('Opening: %s'%target_url)
        f = urllib.request.urlopen(target_url)
    except:
        raise Exception('cannot open url:\n%s'%target_url)

    #pdb.set_trace()
    StationID=[]
    time=[]
    discharge=[]
    # Read and convert the raw data
    for s in f:
        line = s.split()
        if line[0].decode()=='USGS' and line[3].decode() != 'Ice':
            StationID.append(line[1].decode())
            time.append(datetime.strptime(line[2].decode(),'%Y-%m-%d'))
            discharge.append(float(line[3].decode())*scale_fac)
            # time.append(datetime.strptime(line[2].decode()+' '+line[3].decode(),'%Y-%m-%d %H:%M'))
            # if line[5].decode() == '' or line[5].decode() == 'Eqp':
            #     discharge.append(0*scale_fac)
            # else:
            #     discharge.append(float(line[5].decode())*scale_fac)
                
		      
    f.close()
    
    return np.asarray(time), np.asarray(discharge)



def readUSGSmeta(stationid,meta=None):
    """
    Read the station meta data from a web-service

    See here:
        http://waterdata.usgs.gov/nwis/news/?automated_retrieval_info#Examples
    """
    

    target_url = 'http://waterservices.usgs.gov/nwis/site/?format=rdb&sites=%s'%stationid

    try:
        f = urllib.request.urlopen(target_url)
    except:
        raise Exception('cannot open url:\n%s'%target_url)

    if meta==None:
        meta={}
        meta.update({'Station ID':[]})
        meta.update({'StationName':[]})
        meta.update({'Latitude':[]})
        meta.update({'Longitude':[]})

    for s in f:
        line = s.decode().split('\t')
        if line[0]=='USGS':
            meta['Station ID'].append(line[1])
            meta['StationName'].append(line[2])
            meta['Latitude'].append(float(line[4]))
            meta['Longitude'].append(float(line[5]))

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
    ## write to a single netcdf file with multi groups
    netcdfio.writePointData2Netcdf(ncfile,ncdict,globalatts)
    ## write to multi netcdf files
    #netcdfio.writePointData2Netcdf_new(ncfile,ncdict,globalatts)
    # Write to a shape file
    netcdfio.pointNC2shp(ncfile,shpfile)


def getUSGSgh(stationids,starttime,endtime,ncfile, freq='D'):
    """
    Main function for grabbing the USGS gauge height data
    """
    meta={}
    meta.update({'Station ID':[]})
    meta.update({'StationName':[]})
    meta.update({'Latitude':[]})
    meta.update({'Longitude':[]})
    meta.update({'elevation':[]})

    time=[]
    gage_height=[]

    for sid in stationids:
        meta = readUSGSghmeta(sid,meta=meta)
        if freq == 'D':
            tt,dd = readUSGSghtxt(sid,starttime,endtime)
        elif freq == 'H':
            tt,dd = readUSGSghtxt_hourly(sid,starttime,endtime)
        time.append(tt)
        gage_height.append(dd)

    USGSgh2netcdf(ncfile,meta,time,gage_height)

    
def readUSGSghtxt(stationid,starttime,endtime):
    """
    Read the daily station data from a web-service

    See here:
        http://waterdata.usgs.gov/nwis/news/?automated_retrieval_info#Examples
    """

    scale_fac = 0.3048 # feet to metre

    #target_url = 'http://waterservices.usgs.gov/nwis/dv/?format=rdb&sites=%s&startDT=%s&endDT=%s&parameterCd=00060'%(stationid,starttime,endtime)
    #target_url = 'http://nwis.waterdata.usgs.gov/usa/nwis/uv?cb_00060=on&format=rdb&site_no=%s&period=&begin_date=%s&end_date=%s'%(stationid,starttime,endtime)

    #https://waterdata.usgs.gov/nwis/dv?cb_00065=on&format=rdb&site_no=01554000&referred_module=sw&period=&begin_date=2020-02-12&end_date=2021-02-11
    target_url = 'https://waterdata.usgs.gov/nwis/dv?cb_00065=on&format=rdb&site_no=%s&referred_module=sw&period=&begin_date=%s&end_date=%s'%(stationid,starttime,endtime)
    #pdb.set_trace()
    try:
        print('Opening: %s'%target_url)
        f = urllib.request.urlopen(target_url)
    except:
        raise Exception('cannot open url:\n%s'%target_url)

    #pdb.set_trace()
    StationID=[]
    time=[]
    gage_height=[]  ## mean guage height, third column
    # Read and convert the raw data
    for s in f:
        line = s.split()
        #pdb.set_trace()
        if line[0].decode()=='USGS': # and line[5].decode() != 'Ice':
            
            #discharge.append(float(line[3])*scale_fac)
            #time.append(datetime.strptime(line[2].decode()+' '+line[3].decode(),'%Y-%m-%d %H:%M'))
            # if line[5].decode() == '' or line[5].decode() == 'Eqp':
            #     #gauge_height.append(0*scale_fac)
            #     print ('No gauge height data available on')
            # else:
            #     gauge_height.append(float(line[5].decode())*scale_fac)
            #     #except:
            #     #    print (line[5].decode())
            if len(line) == 9 and line[7] != 0:   # only collect mean gauge height data
                StationID.append(line[1].decode())
                time.append(datetime.strptime(line[2].decode(),'%Y-%m-%d'))
                gage_height.append(float(line[7].decode())*scale_fac)
            elif len(line) in [5,7] and line[3] != 0 and stationid != '01482170':
                StationID.append(line[1].decode())
                time.append(datetime.strptime(line[2].decode(),'%Y-%m-%d'))
                gage_height.append(float(line[3].decode())*scale_fac)
            elif len(line) == 7 and line[3] != 0 and stationid == '01482170': ## high tide, low tide gauge height
                StationID.append(line[1].decode())
                time.append(datetime.strptime(line[2].decode(),'%Y-%m-%d'))
                avg = ( float(line[3].decode()) + float(line[5].decode()) ) / 2.
                gage_height.append(avg*scale_fac)
    #pdb.set_trace()
		      
    f.close()

    return np.asarray(time), np.asarray(gage_height)

def readUSGSghtxt_hourly(stationid,starttime,endtime):
    """
    Read the daily station data from a web-service

    See here:
        http://waterdata.usgs.gov/nwis/news/?automated_retrieval_info#Examples
    """

    scale_fac = 0.3048 # feet to metre

    #target_url = 'http://waterservices.usgs.gov/nwis/dv/?format=rdb&sites=%s&startDT=%s&endDT=%s&parameterCd=00060'%(stationid,starttime,endtime)
    #target_url = 'http://nwis.waterdata.usgs.gov/usa/nwis/uv?cb_00060=on&format=rdb&site_no=%s&period=&begin_date=%s&end_date=%s'%(stationid,starttime,endtime)

    #target_url = 'https://waterdata.usgs.gov/nwis/dv?cb_00065=on&format=rdb&site_no=%s&referred_module=sw&period=&begin_date=%s&end_date=%s'%(stationid,starttime,endtime)
    target_url = 'https://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_00065=on&format=rdb&site_no=%s&period=&begin_date=%s&end_date=%s'%(stationid,starttime,endtime)
    
    try:
        print('Opening: %s'%target_url)
        f = urllib.request.urlopen(target_url)
    except:
        raise Exception('cannot open url:\n%s'%target_url)

    StationID=[]
    time=[]
    gage_height=[]  ## mean guage height, third column
    # Read and convert the raw data
    for s in f:
        line = s.split()
        
        if line[0].decode()=='USGS': # and line[5].decode() != 'Ice':
            #pdb.set_trace()
            #discharge.append(float(line[3])*scale_fac)
            #time.append(datetime.strptime(line[2].decode()+' '+line[3].decode(),'%Y-%m-%d %H:%M'))
            # if line[5].decode() == '' or line[5].decode() == 'Eqp':
            #     #gauge_height.append(0*scale_fac)
            #     print ('No gauge height data available on')
            # else:
            #     gauge_height.append(float(line[5].decode())*scale_fac)
            #     #except:
            #     #    print (line[5].decode())
            if len(line) == 7:
                StationID.append(line[1].decode())
                time.append(datetime.strptime(line[2].decode()+' '+line[3].decode(),'%Y-%m-%d %H:%M'))
                gage_height.append(float(line[5].decode())*scale_fac)
		      
    f.close()
    
    return np.asarray(time), np.asarray(gage_height)    
    
    


def readUSGSghmeta(stationid,meta=None):
    """
    Read the station meta data from a web-service

    See here:
        http://waterdata.usgs.gov/nwis/news/?automated_retrieval_info#Examples
    """

    target_url = 'http://waterservices.usgs.gov/nwis/site/?format=rdb&sites=%s'%stationid
    
    try:
        f = urllib.request.urlopen(target_url)
    except:
        raise Exception('cannot open url:\n%s'%target_url)

    if meta==None:
        meta={}
        meta.update({'Station ID':[]})
        meta.update({'StationName':[]})
        meta.update({'Latitude':[]})
        meta.update({'Longitude':[]})
        meta.update({'elevation':[]})

    for s in f:
        line = s.decode().split('\t')
        if line[0]=='USGS':
            meta['Station ID'].append(line[1])
            meta['StationName'].append(line[2])
            meta['Latitude'].append(float(line[4]))
            meta['Longitude'].append(float(line[5]))
            meta['elevation'].append(float(line[8])*0.3048) ## feet to m

    return meta


def USGSgh2netcdf(ncfile,meta,time,gauge_height):
    """
    Convert the USGS files to netcdf4 format
    """


    shpfile = ncfile[:-2]+'shp'

    varname = 'gauge_height'
    longname = 'Gauge height'
    units = 'm'

    ncdict=[]
    ii=-1
    for tt,dd in zip(time,gauge_height):
        ii+=1
        timeout = othertime.MinutesSince(tt,basetime=datetime(1970,1,1))
        #ncdict=netcdfio.createObsDict(varname,longname,units,[dd],[timeout],\
        #              [meta['Latitude'][ii]],[meta['Longitude'][ii]],[0.0],[meta['Station ID'][ii]],[meta['StationName'][ii]],ncdict=ncdict )
        ncdict=netcdfio.createObsDict(varname,longname,units,[dd],[timeout],\
                      [meta['Latitude'][ii]],[meta['Longitude'][ii]],[meta['elevation'][ii]],[meta['Station ID'][ii]],[meta['StationName'][ii]],ncdict=ncdict )

    ## Global atts
    globalatts={'Title':'USGS gauge height data'}
    # Write to a netcdf file
    ## write to a single netcdf file with multi groups
    netcdfio.writePointData2Netcdf(ncfile,ncdict,globalatts)
    ## write to multi netcdf files
    #netcdfio.writePointData2Netcdf_new(ncfile,ncdict,globalatts)
    # Write to a shape file
    netcdfio.pointNC2shp(ncfile,shpfile)


#####
## Example call
#
####
## Input variables
#stationids = ['08066500',\
#            '08078000',\
#            '08067500',\
#            '08076000',\
#            '08075000',\
#            '08074500',\
#            '08076500',\
#            '08073600',\
#            '08075770']
#
#starttime = '2000-01-01'
#endtime = '2012-01-01'
#ncfile = 'C:/Projects/GOMGalveston/DATA/River/USGS_Rivers_20002012.nc'
####
#
#getUSGSnwis(stationids,starttime,endtime,ncfile)
