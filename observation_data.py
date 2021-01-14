# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 17:02:52 2018

@author: fengdongyu
"""

"""
observation data: conductivity, water temperature
"""

import numpy as np
import calendar

from obs_functions import downloadObs
import netcdfio

import pdb

years = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,\
        2011,2012,2013,2014,2015,2016,2017,2018]
months = np.asarray(range(12))+1
starttimes = []
endtimes = []
for year in years:
    for month in months:
        start_day = '01' 
        end_day = calendar.monthrange(year, month)[1]
        starttimes.append('%s-%s-%s'%(year, month, start_day))
        endtimes.append('%s-%s-%s'%(year, month, end_day))

#stationids = ['8775870', \
#            '8775792', \
#            '8775283', \
#            '8775237', \
#            '8775296', \
#            '8776139', \
#            '8775244', \
#            '8775241']
stationids = ['8776604', \
            '8778490', \
            '8779280', \
            '8779770']



vartype = 'temperature'
#vartype = 'salinity'

for staid in stationids:
    
    ff = open('temp/availiability/%s_%s.txt'%(vartype,staid), 'w+')
    print 'Reading %s data at station %s ...'%(vartype,staid)
    ## write to txt file for future reference
    print >>ff, 'Reading %s data at station %s ...'%(vartype,staid)
    data = downloadObs(starttimes, endtimes, staid, ff, vartype)    
    #pdb.set_trace()
    ncfile = 'temp/NOAA_%s_%s.nc'%(vartype,staid)
    globalatts = {'title':'NOAA oceanographic observation data'}
    netcdfio.writePointData2Netcdf(ncfile,data,globalatts)
    print 'end writing tide data...'
