# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:39:57 2018

@author: fengdongyu
"""
import numpy as np
import calendar

from tide_functions import downloadTide
import netcdfio

import pdb

#starttime = ['2010-06-01', '2010-07-01', '2010-08-01']
#endtime = ['2010-06-30', '2010-07-31', '2010-08-31']
#staid = '8771450'
#
#output = []
#t = []
#
#for i in range(len(starttime)):
#    start = starttime[i]
#    end = endtime[i]
#    output_tem,t_tem,atts = getNOAATide(start, end, staid)
#    output += output_tem[:-1]
#    t += t_tem[:-1]



#starttime = '2017-01-01'
#endtime = '2017-01-31'
#starttime = '2002-07-01'
#endtime = '2002-07-30'

#starttimes = ['2010-06-01', '2010-07-01', '2010-08-01']
#endtimes = ['2010-06-30', '2010-07-31', '2010-08-31']
#starttimes = ['2002-06-01', '2002-07-01', '2002-08-01']
#endtimes = ['2002-06-30', '2002-07-31', '2002-08-31']

#years = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,\
#        2011,2012,2013,2014,2015,2016,2017,2018]

years = [1990,1991,1992,1993,1994,1995,1996,1997,1998,1999, \
	2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, \
	2010,2011,2012,2013,2014,2015,2016,2017,2018]
months = np.asarray(range(12))+1
starttimes = []
endtimes = []
for year in years:
    for month in months:
        start_day = '01' 
        end_day = calendar.monthrange(year, month)[1]
        starttimes.append('%s-%s-%s'%(year, month, start_day))
        endtimes.append('%s-%s-%s'%(year, month, end_day))


#stationids = ['8776604', \
#              '8777812', \
#             '8778490', \
#             '8779280', \
#             '8779770', \
#             '8779748' ]

#stationids = ['8775870']  ## Bob Hall Pier
stationids = ['8779770']  ## Port Isabel

            
for staid in stationids:
    
    ff = open('eta/availiability/%s.txt'%staid, 'w+')
    print 'Reading water elevation data at station %s ...'%staid
    ## write to txt file for future reference
    print >>ff, 'Reading water elevation data at station %s ...'%staid
    data = downloadTide(starttimes, endtimes, staid, ff)    
    #pdb.set_trace()
    ncfile = 'eta/NOAATide_%s.nc'%staid
    globalatts = {'title':'NOAA oceanographic observation data'}
    netcdfio.writePointData2Netcdf(ncfile,data,globalatts)
    print 'end writing tide data...'
    
    
