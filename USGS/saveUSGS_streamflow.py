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

## Delaware river and Susquehanna river
# stationids = ['01578310', 
#               '01576000',
#               '01570500',
#               '01554000',
#               '01540500',
#               '01538700',
#               '01536500',
#               '01533400',
#               '01531500',
#               '01515000',
#               '01513831',
#               '01513500',
#               '01503000',
#               '01502731',
#               '01502632',
#               '01500500',
#               '01463500',
#               '01446500',
#               '01438500',
#               '01428500',
#               '01427510',
#               '01427207']

stationids = ['01463500', ## Delaware river
              '01446500',
              '01438500',
              '01428500',
              '01427510',
              '01427207',
              '01578310', ## Susquehanna river
              '01576000',
              '01570500',
              '01554000',
              '01540500',
              '01538700',
              '01536500',
              '01533400',
              '01531500',
              '01515000',
              '01513831',
              '01513500',
              '01503000',
              '01502731',
              '01502632',
              '01500500', 
              '01646502', ## Potomac river
              '01638500',
              '01618000',
              '01594440', ## Patuxent river
              '01592500',
              '01591610',
              '01591000', 
              '01491000', ## Choptank
              '01668000', ## Rappahannock
              '01664000',
              '01674500', ## Mattaponi
              '01674000',
              '01673000', ## Pamunkey
              '02037500', ## James
              '02037000',
              '02035000',
              '02029000',
              '02026000',
              '02025500',
              '02024752',
              '02019500',
              '02016500',
              '02041650', ## Appomattox
              '02040892', 
              '02040000', 
              '02039500']



vname = 'discharge'

starttime = '1990-01-01'
endtime = '2021-01-01'
      
wdir = os.getcwd() ## work directory
ddir = '/Users/feng779/OneDrive - PNNL/Documents/DATA/USGS' ## data directory

if not os.path.isdir(ddir):
    os.makedirs(ddir)

ncfile = os.path.join(ddir, 'USGS_streamflow.nc' )
if os.path.exists(ncfile):
     os.remove(ncfile) #this deletes the file
else:
     print("creating file:", ncfile)

getUSGSnwis(stationids,starttime,endtime,ncfile)   

# for staid in stationids:
    
#     #ncfile = os.path.join(ddir, 'USGS_{}.nc'.format(staid) )
#     ncfile = os.path.join(ddir, 'USGS_{}_{}.nc'.format(vname, staid) )
    
#     if os.path.exists(ncfile):
#         os.remove(ncfile) #this deletes the file
#     else:
#         print("creating file:", ncfile)#add this to prevent errors
    
#     if vname == 'discharge':
#         ## dishcarge
#         getUSGSnwis([staid],starttime,endtime,ncfile)   
#     else:
#         print ('Look for other scripts for downloading gauge height data ...')