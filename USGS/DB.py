# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 13:58:26 2018

@author: fengdongyu
"""

import netcdfio
import os

import pdb

def write2db(ncfiles, dbfile, create=True):
    '''
    The function below aims at writing the resulting river inflow data
    and tide data into the database file
    ''' 
    #Update and write the resulting data into the database file
#    create = True
    #dbfile = 'GalvestonObs.db'
    #dire1=os.path.dirname(os.path.abspath(__file__))
#    dire1=os.getcwd()
    
#    folder = dire1+'/DATA/'
#    #pdb.set_trace()
#    ncfiles = [folder+'USGS_Rivers.nc',folder+'TCOONTide.nc']

    if create:
        print 'Creating database: %s'%dbfile
        netcdfio.createObsDB(dbfile)
    
    for nc in ncfiles:
        print 'Inserting metadata from: %s'%nc    
        netcdfio.netcdfObs2DB(nc,dbfile)

    print 'Done.'
#    if os.path.isfile(folder+dbfile):
#        os.remove(folder+dbfile) 
#    copy_file(dbfile, folder)

#ncfiles = ['DATA/NOAATide_8775283.nc', 'DATA/NOAATide_8775237.nc', 'DATA/NOAATide_8775792.nc', \
#            'DATA/NOAATide_8775870.nc', 'DATA/NOAATide_8775241.nc', 'DATA/NOAATide_8775244.nc', \
#            'DATA/NOAATide_8775296.nc', 'DATA/NOAATide_8776139.nc']

#ncfiles = ['DATA/OBS/NOAA_temperature_8775283.nc', 'DATA/OBS/NOAA_temperature_8775237.nc', 'DATA/OBS/NOAA_temperature_8775792.nc', \
#            'DATA/OBS/NOAA_temperature_8775870.nc', 'DATA/OBS/NOAA_temperature_8775241.nc', 'DATA/OBS/NOAA_temperature_8775244.nc', \
#            'DATA/OBS/NOAA_temperature_8775296.nc', 'DATA/OBS/NOAA_temperature_8776139.nc']            
ncfiles = [os.getcwd()+'/DATA/TWDB_rivers.nc']
            
dbfile = 'DATA/CCObs.db'
write2db(ncfiles, dbfile, create=False)
