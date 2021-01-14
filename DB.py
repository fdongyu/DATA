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

basedir = os.getcwd()
#filenames = ['eta/NOAATide_8775870.nc', 'eta/NOAATide_8776604.nc', 'eta/NOAATide_8777812.nc', \
#            'eta/NOAATide_8778490.nc', 'eta/NOAATide_8779280.nc', 'eta/NOAATide_8779748.nc', \
#            'eta/NOAATide_8779770.nc']

filenames = ['eta/NOAATide_8776139.nc']

#filenames = ['temp/NOAA_temperature_8776604.nc', 'temp/NOAA_temperature_8778490.nc', 'temp/NOAA_temperature_8779280.nc', \
#            'temp/NOAA_temperature_8779770.nc']            

#filenames = ['river/DATA/TWDB_rivers.nc']

ncfiles = []        
for i in range(len(filenames)):
    filename = '%s/%s'%(basedir, filenames[i])
    ncfiles.append(filename)


dbfile = 'LMObs.db'
write2db(ncfiles, dbfile, create=False)
