# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 10:38:44 2018

@author: fengdongyu
"""

from getTWDBriver import getTWDBriver

import pdb


#station_info = {'aransas':08189700, 'cavasso':12100405, 'copano':08189200, \
#                'mission':08189500, 'nueces':08211000, 'oso':08211520}

name_info = {'08189700':'aransas', '12100405':'cavasso', '08189200':'copano', \
                '08189500':'mission', '08211000':'nueces', '08211520':'oso'}
                
##cavasso creek: 28.21444, -96.98306 
            
lat_info = {'08189700':'28.28250426', '12100405':'28.21444', '08189200':'28.30361693', \
            '08189500':'28.29195088', '08211000':'28.03834719', '08211520':'27.71141879'}
            
lon_info = {'08189700':'-97.620829', '12100405':'-96.98306', '08189200':'-97.1124913', \
            '08189500':'-97.2791593', '08211000':'-97.8602769', '08211520':'-97.5019377'}
            

stationids = ['08189700', \
            '12100405', \
            '08189200', \
            '08189500', \
            '08211000', \
            '08211520']
            
ncfile = 'DATA/TWDB_rivers.nc'
getTWDBriver(stationids, name_info, lat_info, lon_info, ncfile)



            
