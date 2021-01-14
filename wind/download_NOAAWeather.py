# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 15:20:49 2018

@author: dongyu
"""

from datetime import datetime
import getNOAAWeatherStation as noaa

## 25.943583, -97.844853
## 27.501967, -97.234824
#latlon = [-97.6,-96.88,27.32,28.24]
latlon = [-97.844853,-97.234824, 25.943583,27.501967]
timestart = datetime.strptime('2000', '%Y') 
timeend = datetime.strptime('2018', '%Y')
dt = 1.0/24.0
localdir = 'rawdata'
showplot = False
ncfile = 'NCDCNWS_AirObs_20002018.nc'
shpfile = 'NCDCNWS_AirObs_20002018.shp'
data = noaa.noaaish2nc(latlon,[timestart.year,timeend.year],localdir,ncfile,shpfile)
