# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:32:36 2018

@author: dongyu
"""

from WindTool import NARR_TAMU_wind

import pdb

starttime = '2011-12-01-00'
endtime = '2013-05-31-00'

NTW = NARR_TAMU_wind(starttime, endtime)
NTW.writeSUNTANS('LM_NARR_20122013.nc')
