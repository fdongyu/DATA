# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:22:56 2018

@author: fengdongyu
"""

import numpy as np
import datetime
import urllib
import string
from lxml import etree
import netcdfio
from othertime import MinutesSince

import pdb


def downloadTide(starttimes,endtimes,staid, ff):
    """ Main function to construct data in a uniform format"""
    lon=-95.3025
    lat=28.9433
    #ID=8772447
    nn='Galveston Freeport'
    vv='waterlevel'
    # Build up the output data as a list of dictionaries
    meta=[]
    coords = [{'Name':'longitude','Value':lon,'units':'degrees East'},\
        {'Name':'latitude','Value':lat,'units':'degrees North'},\
        {'Name':'time','Value':[],'units':'minutes since 1970-01-01 00:00:00'}]
    attribs = {'StationID':str(staid),'StationName':nn,'Data':[],'coordinates':'time, longitude, latitude','coords':coords} 
    meta.append(attribs)

    data=[]
    tmp = {vv:meta[0]}
    ctr=0

    output = []
    t = []
    for i in range(len(starttimes)):
        starttime = starttimes[i]
        endtime = endtimes[i]
        output_tem,t_tem,atts,coors = getNOAATide(starttime,endtime,staid, ff)
        output += output_tem[:-1]
        t += t_tem[:-1]
        
    ## Update coordination and station info
    tmp[vv]['coords'][0]['Value'] = coors['lon']
    tmp[vv]['coords'][1]['Value'] = coors['lat']        
    if ctr==0:                
        tmp[vv].update(atts)
        
    # Append the data to the list array
    tmp[vv]['Data'] += output
    # Append the time data
    ctr=-1
    for cc in tmp[vv]['coords']:
        ctr+=1
        if cc['Name']=='time':
            tmp[vv]['coords'][ctr]['Value'] += t 
    if np.size(tmp[vv]['Data']) > 0:
            data.append(tmp)
    
    #pdb.set_trace()        
    return data
    
def getNOAATide(starttime, endtime, staid, ff):

    attribs = {'long_name':'Water surface elevation','units':'m'}    
    timeformat = "%Y-%m-%d"
    start = datetime.datetime.strptime(starttime, timeformat)
    end = datetime.datetime.strptime(endtime, timeformat) + datetime.timedelta(days=1)
    
    start = datetime.datetime.strftime(start, "%Y%m%d")
    end = datetime.datetime.strftime(end, "%Y%m%d")
    
    url_water_level = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date=%s 00:00&end_date=%s 00:00&station=%s&product=water_level&datum=msl&units=metric&time_zone=gmt&application=web_services&format=xml'%(start, end, staid)
    url_predictions = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date=%s 00:00&end_date=%s 00:00&station=%s&product=predictions&datum=msl&units=metric&time_zone=gmt&application=web_services&format=xml'%(start, end, staid)
    
    
    fp = urllib.urlopen(url_water_level)
    tree = etree.parse(fp)
    root = tree.getroot()
    fp.close()
    
    
    att = root[0].attrib
    attribs = {'StationID':att['id'],'StationName':att['name']}
    coors = {'lon': string.atof(att['lon']), 'lat':string.atof(att['lat'])}
    
    #pdb.set_trace()
    
#    print root[1][0].attrib['t']
    if len(root[1]) != 0:      
        print 'Available observation data from %s to %s ...'%(str(root[1][0].attrib['t']), str(root[1][-1].attrib['t']))
        print >>ff, 'Available observation data from %s to %s ...'%(str(root[1][0].attrib['t']), str(root[1][-1].attrib['t']))
        time = []
        water_level = []
        for i in range(len(root[1])):
            time.append(datetime.datetime.strptime(root[1][i].attrib['t'], "%Y-%m-%d %H:%M"))
            water_level.append(root[1][i].attrib['v'])
        
        water_level = np.asarray(water_level)
        water_level[water_level==''] = '0'
        ele = []
        for i in range(len(water_level)):
            ele.append(string.atof(water_level[i]))
        
        time = np.asarray(time)
        ele = np.asarray(ele)
        #    #### predictions ####
        #    fp2 = urllib.urlopen(url_predictions)
        #    tree2 = etree.parse(fp2)
        #    root2 = tree2.getroot()
        #    fp2.close()
        #        
        #    print root2[0].attrib['t']
        #    
        #    time2 = []
        #    predictions = []
        #    for i in range(len(root2)):
        #        time2.append(datetime.datetime.strptime(root2[i].attrib['t'], "%Y-%m-%d %H:%M"))
        #        predictions.append(root2[i].attrib['v'])
        #        
        #    predictions = np.asarray(predictions)
        #    predictions[predictions==''] = '0'
        #    ele2 = []
        #    for i in range(len(predictions)):
        #        ele2.append(string.atof(predictions[i])) 
        #
        #    time2 = np.asarray(time2)
        #    ele2 = np.asarray(ele2)    
        #    
        #    nn = len(ele)
        #    ele2[:nn] = ele
        #    
        #    tt=[] #time
        #    for kk in range(len(time2)):
        #        outTime=parseTime2(time2[kk])
        #        tt.append(outTime)
        #    pdb.set_trace()
        tt = []
        for kk in range(len(time)):
            outTime=parseTime2(time[kk])
            tt.append(outTime)
    
    else:
        print "observation data unavailable ..."
        print >>ff, "observation data unavailable ..."
        tt = []
        ele = np.asarray([])
    
    return ele.tolist(), tt, attribs, coors
    
def parseTime2(inTime):

    return MinutesSince(inTime)[0]  