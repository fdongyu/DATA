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


def downloadObs(starttimes,endtimes,staid, ff, vartype):
    """ Main function to construct data in a uniform format"""
    lon=-95.3025
    lat=28.9433
    #ID=8772447
    nn='Galveston Freeport'
    
    #vv='waterlevel'
    if vartype == 'temperature':
        vv = 'watertemp'
    elif vartype == 'salinity':
        vv = 'watersali'
    
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
        output_tem,t_tem,atts,coors = getNOAAObs(starttime,endtime,staid, ff, vartype)
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
    
def getNOAAObs(starttime, endtime, staid, ff, vartype):
        
    timeformat = "%Y-%m-%d"
    start = datetime.datetime.strptime(starttime, timeformat)
    end = datetime.datetime.strptime(endtime, timeformat) + datetime.timedelta(days=1)
    
    start = datetime.datetime.strftime(start, "%Y%m%d")
    end = datetime.datetime.strftime(end, "%Y%m%d")
    
    if vartype == 'temperature':
        attribs = {'long_name':'Water temperature','units':'m'}
        url_water_obs = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date=%s 00:00&end_date=%s 00:00&station=%s&product=water_temperature&units=metric&time_zone=gmt&application=web_services&format=xml'%(start, end, staid)
    elif vartype == 'salinity':
        attribs = {'long_name':'Salinity','units':'m'}
        url_water_obs = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date=%s 00:00&end_date=%s 00:00&station=%s&product=salinity&units=metric&time_zone=gmt&application=web_services&format=xml'%(start, end, staid)
    
    #pdb.set_trace()
    
    fp = urllib.urlopen(url_water_obs)
    tree = etree.parse(fp)
    root = tree.getroot()
    fp.close()
    
    
    att = root[0].attrib
    attribs = {'StationID':att['id'],'StationName':att['name']}
    coors = {'lon': string.atof(att['lon']), 'lat':string.atof(att['lat'])}
    
    #pdb.set_trace()
    
#    print root[1][0].attrib['t']
    if len(root[1]) != 0:      
        print 'Available %s data from %s to %s ...'%(vartype, str(root[1][0].attrib['t']), str(root[1][-1].attrib['t']))
        print >>ff, 'Available %s data from %s to %s ...'%(vartype, str(root[1][0].attrib['t']), str(root[1][-1].attrib['t']))
        time = []
        water_obs = []
        for i in range(len(root[1])):
            if root[1][i].attrib['v'] != '':
                time.append(datetime.datetime.strptime(root[1][i].attrib['t'], "%Y-%m-%d %H:%M"))
                water_obs.append(root[1][i].attrib['v'])
        
        water_obs = np.asarray(water_obs)
        #water_tem[water_tem==''] = '0'
        obs = []
        for i in range(len(water_obs)):
            obs.append(string.atof(water_obs[i]))
        
        time = np.asarray(time)
        obs = np.asarray(obs)
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
        obs = np.asarray([])
    
    return obs.tolist(), tt, attribs, coors
    
def parseTime2(inTime):

    return MinutesSince(inTime)[0]  