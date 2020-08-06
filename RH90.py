# -*- coding: utf-8 -*-
"""
Created on Tue May 26 12:35:58 2020
@author: Jake
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from netCDF4 import Dataset
import datetime
import xarray as xr
#import netCDF4
from CART import *

def RH90(RH, time, RH_cutoff):
    '''
    Model predicts leaf wetness duration (LWD) based upon the assumption the
    leaf surface is wet for RH > cutoff (%) and dry for RH < cutoff (%).
    Inputs: 
        RH:         Time series of RH values 
        time:       Time series for which RH was measured
        RH_cutoff:  The threshold for cutoff
    Outputs:
        wetness: Defined as '1' if the RH >= threshold (wet) and '0' if the RH 
        < threshold (dry)
    '''
    if len(RH) != len(time):
        print('ERROR: RH and time series are not the same length')
    # Set up a wetness array
    wetness = np.zeros_like(time)
    
    # Loop over each time step
    for i in range(len(RH)):
        if RH[i] >= RH_cutoff:
            wetness[i] = 1
        else:
            wetness[i] = 0
            
    DUR = np.sum(wetness==1)
    
    return wetness, DUR

    
#################################

def convert_npdt2dtdt(in_xarray):
    '''
    Converts numpy datetime into datetime.datetime
    Input:
        in_xarray:    e.g,in_xarray = xr.open_dataset(fpath)
    Output:
        Time:         in dtdt format
    '''
        
    time = np.array([pd.to_datetime(in_xarray['time'].values[i]).to_pydatetime() \
                             for i in range(len(in_xarray['time'].values))])
    
    return time

#######################


#data = pd.read_csv('PERRYfullCSV1.csv')
#PERRY_Temp = pd.DataFrame(data, columns = ['Temp'])
#PERRY_Time = pd.DataFrame(data, columns = ['Date & Time'])
#PERRY_RH = np.rot90(pd.DataFrame(data, columns = ['RH']))
#PERRY_Wet = np.rot90(pd.DataFrame(data, columns = ['Wet']))
#
#T = len(PERRY_Time)
###############################################
#fpath = 'oregon_trial_data.nc'
#data_xr = xr.open_dataset(fpath)
#
#start_date = datetime.date(2016,3,1)
#end_date = datetime.date(2016,6,30)
#print(data_xr)
# To align with PERRYfullCSV.csv data start and end points
#data_xr_slice = data_xr.sel({'time': slice(start_date,end_date)})

#print(data_xr_slice)
#
#site_lat = 51.1872
#site_lon = 0.9155
#
#datetime_time = convert_npdt2dtdt(data_xr_slice)
########print(type(datetime_time[0]))
#
## Data for just the given site_lat and site_lon
#site_data = data_xr_slice.sel( {'latitude': site_lat,
#                     'longitude': site_lon},
#                   method='nearest')
#
## Data direct from netCDF file
#T_d = data_xr_slice['d2m'].values - 273.15
#T = data_xr_slice['t2m'].values - 273.15
#
## Site specific data
#T_ds = site_data['d2m'].values - 273.15
#T_s = site_data['t2m'].values - 273.15
#
#e1 = vapour_pressure(T_ds)
#e_s1 = vapour_pressure(T_s)
#RHERA = 100 * (e1 / e_s1)
#print('RH_ERA =', RHERA)
##RH = 100 * (vapour_pressure(T_d) / vapour_pressure(T))
#
#PERRY_RH_Data = list(PERRY_RH[0])
#PERRY_WET_Data = list(PERRY_Wet[0])
#T1 = range(len(PERRY_RH_Data))
#
#mod_wetness, Mod_DUR = RH90(RHERA,datetime_time,83)
#
#plt.figure(1)
#plt.plot(T1,PERRY_RH_Data, label = 'Perry')
#plt.ylabel('Relative Humidity (%)')
#plt.xlabel('Date')
#plt.title('Recorded RH at Perry Site')
#plt.show()
#
#plt.figure(2)
#plt.plot(datetime_time,RHERA, label = 'Model')
#plt.ylabel('Relative Humidity (%)')
#plt.xlabel('Date')
#plt.title('Model Relative Humidity from ERA5')
#plt.show()
#
#print('Max ERA RH =', np.max(RHERA))
#print('Min ERA RH =', np.min(RHERA))
#print('Max Perry RH =', np.max(PERRY_RH_Data))
#print('Min Perry RH =', np.min(PERRY_RH_Data))
#print('Mean Perry =', np.mean(PERRY_RH_Data))
#print('Mean ERA =', np.mean(RHERA))
#
#
#
#
#Pure_DUR = np.sum(PERRY_Wet==1)
#print('perry data duration =', Pure_DUR)
#print('mod dur / pure dur =', Mod_DUR / Pure_DUR)
#
#RH_cutoff = np.array([78,79,80,81,82,83,84,85])
#Mod_Pure = np.array([1.2879,1.2361,1.1772,1.1166,1.0506,0.9823,0.9076,0.83235])
#closeness = np.absolute(Mod_Pure - 1) 
#
#
#plt.figure()
#plt.plot(RH_cutoff, closeness)
#plt.xlabel('RH cutoff value (%)')
#plt.ylabel('Wetness duration (model / data)')
#plt.axhline(y=0, color='r', linestyle='-')
#plt.ylim(-0.05,0.35)
#plt.show()
#




#T = range(len(PERRY_Time))
#print(T)
#RH90(PERRY_RH,T)
#print(PERRY_Time)
#datetime_time = convert_npdt2dtdt(data)
#print(type(datetime_time[0]))
#
#date = pd.to_datetime(PERRY_Time, 10-3-1994)
#
#dates = plt.dates.date2num(PERRY_Time)
#plt.plot_date(dates, values)
#
#Ttime = pd.DataFrame(data, columns = ['Date & Time'])
#print(Ttime)
#T = range(len(Ttime))
#
#plt.figure()
#plt.plot(T,PERRY_Temp)
#plt.show()
#
#rootgrp = Dataset("oregon_trial_data.nc", "w", format="NETCDF4")
#print(rootgrp.data_model)
#dset = netCDF4.Dataset('oregon_trial_data.nc')
#dset.variables.keys()


