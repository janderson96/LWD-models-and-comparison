# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 12:16:18 2020

@author: Jake
"""
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import datetime
import matplotlib.pyplot as plt
import pandas as pd

from netCDFread import *

def CRC_threshold(time, CRC, threshold):
    '''
    CRC_threshold inputs ERA5 data of time and Canopy Reserviour Content (CRC)
    and returns a wetness array of length time based off of a threshold
    depth of water on leaf surface.
    Inputs:
        time:      time array of data
        CRC:       Canopy reserviour content (m)
        threshold: Given threshold for canopy (m)
    Output:
        Wetness:   Wetness[i] = 1 if time step is over given threshold and leaf
                   if wet. Wetness[i] = 0 if time step gives lower value than
                   threshold and is dry.
    '''
    
    wetness = np.zeros(len(time), dtype = int)
    
    for i in range(len(time)):
        if CRC[i] >= threshold:
            wetness[i] = wetness[i] + 1
            
   
    return wetness


#if __name__ == '__main__':
#    
#    '''
#    Main program script.
#    '''
#    fpath = 'oregon.nc'
#    data_xr = xr.open_dataset(fpath)
#    
#    start_date = datetime.date(2016,1,1)
#    end_date = datetime.date(2016,12,31)
#    
#    data_xr_slice = data_xr.sel({'time': slice(start_date,end_date)})
#    print(data_xr_slice)
#    
#    site_lat = 51.1872
#    site_lon = 0.9155
#    
#    datetime_time = convert_npdt2dtdt(data_xr_slice)
#    print(type(datetime_time[0]))
#    
#    # Data for just the given site_lat and site_lon
#    site_data = data_xr_slice.sel( {'latitude': site_lat,
#                         'longitude': site_lon},
#                       method='nearest')
#    
#    CRC = site_data['src'].values
#    
#    wetness = CRC_threshold(datetime_time, CRC, 0.0000129)
#    
#    plt.figure(figsize=(12,5))
#    plt.plot(datetime_time,CRC*1000, label = 'Skin reserviour content') 
#    plt.axhline(y=0.0000129*1000, color='r', linestyle='-', label ='Threshold = 0.0129mm')
#    plt.legend()
#    plt.ylabel('Skin reserviour content (mm)')
#    plt.xlabel('Date')
#    
#    
#    ## For comparasion between CRC_threshold and PERRY data ##
#    
#    data = pd.read_csv('PERRYfullCSV2.csv')
#    PERRY_Wet = np.rot90(pd.DataFrame(data, columns = ['Wet']))
#    
#    
#    
#    
#    start = datetime.date(1994,3,11)
#    day_rng = 233
#
#    for i in range(day_rng):
#    
#        start_seg = start + datetime.timedelta(hours = 24*i)
#        end_seg = start + datetime.timedelta(hours = 24*i + 24)
#        
#    lwd_m = wetness # Hourly
#            
#    lwd_p = PERRY_Wet # 12 mins (5 in 1 hour)
#    
##    for j in range(start_seg,end_seg):
##        lwd_m = CRC_threshold(datetime_time, CRC, 0.0000129)
##            
##        lwd_p = PERRY_Wet
#        
#    print(lwd_p)
#    print(len(lwd_p[0]))
#    print(lwd_m)
#    print(len(lwd_m))
#        plt.scatter(lwd_m,lwd_p)
#    for i in range(235):
#        PerryHourWet[i] = np.sum()
#        
        
        
    
    
    
    
    
    
    
    
    