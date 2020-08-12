# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 10:36:41 2020
@author: Jake
"""
## To run for different sites, note locations of double hashtags ##
## Change site specific data name AND
## Threshold varibale in model_data AND
## valiation_data varibale name

import numpy as np
import xarray as xr
import pandas as pa
import datetime as dt
import matplotlib.pyplot as plt
from scipy.stats.stats import linregress

from RH90 import *
from CART import *

# Data from oregon_trial_data.csv
df = pa.read_csv('leafwetness2016.csv')
gdf = df.groupby('site')

# Site specific data
moose_mount_df = gdf.get_group('MMO')
falls_creek_df = gdf.get_group('FCO')
soapgrass_mount_df = gdf.get_group('SGO')
woods_creek_df = gdf.get_group('WCF')

# Moose Mountain LW average data and array of equal length for wetness
moose_mount_LW = np.rot90(pa.DataFrame(moose_mount_df, columns = \
                                       ['LW_average']))
moose = list(moose_mount_LW[0])
moose_wet = np.zeros_like(moose)

# Falls Creek LW average data and array of equal length for wetness
falls_creek_LW = np.rot90(pa.DataFrame(falls_creek_df, columns = \
                                       ['LW_average']))
falls = list(falls_creek_LW[0])
falls_wet = np.zeros_like(falls)

# Soapgrass Mountain LW average data and array of equal length for wetness
soapgrass_mount_LW = np.rot90(pa.DataFrame(soapgrass_mount_df, columns = \
                                           ['LW_average']))
soapgrass = list(soapgrass_mount_LW[0])
soapgrass_wet = np.zeros_like(soapgrass)

# Woods Creek LW average data and array of equal length for wetness
woods_creek_LW = np.rot90(pa.DataFrame(woods_creek_df, columns = \
                                       ['LW_average']))
woods = list(woods_creek_LW[0])
woods_wet = np.zeros_like(woods)

# For all sites, exceding 280mV determines wetness
for i in range(len(moose)):
    if moose[i] >= 280:
        moose_wet[i] = 1
    else:
        moose_wet[i] = 0
        
for j in range(len(falls)):
    if falls[j] >= 280:
        falls_wet[j] = 1
    else:
        falls_wet[j] = 0
        
for k in range(len(soapgrass)):
    if soapgrass[k] >= 280:
        soapgrass_wet[k] = 1
    else:
        soapgrass_wet[i] = 0
        
for l in range(len(woods)):
    if woods[l] >= 280:
        woods_wet[l] = 1
    else:
        woods_wet[l] = 0
      
print('moose LWD duration (h) =', np.sum(moose_wet==1))
print('falls LWD duration (h) =', np.sum(falls_wet==1))
print('soap LWD duration  (h) =', np.sum(soapgrass_wet==1))
print('woods LWD duration (h) =', np.sum(woods_wet==1))

################### Oregon ERA5 data ##################

fpath = 'oregon.nc'
data_xr = xr.open_dataset(fpath)

start_date = dt.date(2016,1,1)
end_date = dt.date(2016,12,31)

# To align with PERRYfullCSV.csv data start and end points
data_xr_slice = data_xr.sel({'time': slice(start_date,end_date)})
print(data_xr_slice)

# Site coordinates 
moose_lat = 44.41469
moose_lon = -122.3943

woods_lat = 44.5335
woods_lon = -123.5498

falls_lat = 44.3956
falls_lon = -122.37374

soap_lat = 44.34796
soap_lon = -122.032680

datetime_time = convert_npdt2dtdt(data_xr_slice)
print(type(datetime_time[0]))

# Data for just the given site_lat and site_lon
moose_data = data_xr_slice.sel( {'latitude': moose_lat,
                     'longitude': moose_lon},
                   method='nearest')

woods_data = data_xr_slice.sel( {'latitude': woods_lat,
                     'longitude': woods_lon},
                   method='nearest')

falls_data = data_xr_slice.sel( {'latitude': falls_lat,
                     'longitude': falls_lon},
                   method='nearest')

soap_data = data_xr_slice.sel( {'latitude': soap_lat,
                     'longitude': soap_lon},
                   method='nearest')
# Site specific data 
## Change variable 'soap_data' to 'moose_data', 'falls_data' etc ##
U_xs = falls_data['u10'].values
U_ys = falls_data['v10'].values
T_ds = falls_data['d2m'].values - 273.15
T_s  = falls_data['t2m'].values - 273.15

## Final 3 variables in model_data are DPD, wind speed and RH thresholds ##
## Threshold values (DPD, wind speed, RH) ##
## Moose Mountain = (3.50,0.25,70)
## Woods Creek = (3.75,0.25,70)
## Falls Creek = (3.00,0.25,70)
## Soapgrass Mountain = (2.50,0.25,70)
model_data, DPD, RH, U = CART(T_s, T_ds, U_xs, U_ys, datetime_time,3,0.25,70)

############ Comparison ############
model_time = datetime_time
validation_time = datetime_time

## Change 'soapgrass_wet' to 'moose_wet', 'falls_wet' or 'woods_wet' ##
validation_data = falls_wet

model_dump = []
validation_dump = []

# loop through days
for i in range(366):
    
    # create daily sum
    model_sum = []
    validation_sum = []
    
    # loop through through hours of the day
    for j in range(24):
        
        # create segments, set both to start at start_date
        start_seg = dt.datetime(2016,1,1,0,0) + dt.timedelta(hours=(i*24)+j )
        end_seg = dt.datetime(2016,1,1,0,0) + dt.timedelta(hours=(i*24)+j + 1)        
        
        # creat index of where the time array is between the start and end seg
        model_ind = np.where(((model_time >= start_seg) & \
                              (model_time < end_seg)) == True)
        validation_ind = np.where(((validation_time >= start_seg) & \
                                   (validation_time < end_seg)) == True)
        
        # slice the data to the indexes and sum it up
        model_hour_sum = np.sum(model_data[model_ind])
        # if there is no data then make the sum a nan so it is not counted
        # len(model_ind[0]) will be 0 if there is no data between segs
        if len(model_ind[0]) == 0:
            model_hour_sum = np.nan
            
        validation_hour_sum = np.sum(validation_data[validation_ind])
        if len(validation_ind[0]) == 0:
            validation_hour_sum = np.nan
            
        model_sum.append(model_hour_sum)
        validation_sum.append(validation_hour_sum)
    
    model_dump.append(np.nansum(model_sum))
    validation_dump.append(np.nansum(validation_sum))
    
err = np.abs(np.array([validation_dump]) - np.array([model_dump]))
err_av = np.mean(err)
sum1 = np.sum(model_dump)
print('Average daily error =', err_av)
print('Sum of model LWD =', sum1)
print('Sum of valid LWD =', np.sum(validation_dump))

slope, intercept, r_value, p_value, std_err = \
linregress(model_dump,validation_dump)
print('slope =', slope)
print('intercept =', intercept)

counter1 = np.zeros_like(model_dump)
counter2 = np.zeros_like(model_dump)
counter5 = np.zeros_like(model_dump)

# Counters for number of model LWD that are between 1, 2 and 5 hours of 
# validation LWD
for i in range(len(model_dump)):
    if (model_dump[i]<=validation_dump[i]+1 and \
        model_dump[i]>=validation_dump[i]-1):
        counter1[i] = 1
        
    if (model_dump[i]<=validation_dump[i]+2 and \
        model_dump[i]>=validation_dump[i]-2):
        counter2[i] = 1
        
    if (model_dump[i]<=validation_dump[i]+5 and \
        model_dump[i]>=validation_dump[i]-5):
        counter5[i] = 1


print('Number of points within 1 hour of daily LWD =', np.sum(counter1==1))
print('Percentage of points within 1 hour =', np.sum(counter1==1)*100/366)
print('Number of points within 2 hours of daily LWD =', np.sum(counter2==1))
print('Percentage of points within 2 hours =', np.sum(counter2==1)*100/366)
print('Number of points within 5 hours of daily LWD =', np.sum(counter5==1))
print('Percentage of points within 5 hours =', np.sum(counter5==1)*100/366)

print('r^2 =', r_value**2)

x_ax = np.linspace(-1,27,5)

plt.figure(figsize=(8,8))
plt.scatter(model_dump,validation_dump)
plt.grid()
plt.plot((-0.5,24),(-0.5,24),c='k')
plt.plot(x_ax, x_ax+1,':b')
plt.plot(x_ax, x_ax-1,':b')
plt.plot(x_ax, x_ax+2, '-.g')
plt.plot(x_ax, x_ax-2, '-.g')
plt.plot(x_ax, x_ax+5,'--m')
plt.plot(x_ax, x_ax-5,'--m')
plt.ylabel('Observed Daily LWD (hours)')
plt.xlabel('Predicted Daily LWD (hours)')
plt.title('CART model daily LWD for given site')
plt.xlim(-0.5,25)
plt.ylim(-0.5,25)