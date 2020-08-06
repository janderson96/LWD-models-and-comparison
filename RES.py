# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 11:09:49 2020

@author: Jake
"""
## Aerodynamic Resistance model ##
import numpy as np
import datetime as dt
from CART import *
from RH90 import *
import xarray as xr
import matplotlib.pyplot as plt

def RES(time, U_x, U_y, R_n, T, T_d, precip,tot):
    '''
    RES uses LE function to calculate Latent Heat Energy (LE), iterating over 
    each time step.
    
    Inputs:
        time:   Time series of the data
        U_x:    Zonal wind speed (m/s)
        U_y:    Meridional wind speed (m/s)
        R_n:    Net radiation (Wm^-2)
        T:      Temperature (C)
        T_d:    Dew point temperature (C)
        precip: Precipitation (m)
        tot:    Threshold for water content (mm)
        
    Output:
        wetness: Array of wetness scores (=1 if wet) for each time
        LE_r:    Latent heat flux denisty (Wm^-2)
    '''

    z_c = 2.0 # Crop height (m)
    z_T = 2.0 # Height of input temperature (m)
    diam = 0.1 # Diameter of crop (m)
    
    zom = 0.123 * z_c # Surface Roughness (m)
    d = 0.67 * z_c # Zero plane displacment (m)
    kappa = 0.41 #  Von Kármán constant
    
    # Absolute mean weind speed
    u_zt = np.abs((U_x + U_y) / 2)
    
    # Set up arrays
    wetness = np.zeros(len(time))
    u_star = np.zeros(len(time))
    r_a = np.zeros(len(time))
    e_s = np.zeros(len(time))
    e = np.zeros(len(time))
    s = np.zeros(len(time))
    r_b = np.zeros(len(time))
    LE_r = np.zeros(len(time))
    u_2m = np.zeros(len(time))
    LEmm = np.zeros(len(time))
    precipmm = np.zeros(len(time))
    tot_wet = np.zeros(len(time))
    
    # Convert ERA5 (Jm^-2) into (Wm^-2). Divide by seconds in an hour 
    R_nWm = -R_n / (60 * 60)
    
    for i in range(len(time)):
        
        # Friction velocity [ms^-1]
        u_star[i] = 0.4 * u_zt[i] / (np.log((z_T - 0.65*z_c)/ (0.13 * z_c)))

        # Conversion of 10m wind speed into 2m wind speed [ms^-1]
        # SQ (MARTRE, 2006) and Oke (1987)
        u_2m[i] = (u_star[i] / kappa) * np.log((z_c - d) / zom)
        
        # Aerodynamic resistance [sm^-1]
        r_a[i] = np.log((z_T - (0.65 * z_c)) / (0.13 * z_c))/ (0.4 * u_star[i])
        
        # Convert kPa to hPa --> x 10
        e_s[i] = vapour_pressure(T[i]) * 10
        e[i] = vapour_pressure(T_d[i]) * 10
        
        # Slope of saturation vapour pressure v temperature relationship
        # [hPa C^-1]
        s[i] = (4098 * e_s[i]) / (T[i] + 273.3)**2
        
        # [sm^-1]
        r_b[i] = 314 * (diam / u_2m[i])**0.5
        
        # Confirmed to be [Wm^-2]
        LE_r[i] = - (s[i] * R_nWm[i] + (1200 * (e_s[i] - e[i]))\
            / (r_a[i] + r_b[i])) / (s[i] + 0.64)
        
        # Convert LE [Wm^-2] into [mm m^-2 s^-1]
        # SQ paper. 3600 gets W m^-2 into J m^-2. 
        # 2.45E6 is latent heat of vaporization [J kg^-1]
        LEmm[i] = (LE_r[i] * 3600)  / 2.45E6
        
        # Convert ERA5 precip [m] into [mm]
        precipmm[i] = precip[i] * 1000
        
        # Total depth of water on leaf is LE + precipitation [mm]
        tot_wet[i] = LEmm[i] + precipmm[i]
        # Rao paper, max reservoir is 0.8 mm
        if tot_wet[i] > 0.8:
            tot_wet[i] = 0.8
        
        # If depth of water on leaf is above a threshold (to be determined)
        # then wetness array labels time step as wet (=1)
        if tot_wet[i] >= tot:
            wetness[i] = wetness[i] + 1
         
    
    return wetness, LE_r
        

if __name__ == '__main__':
    
    '''
    Main program script.
    '''
#    fpath = 'RES_precip.nc'
#    data_xr = xr.open_dataset(fpath)
#    
#    start_date = dt.date(1994,3,11)
#    end_date = dt.date(1994,10,31)
#    
#    # To align with PERRYfullCSV.csv data start and end points
#    data_xr_slice = data_xr.sel({'time': slice(start_date,end_date)})
#    
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
#    U_x = data_xr_slice['u10'].values
#    U_y = data_xr_slice['v10'].values
#    T_d = data_xr_slice['d2m'].values - 273.15
#    T = data_xr_slice['t2m'].values - 273.15
#    ssr = data_xr_slice['ssr'].values
#    
#    # Site specific data
#    U_xs = site_data['u10'].values
#    U_ys = site_data['v10'].values
#    T_ds = site_data['d2m'].values - 273.15
#    T_s = site_data['t2m'].values - 273.15
#    ssr_s = site_data['ssr'].values
#    
#    #print('ssr_s =', ssr_s)
#    R_nDur = np.sum(ssr_s > 0)
#    R_nDur1 = np.sum(ssr_s <= 0)
#    #print('R_n above 0 =', R_nDur)
#    #print('R_n below 0 =', R_nDur1)
#    #print('temp =', T_s)
#    wetness, LE_r = RES(datetime_time,U_xs,U_ys,ssr_s,T_s,T_ds)
  

    
    #wetness, LE_r = RES(time,U_x,U_y,R_n,T,T_d)

#time = np.array([0,1,2,3,4])
#U_x = np.array([5,4,7,3,5])
#U_y = np.array([5,3,7,4,5])
#R_n = np.array([25,30,40,35,20])
#T = np.array([5,7,8,10,11])
#T_d = np.array([4,5,6,5,8])

#plt.figure()
#plt.plot(time,LE_r)
#plt.show()