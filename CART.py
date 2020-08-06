# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 12:22:17 2020

@author: Jake
"""
import numpy as np

def CART(T, T_d, U_x, U_y, time, DPDt, Ut, RHt): 
    '''
    Based off of Kim et al (2004) model, decision tree determines if 
    condistions are wet (=1) or dry (=0)
    
    Inputs:
        T:      Temperature array (C)
        T_d:    Dew point temperature array (C)
        U_x:    Wind speed in x direction array (m/s)
        U_y:    Wind speed in y direction array (m/s)
        time:   Time array (YYYY,M,D,H)
        
    Output:
        Wetness: For each step in time, dtermines if conditions satisfy 
                 conditions for wetness (=1) or dry (=0)
        DPD:     dewpoint depression array (C)
        RH:      Relative humidity array (%)
        U:       The mean absolute wind speed (m/s)
    '''
    
    # Set up wetness array, initially all zeros, same length as time array
    wetness = np.zeros(len(time), dtype = int)
    
    # Absolute mean wind speed
    U = np.abs((U_x + U_y) / 2)
    
    # Calculate RH using 'vapour_pressure' function which uses Tetens' formula
    RH = 100 * (vapour_pressure(T_d) / vapour_pressure(T))
    
    # Calculate the dew point depression (DPD)
    DPD = [T - np.abs(T_d) for T, T_d in zip(T,T_d)]
    
    # Loop over entire time series
    for i in range(len(time)):
        # Condition 1
        if DPD[i] <= DPDt: 
            wetness[i] = wetness[i] + 1
        else:
            # Condition 2
            if U[i] <= Ut: 
                test1 = 1.6064 * np.sqrt(T[i]) + 0.0036 * T[i]**2 + 0.1531 \
                * RH[i] - 0.4599 * U[i] * DPD[i] - 0.0035 * T[i] * RH[i]
                if test1 >= 14.4674: 
                    wetness[i] = wetness[i] + 1 
                else:
                    wetness[i] == 0
            else:
                # Condition 3
                if RH[i] <= RHt:  # Kim et al uses RH = 87.8 mine: 63
                    wetness[i] == 0
                else:
                    test2 = 0.7921 * np.sqrt(T[i]) + 0.0046 * RH[i] - 2.3889 \
                    * U[i] - 0.0396 * T[i] + 1.0613 * U[i] * DPD[i]
                    if test2 >= 37.0:
                        wetness[i] = wetness[i] + 1
                    else:
                        wetness[i] == 0
    
    return wetness, DPD, RH, U

def vapour_pressure(x):
    '''
    Using Tetens' formula, vapour pressure is calculated
    
    Inputs:
        T:   Input of temperature for x gives the saturated vapour pressure (C)
        T_d: Input of the dewpoint temperature for x gives the vapour 
             pressure (C)
        
    Output:
        e_s: Saturated vapour pressure if given input is T (kPa)
        e:   Vapour pressure if given input is T_d (kPa)
    '''
    
    e = 0.6112 * np.exp((17.67 * x) / (x + 243.5))
    
    return e


#T = np.array([10,12.3,13.6,15.0,10])
#T_d = np.array([8,9,4,5,1])
#time1 = np.array([0,1,2,3,4])
#U_x = np.array([0,5,10,15,0])
#U_y = np.array([9,2,7,5,4])

#wetness, DPD = CART(Ts, Tds, Uxs, Uys, datetime_time)
