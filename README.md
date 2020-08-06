# LWD-models-and-comparison
The python code for the dissertation titled "Simulation of Leaf Surface Wetness Duration with Threshold and Physical Models using ERA5 Reanalysis Data"
To run, the user will need to download all files to the same folder, including the netCDF files containing ERA5 reanalysis data

Scripts:

RH_90_neat.py

RH_90

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
    
CRC_threshold_neat.py

CRC_threshold

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
    
CART_neat.py

CART

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
    
vapour_pressure

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
    
 
