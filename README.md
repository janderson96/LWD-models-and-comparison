# LWD-models-and-comparison
The python code for the dissertation titled "Simulation of Leaf Surface Wetness Duration with Threshold and Physical Models using ERA5 Reanalysis Data"
To run, the user will need to download all files to the same folder, including the netCDF files containing ERA5 reanalysis data

Scripts:
RH_90_neat

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
    
CRC_threshold_neat

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
    
 
