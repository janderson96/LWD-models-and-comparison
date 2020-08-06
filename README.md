# LWD-models-and-comparison
The python code for the dissertation titled "Simulation of Leaf Surface Wetness Duration with Threshold and Physical Models using ERA5 Reanalysis Data"
To run, the user will need to download all files to the same folder, including the netCDF files containing ERA5 reanalysis data

Scripts:

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
 
