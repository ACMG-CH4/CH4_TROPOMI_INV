import xarray as xr
import matplotlib.pyplot as plt
import math
import numpy as np
from os.path import join
from os import listdir
import sys

def make_cluster_file(hemco_diag_pth, land_cover_pth, save_pth, lat_min, lat_max, lon_min, lon_max, emis_threshold=0.1, land_threshold=0.25):
    '''
    Generates the cluster file for an analytical inversion.
    
    Arguments
        hemco_diag_pth [str]   : Path to HEMCO diagnostics file containing prior emissions data
        land_cover_pth [str]   : Path to land cover file
        save_pth       [str]   : Where to save the cluster file
        lat_min        [float] : Minimum latitude 
        lat_max        [float] : Maximum latitude 
        lon_min        [float] : Minimum longitude 
        lon_max        [float] : Maximum longitude
        emis_threshold [float] : Minimum prior emission value to include pixel as cluster [Mg/km2/yr]
        land_threshold [float] : Minimum land fraction to include pixel as cluster

    Returns
        ds_clusters    []      : xarray dataset containing clusters field formatted for HEMCO

    Notes
        - What should we choose for the HEMCO diagnostics default? Should it come from the
          spinup run?
        - Land cover file default can be something like 'GEOSFP.20200101.CN.025x03125.NA.nc'
    '''

    # Load prior emissions data
    emis = xr.open_dataset(hemco_diag_pth)

    # Average over time
    if 'time' in emis.dims:
        emis = emis.mean(dim='time')

    # Separate out anthropogenic methane emissions
    emis['EmisCH4_Anthro'] = (emis['EmisCH4_OtherAnth']  + emis['EmisCH4_Rice'] +
                              emis['EmisCH4_Wastewater'] + emis['EmisCH4_Coal'] +
                              emis['EmisCH4_Landfills']  + emis['EmisCH4_Gas']  +
                              emis['EmisCH4_Livestock']  + emis['EmisCH4_Oil'])
    emis = emis['EmisCH4_Anthro'].copy()

    # Adjust units to Mg/km2/yr
    emis *= 0.001*60*60*24*365*1000*1000

    # Load land cover data
    lc = xr.open_dataset(land_cover_pth)

    # Group together
    lc = (lc['FRLAKE'] + lc['FRLAND'] + lc['FRLANDIC']).drop('time').squeeze()

    if emis_threshold:
        # Where emissions are above threshold, replace with NaN (to be filled later)
        emis = emis.where(emis < emis_threshold)
    
    if land_threshold:
        # Where there is land, replace with NaN (to be filled later)
        emis = emis.where(lc < land_threshold)

    if not emis_threshold and not land_threshold:
        # Replace all values with NaN (to be filled later)
        emis = emis.where(emis == -9999.)
        
    # Anywhere that is not NaN, replace with 0
    emis = emis.where(emis.isnull(), 0)

    # Fill in the NaNs with cluster values
    da_clusters = emis.copy()
    da_clusters.values[da_clusters.isnull()] = np.arange(1, da_clusters.isnull().sum()+1)[::-1]

    # Make dataset
    ds_clusters = da_clusters.to_dataset()
    ds_clusters = ds_clusters.rename({'EmisCH4_Anthro':'Clusters'})

    # Add attribute metadata
    ds_clusters.lat.attrs['units'] = 'degrees_north'
    ds_clusters.lat.attrs['long_name'] = 'Latitude'
    ds_clusters.lon.attrs['units'] = 'degrees_east'
    ds_clusters.lon.attrs['long_name'] = 'Longitude'
    ds_clusters.Clusters.attrs['units'] = 'none'

    # Save
    ds_clusters.to_netcdf(save_pth)
    
    return ds_clusters