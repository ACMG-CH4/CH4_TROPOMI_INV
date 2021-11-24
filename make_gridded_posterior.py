import datetime
import xarray as xr
import os
import numpy as np

def make_gridded_posterior(posterior_SF_path, state_vector_path, save_path):
    '''
    The IMI code outputs the posterior scaling factors as a vector (in a .nc file).
    HEMCO wants the scaling factors as a gridded product, by latitude/longitude.
    This script uses the posterior vector file and the state vector file to generate a gridded
    version of the posterior scaling factors.

    Arguments
       posterior_SF_path [str] : path to the posterior scaling factors from an inversion
       state_vector_path [str] : path to the state vector file, from which we will take coords information
       save_path         [str] : path where the gridded posterior should be saved

    '''

    # Load state vector and scaling factor data
    statevector = xr.load_dataset(state_vector_path)
    scalefactor = xr.load_dataset(posterior_SF_path)

    # Map the scaling factors to the state vector elements
    nlat = len(statevector['lat'])
    nlon = len(statevector['lon'])
    scfac_arr = np.zeros(statevector['StateVector'].shape)
    for ilat in range(nlat):
        for ilon in range(nlon):
            element_id = int(statevector['StateVector'].values[ilat,ilon])
            scfac_arr[ilat,ilon] = scfac['xhat'].values[element_id-1]

    # Convert to data array
    lat = statevector['lat'].values
    lon = statevector['lon'].values
    scfac_arr = xr.DataArray(scfac_arr, [("lat", list(lat)), ("lon", list(lon))], attrs={'units': "none"})
    
    # Create dataset
    ds_scfac = xr.Dataset({"ScaleFactor": (["lat", "lon"], scfac_arr)},
                          coords={"lon": ("lon", lon), "lat": ("lat", lat)})

    # Add attribute metadata
    ds_scfac.lat.attrs['units'] = 'degrees_north'
    ds_scfac.lat.attrs['long_name'] = 'Latitude'
    ds_scfac.lon.attrs['units'] = 'degrees_east'
    ds_scfac.lon.attrs['long_name'] = 'Longitude'
    ds_scfac.ScaleFactor.attrs['units'] = '1'

    # Create netcdf
    ds_scfac.to_netcdf(save_path)

    print(f'Saved gridded file to {save_path}')

if __name__ == '__main__':
    import sys

    posterior_SF_path = sys.argv[1]
    state_vector_path = sys.argv[2]
    save_path = sys.argv[3]

    make_gridded_posterior(posterior_SF_path, state_vector_path, save_path) 
