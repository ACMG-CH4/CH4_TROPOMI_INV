import datetime
import xarray as xr
import os
import numpy as np

def do_gridding(vector, clust):
    '''
    Project input vector on the inversion grid using information from the state vector file.
    Input vector should be a numpy array of scale factors (SF), diagonal elements of the posterior 
    error covariance matrix (S_post), or diagonal elements of the averaging kernel matrix (A).
    '''
    
    # Map the vector (e.g., scale factors) to the cluster grid
    nlat = len(clust['lat'])
    nlon = len(clust['lon'])
    target_array = np.zeros(clust['Clusters'].shape)
    for ilat in range(nlat):
        for ilon in range(nlon):
            cluster_id = int(clust['Clusters'].values[ilat,ilon])
            target_array[ilat,ilon] = vector[cluster_id-1]

    # Convert to data array
    lat = clust['lat'].values
    lon = clust['lon'].values
    target_array = xr.DataArray(target_array, 
                                [("lat", list(lat)), ("lon", list(lon))], 
                                attrs={'units': "none"})

    return target_array


def make_gridded_posterior(posterior_SF_path, clusters_path, save_path):
    '''
    The IMI code outputs the inversion results as vectors and matrices (in a .nc file).
    We (and HEMCO, for scale factors) want the results as a gridded product, by latitude/longitude.
    This script uses the inversion results file and the clusters file to generate a gridded
    version of the posterior scaling factors, posterior errors, and averaging kernel sensitivities.

    Arguments
       posterior_SF_path [str] : path to the posterior scaling factors from an inversion
       clusters_path     [str] : path to the cluster file, from which we will take coords information
       save_path         [str] : path where the gridded posterior should be saved

    '''

    # Load clusters and inversion results data
    clust = xr.load_dataset(clusters_path)
    inv_results = xr.load_dataset(posterior_SF_path)

    # Get the scale factors and the diagonals of the S_post and A matrices
    SF = inv_results['xhat'].values
    S_post = np.diagonal(inv_results['S_post'].values)
    A = np.diagonal(inv_results['A'].values)

    # Do gridding
    gridded_SF = do_gridding(SF, clust)
    gridded_S_post = do_gridding(S_post, clust)
    gridded_A = do_gridding(A, clust)

    # Create dataset
    ds = xr.Dataset({'SF_Nonwetland': (["lat", "lon"], gridded_SF),
                     'S_post': (["lat", "lon"], gridded_S_post),
                     'A': (["lat", "lon"], gridded_A)},
                    coords={"lon": ("lon", lon), "lat": ("lat", lat)})

    # Add attribute metadata
    ds.lat.attrs['units'] = 'degrees_north'
    ds.lat.attrs['long_name'] = 'Latitude'
    ds.lon.attrs['units'] = 'degrees_east'
    ds.lon.attrs['long_name'] = 'Longitude'
    ds.SF_Nonwetland.attrs['units'] = '1'
    ds.S_post.attrs['units'] = '1'
    ds.A.attrs['units'] = '1'
    
    # Create netcdf
    ds.to_netcdf(save_path)

    print(f'Saved gridded file to {save_path}')

if __name__ == '__main__':
    import sys

    posterior_SF_path = sys.argv[1]
    clusters_path = sys.argv[2]
    save_path = sys.argv[3]

    make_gridded_posterior(posterior_SF_path, clusters_path, save_path) 
