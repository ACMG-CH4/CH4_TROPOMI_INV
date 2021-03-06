#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import numpy as np
from netCDF4 import Dataset
import xarray as xr
import pickle
import os

# Notes:
# ======
# - Manual step: generate observational error data file "mean_error.nc"
# - Inversion steps at the end are not clear. Need new variable names. Possible to rewrite that
#   section in a more mathematically transparent way, or will that slow things down?

# ==================================================================================================
#
#                                      Define functions
#
# ==================================================================================================

def save_obj(obj, name ):
    """ Save something with Pickle. """

    with open(name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


# --------------------------------------------------------------------------------------------------

def load_obj(name):
    """ Load something with Pickle. """

    with open( name, 'rb') as f:
        return pickle.load(f)


# --------------------------------------------------------------------------------------------------

def nearest_loc(loc_query, loc_grid, tolerance=1):
    """ Find the index of the nearest grid location to a query location, with some tolerance. """

    distances = np.abs(loc_grid - loc_query)
    ind = distances.argmin()
    if distances[ind] >= tolerance:
        return np.nan
    else:
        return ind


# ==================================================================================================
#
#                                      Run the code
#
# ==================================================================================================

# Configuration
n_clust = 1199+8
xlim = [-106.5625,-82.1875];
ylim = [11,35]
gamma = 0.25

# Read observational error data
filename = "/n/holyscratch01/jacob_lab/lshen/CH4/GEOS-Chem/Flexgrid/CPU_mexico_inversion/Step4_inversion_pixel/mean_error.nc"
data = xr.open_dataset(filename)
lon_GC = data['lon'].values
lat_GC = data['lat'].values
mean_error_std = data['error'].values
mean_error = mean_error_std**2                  # Error variance from standard deviation
mean_error = np.einsum('ij->ji', mean_error)
data.close()

# Read output data from Step1 (virtual TROPOMI column XCH4, Jacobian matrix)
os.chdir("/n/holyscratch01/jacob_lab/lshen/CH4/GEOS-Chem/Flexgrid/CPU_mexico_inversion/data_converted/")
files = glob.glob("*.pkl")
files.sort()

# Initialize ________
all_part1 = np.zeros([n_clust,n_clust], dtype=float)
all_part2 = np.zeros([n_clust], dtype=float)

# For each .pkl file from Step1:
for fi in files:
    
    # Load TROPOMI/GEOS-Chem and Jacobian matrix data from the .pkl file
    print(fi)
    met = load_obj(fi)
    # If there aren't any TROPOMI observations on this day, skip 
    # [****Shouldn't we have no files for those cases anyway?]
    if met['obs_GC'].shape[0] == 0:
        continue
    # Otherwise, grab the TROPOMI/GEOS-Chem data
    obs_GC = met['obs_GC']
    # Only consider data within latitude and longitude bounds
    ind = np.where((obs_GC[:,2]>=xlim[0]) & (obs_GC[:,2]<=xlim[1]) & (obs_GC[:,3]>=ylim[0]) & (obs_GC[:,3]<=ylim[1]))
    if (len(ind[0]) == 0):          # Skip if no data in bounds
        continue
    obs_GC = obs_GC[ind[0],:]       # TROPOMI and GEOS-Chem data within bounds
    KK = 1e9 * met['KK'][ind[0],:]  # Jacobian entries for observations within bounds [ppb]
    NN = obs_GC.shape[0]            # Number of observations

    # Now lower the sensitivity to BC by 50%
    # ****This makes the k-means clusters less sensitive
    #KK[:,1199:] = KK[:,1199:]/2
    
    # Initialize observation error matrix diagonal entries
    obs_error = np.zeros((NN,))
    # For each TROPOMI observation:
    for iNN in range(NN):
        # Get the closest GC pixel by latitude and longitude
        iGC = nearest_loc(obs_GC[iNN,2], lon_GC)
        jGC = nearest_loc(obs_GC[iNN,3], lat_GC)
        # Get the observational error for that location
        obs_error[iNN] = mean_error[iGC,jGC]
    
    # Measurement-model mismatch: TROPOMI columns minus GEOS-Chem virtual TROPOMI columns 
    deltaY = obs_GC[:,0] - obs_GC[:,1] # [ppb]
    # If there are any nan's in the data, abort 
    if (np.any(np.isnan(deltaY)) or np.any(np.isnan(KK)) or np.any(np.isnan(obs_error))):
        print('missing values', fi)
        break
    
    # [****What are these steps?]
    # [****We're clearly building some terms for the inversion here]
    KK_t = KK.transpose() 
    KK_t2 = np.zeros(KK_t.shape, dtype=float)
    for k in range(KK_t.shape[1]):
        KK_t2[:,k] = KK_t[:,k]/obs_error[k]        

    # [****What are part1 and part2?]
    # [****What is happening here?]
    part1 = KK_t2@KK
    part2 = KK_t2@deltaY
   
    # Add part1 & part2 to sums 
    all_part1 += part1
    all_part2 += part2
        
# Prior covariance matrix
emis_error = np.zeros(n_clust);
emis_error.fill(0.5**2)          # 50% error, 1-sigma
inv_Sa = np.diag(1/emis_error)   # Inverse of prior covariance matrix

# Solve for posterior scaling factors
# [****Is that what "ratio" is?]
ratio = np.linalg.inv(gamma*all_part1 + inv_Sa)@(gamma*all_part2)
xhat = 1 + ratio # xhat = x_A + ratio = 1 + ratio

# Print some statistics
print('Min:',xhat.min(),'Mean:',xhat.mean(),'Max',xhat.max())

# Dictionary to store results
# [****But we don't save it out? Delete these lines?]
#met = {}
#met['all_part1'] = all_part1
#met['all_part2'] = all_part2
#met['ratio'] = ratio

# Save results
outputname = '/n/holyscratch01/jacob_lab/lshen/CH4/GEOS-Chem/Flexgrid/CPU_mexico_inversion/Step4_inversion_pixel/inversion_result.nc'
dataset = Dataset(outputname, 'w', format='NETCDF4_CLASSIC')
nvar = dataset.createDimension('nvar', n_clust)
nc_all_part1 = dataset.createVariable('all_part1', np.float32,('nvar','nvar'))
nc_all_part2 = dataset.createVariable('all_part2', np.float32,('nvar'))
nc_ratio = dataset.createVariable('ratio', np.float32,('nvar'))
nc_xhat = dataset.createVariable('xhat', np.float32, ('nvar'))
nc_all_part1[:,:] = all_part1    # [****What is this?]
nc_all_part2[:] = all_part2      # [****What is this?]
nc_ratio[:] = ratio              # [****What is this?]
nc_xhat[:] = xhat
dataset.close()
