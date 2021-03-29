import numpy as np
import xarray as xr
from sklearn.cluster import KMeans

def make_cluster_file(land_cover_pth, save_pth, lat_min, lat_max, lon_min, lon_max, buffer_deg=5, land_threshold=0.25, n_clust=8):
    '''
    Generates the cluster file for an analytical inversion.
    
    Arguments
        land_cover_pth [str]   : Path to land cover file
        save_pth       [str]   : Where to save the cluster file
        lat_min        [float] : Minimum latitude 
        lat_max        [float] : Maximum latitude 
        lon_min        [float] : Minimum longitude 
        lon_max        [float] : Maximum longitude
        buffer_deg     [float] : Width of k-means buffer area in degrees
        land_threshold [float] : Minimum land fraction to include pixel as cluster
        nclust         [int]   : Number of clusters for k-means

    Returns
        ds_clusters    []      : xarray dataset containing clusters field formatted for HEMCO

    Notes
        - Land cover file default can be something like 'GEOSFP.20200101.CN.025x03125.NA.nc'
    '''

    # Load land cover data
    lc = xr.open_dataset(land_cover_pth)

    # Group together
    lc = (lc['FRLAKE'] + lc['FRLAND'] + lc['FRLANDIC']).drop('time').squeeze()
    
    # Subset the area of interest
    lc = lc.isel(lon=lc.lon>=lon_min, lat=lc.lat>=lat_min)
    lc = lc.isel(lon=lc.lon<=lon_max, lat=lc.lat<=lat_max)

    # Replace all values with NaN (to be filled later)
    clust = lc.where(lc == -9999.)

    # Set pixels in buffer areas to 0
    clust[:, (clust.lon < lon_min+buffer_deg) | (clust.lon > lon_max-buffer_deg)] = 0
    clust[(clust.lat < lat_min+buffer_deg) | (clust.lat > lat_max-buffer_deg), :] = 0

    # Fill in the NaNs with cluster values
    clust.values[clust.isnull()] = np.arange(1, clust.isnull().sum()+1)[::-1]

    # Assign buffer pixels to clusters
    # -----------------------------------------------------------------------------
    buffer_area = np.abs((clust.values > 0) - 1)
    
    # Get image coordinates of all pixels in buffer area
    irows = np.arange(buffer_area.shape[0])
    icols = np.arange(buffer_area.shape[1])
    irows = np.transpose(np.tile(irows,(len(icols),1)))
    icols = np.tile(icols,(len(irows),1)) * (buffer_area > 0)
    irows_good = irows[buffer_area > 0]
    icols_good = icols[buffer_area > 0]
    coords = [[icols_good[j], irows_good[j]] for j in range(len(irows_good))]
    
    # K-means
    X = np.array(coords)
    kmeans = KMeans(n_clusters=n_clust, random_state=0).fit(X)
    
    # Assign pixels to clusters
    highres_cluster_max = np.max(clust.values)
    n_rows = clust.shape[0]
    n_cols = clust.shape[1]
    for r in range(n_rows):
        for c in range(n_cols):
            if clust[r,c].values == 0:
                clust[r,c] = kmeans.predict([[c,r]])[0] + 1 + highres_cluster_max
    # -----------------------------------------------------------------------------

    if land_threshold:
        # Where there is no land, replace with NaN
        clust = clust.where(lc > land_threshold)

    # Make dataset
    da_clusters = clust.copy()
    ds_clusters = da_clusters.to_dataset(name='Clusters')

    # Add attribute metadata
    ds_clusters.lat.attrs['units'] = 'degrees_north'
    ds_clusters.lat.attrs['long_name'] = 'Latitude'
    ds_clusters.lon.attrs['units'] = 'degrees_east'
    ds_clusters.lon.attrs['long_name'] = 'Longitude'
    ds_clusters.Clusters.attrs['units'] = 'none'

    # Save
    if save_pth:
        ds_clusters.to_netcdf(save_pth)
    
    return ds_clusters