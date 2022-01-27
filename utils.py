import numpy as np
import xarray as xr
from functools import partial
import pyproj
from shapely.geometry.polygon import Polygon
import shapely.ops as ops
import cartopy
import cartopy.crs as ccrs

def calculate_gridcell_areas(state_vector, mask, dlat, dlon):
    '''
    Compute the surface areas of grid cells in the region of interest, in m2.
    '''

    xgrid = range(len(state_vector.lon.values))
    ygrid = range(len(state_vector.lat.values))
    areas = []
    for j in xgrid:
        for i in ygrid:
            if mask.values[i,j] == 1:
                lat_top = state_vector.lat.values[i] + dlat
                lat_bot = state_vector.lat.values[i] - dlat
                lon_left = state_vector.lon.values[j] - dlon
                lon_righ = state_vector.lon.values[j] + dlon
                geom = Polygon([(lon_left, lat_bot), (lon_left, lat_top), (lon_righ, lat_top), 
                                (lon_righ, lat_bot), (lon_left, lat_bot)])
                geom_area = ops.transform(
                                        partial(
                                            pyproj.transform,
                                            pyproj.Proj(init='EPSG:4326'),
                                            pyproj.Proj(
                                                proj='aea',
                                                lat_1=geom.bounds[1],
                                                lat_2=geom.bounds[3])),
                                        geom)
                areas.append(geom_area.area)
    return areas


def sum_total_emissions(emissions, areas, state_vector_labels, last_ROI_element):
    '''
    Function to sum total emissions across the region of interest.
    
    Arguments:
        emissions           : emissions dataarray
        areas               : list of pixel areas (in m2) for region of interest
        state_vector_labels : state vector element IDs, dataarray
        last_ROI_element    : ID of last state vector element in the region of interest
        
    Returns:
        Total emissions in Tg/y
    '''

    s_per_d = 86400
    d_per_y = 365
    tg_per_kg = 1e-9
    xgrid = range(len(state_vector_labels.lon.values))
    ygrid = range(len(state_vector_labels.lat.values))    
    mask = (state_vector_labels <= last_ROI_element)
    emiss = []
    for j in xgrid:
        for i in ygrid:
            if mask.values[i,j] == 1:
                emiss.append(emissions.values[i,j])
    total = np.sum(np.asarray([areas[r] * emiss[r] for r in range(len(areas))]))
    return total * s_per_d * d_per_y * tg_per_kg


def read_tropomi(filename):
    """
    Read TROPOMI data and save important variables to dictionary.

    Arguments
        filename [str]  : TROPOMI netcdf data file to read

    Returns
        met      [dict] : Dictionary of important variables from TROPOMI:
                            - CH4
                            - Latitude
                            - Longitude
                            - QA value
                            - UTC time
                            - Time (utc time reshaped for orbit)
                            - Averaging kernel
                            - SWIR albedo
                            - NIR albedo
                            - Blended albedo
                            - CH4 prior profile
                            - Dry air subcolumns
                            - Latitude bounds
                            - Longitude bounds
                            - Vertical pressure profile
    """

    # Initialize dictionary for TROPOMI data
    met = {}
    
    # Store methane, QA, lat, lon
    data = xr.open_dataset(filename, group='PRODUCT')
    met['methane'] = data['methane_mixing_ratio_bias_corrected'].values[0,:,:]
    met['qa_value'] = data['qa_value'].values[0,:,:]
    met['longitude'] = data['longitude'].values[0,:,:]
    met['latitude'] = data['latitude'].values[0,:,:]
    
    # Store UTC time
    referencetime = data['time'].values
    delta_time = data['delta_time'][0].values
    strdate = []
    if delta_time.dtype == '<m8[ns]':
        strdate = referencetime+delta_time
    elif delta_time.dtype == '<M8[ns]':
        strdate = delta_time
    else:
        print(delta_time.dtype)
        pass
    met['utctime'] = strdate
    
    # Store time for whole orbit
    times = np.zeros(shape=met['longitude'].shape, dtype='datetime64[ns]')
    for kk in range(met['longitude'].shape[0]):
        times[kk,:] = strdate[kk]
    met['time'] = times
    data.close()
 
    # Store column averaging kernel and SWIR, NIR surface albedo
    data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/DETAILED_RESULTS')
    met['column_AK'] = data['column_averaging_kernel'].values[0,:,:,::-1]
    met['swir_albedo'] = data['surface_albedo_SWIR'].values[0,:,:]
    met['nir_albedo'] = data['surface_albedo_NIR'].values[0,:,:]
    met['blended_albedo'] = 2.4*met['nir_albedo'] - 1.13*met['swir_albedo']
    data.close()
    
    # Store methane prior profile, dry air subcolumns
    data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
    met['methane_profile_apriori'] = data['methane_profile_apriori'].values[0,:,:,::-1] # mol m-2
    met['dry_air_subcolumns'] = data['dry_air_subcolumns'].values[0,:,:,::-1]           # mol m-2
    
    # Also get pressure interval and surface pressure for use below
    pressure_interval = data['pressure_interval'].values[0,:,:]/100                     # Pa -> hPa
    surface_pressure = data['surface_pressure'].values[0,:,:]/100                       # Pa -> hPa
    data.close()

    # Store lat, lon bounds for pixels
    data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/GEOLOCATIONS')
    met['longitude_bounds'] = data['longitude_bounds'].values[0,:,:,:]
    met['latitude_bounds'] = data['latitude_bounds'].values[0,:,:,:]
    data.close()
    
    # Store vertical pressure profile
    N1 = met['methane'].shape[0]
    N2 = met['methane'].shape[1]
    pressures = np.zeros([N1,N2,13], dtype=np.float32)
    pressures.fill(np.nan)
    for i in range(12+1):
        pressures[:,:,i] = surface_pressure - i*pressure_interval
    met['pressures'] = pressures
    
    return met


def count_obs_in_mask(mask, df):
    """
    Count the number of observations in a boolean mask
    mask is boolean xarray data array
    df is pandas dataframe with lat, lon, etc.
    """
    
    reference_lat_grid = mask['lat'].values
    reference_lon_grid = mask['lon'].values
    
    # Query lats/lons
    query_lats = df['lat'].values
    query_lons = df['lon'].values
    
    # Loop
    bad_ind = []
    for k in range(len(df)):
        # Find closest reference coordinates to selected lat/lon bounds
        ref_lat_ind = np.abs(reference_lat_grid - query_lats[k]).argmin()
        ref_lon_ind = np.abs(reference_lon_grid - query_lons[k]).argmin()
        # If not in mask, save as bad index
        if mask[ref_lat_ind,ref_lon_ind] == 0:
            bad_ind.append(k)
    
    # Drop bad indexes and count remaining entries
    df_copy = df.copy()
    df_copy = df_copy.drop(df_copy.index[bad_ind])
    n_obs = len(df_copy)
    
    return n_obs


def plot_field(ax, field, cmap, plot_type='pcolormesh', lon_bounds=None, lat_bounds=None, 
               levels=21, vmin=None, vmax=None, title=None, cbar_label=None, mask=None, only_ROI=False,
               state_vector_labels=None, last_ROI_element=None):
    '''
    Function to plot inversion results.
    
    Arguments
        ax         : matplotlib axis object
        field      : xarray dataarray
        cmap       : colormap to use, e.g. 'viridis'
        plot_type  : 'pcolormesh' or 'imshow'
        lon_bounds : [lon_min, lon_max]
        lat_bounds : [lat_min, lat_max]
        levels     : number of levels for pcolormesh option
        vmin       : colorbar lower bound
        vmax       : colorbar upper bound
        title      : plot title
        cbar_label : colorbar label
        mask       : mask for region of interest, boolean dataarray
        only_ROI   : zero out data outside the region of interest, true or false
    '''
    
    # Select map features
    oceans_50m = cartopy.feature.NaturalEarthFeature('physical', 'ocean', '50m')
    lakes_50m = cartopy.feature.NaturalEarthFeature('physical', 'lakes', '50m')
    states_provinces_50m = cartopy.feature.NaturalEarthFeature('cultural','admin_1_states_provinces_lines', '50m')
    ax.add_feature(cartopy.feature.BORDERS, facecolor='none')
    ax.add_feature(oceans_50m, facecolor=[1,1,1], edgecolor='black')
    ax.add_feature(lakes_50m, facecolor=[1,1,1], edgecolor='black')
    ax.add_feature(states_provinces_50m, facecolor='none', edgecolor='black')
    
    # Show only ROI values?
    if only_ROI:
        field = field.where((state_vector_labels <= last_ROI_element))  
    
    # Plot
    if plot_type == 'pcolormesh':
        field.plot.pcolormesh(cmap=cmap, levels=levels, ax=ax,
                              vmin=vmin, vmax=vmax, cbar_kwargs={'label':cbar_label,
                                                                 'fraction':0.041, 
                                                                 'pad':0.04});
    elif plot_type == 'imshow':
        field.plot.imshow(cmap=cmap, ax=ax,
                          vmin=vmin, vmax=vmax, cbar_kwargs={'label':cbar_label,
                                                             'fraction':0.041, 
                                                             'pad':0.04});
    else:
        raise ValueError('plot_type must be "pcolormesh" or "imshow"')
    
    # Zoom on ROI?
    if lon_bounds and lat_bounds:
        extent = [lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]]
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    # Show boundary of ROI?
    if mask is not None:
        mask.plot.contour(levels=1,colors='k',linewidths=2,ax=ax)
    
    # Remove duplicated axis labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0)
    gl.right_labels = False
    gl.top_labels = False
    
    # Title
    if title:
        ax.set_title(title)