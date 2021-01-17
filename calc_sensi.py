import numpy as np
import xarray as xr
import datetime

def zero_pad_num(n):
    nstr = str(n)
    if len(nstr) == 1:
        nstr = '000'+nstr
    if len(nstr) == 2:
        nstr = '00'+nstr
    if len(nstr) == 3:
        nstr = '0'+nstr
    return nstr

def zero_pad_num_hour(n):
    nstr = str(n)
    if len(nstr) == 1:
        nstr = '0'+nstr
    return nstr

def calc_sensi(startday, endday, run_dirs_pth, run_name, sensi_save_pth):
    '''
    nclust = count the number of clusters
    nlon = count the number of longitudes
    nlat = count the number of latitudes
    nlev = count the number of vertical levels
    
    for each day:
        load the base run SpeciesConc file
        for each hour:
            base = extract the base run data for the hour
            Sensi = np.empty((nlon, nlat, nlev, nclust))
            Sensi.fill(np.nan)
            for each cluster:
                load the SpeciesConc .nc file for the cluster and day
                pert = extract the data for the hour
                sens = pert - base
                Sensi[:,:,:,cluster] = sens
            save Sensi as netcdf with appropriate coordinate variables

    e.g. Lu's Sensi files look like:
    
    <xarray.Dataset>
    Dimensions:  (grid: 1207, lat: 105, lev: 47, lon: 87)
    Coordinates:
      * lon      (lon) float64 -107.8 -107.5 -107.2 -106.9 ... -81.56 -81.25 -80.94
      * lat      (lat) float64 10.0 10.25 10.5 10.75 11.0 ... 35.25 35.5 35.75 36.0
      * lev      (lev) int32 1 2 3 4 5 6 7 8 9 10 ... 38 39 40 41 42 43 44 45 46 47
      * grid     (grid) int32 1 2 3 4 5 6 7 8 ... 1201 1202 1203 1204 1205 1206 1207
    Data variables:
        Sensi    (grid, lev, lat, lon) float32 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0
    
    '''

    nlon = 52
    nlat = 61
    nlev = 47
    nclust = 235+8

    # Make date range
    days = []
    dt = datetime.datetime.strptime(startday, '%Y%m%d')
    dt_max = datetime.datetime.strptime(endday, '%Y%m%d')
    while dt < dt_max:
        dt_str = str(dt)[0:10].replace('-','')
        days.append(dt_str)
        delta = datetime.timedelta(days=1)
        dt += delta

    hours = range(24)
    clust = range(nclust)

    for d in days:
        base_data = xr.load_dataset(run_dirs_pth+'/'+run_name+'_0000/OutputDir/GEOSChem.SpeciesConc.'+d+'_0000z.nc4')
        for h in hours:
            base = base_data['SpeciesConc_CH4'][h,:,:,:]
            Sensi = np.empty((nclust, nlev, nlat, nlon))
            Sensi.fill(np.nan)
            for c in clust:
                cstr = zero_pad_num(c+1) # Because clusters are numbered 1..nclust
                pert_data = xr.load_dataset(run_dirs_pth+'/'+run_name+'_'+cstr+'/OutputDir/GEOSChem.SpeciesConc.'+d+'_0000z.nc4')
                pert = pert_data['SpeciesConc_CH4'][h,:,:,:]
                sens = pert.values - base.values
                Sensi[c,:,:,:] = sens
            Sensi = xr.DataArray(Sensi, 
                                 coords=(np.arange(1,nclust+1),np.arange(1,nlev+1), base.lat, base.lon), 
                                 dims=['clust','lev','lat','lon'],
                                 name='Sensitivities')
            Sensi = Sensi.to_dataset()
            Sensi.to_netcdf(sensi_save_pth+'/Sensi_'+d+'_'+zero_pad_num_hour(h)+'.nc')

    print("Saved GEOS-Chem sensitivity files to {}".format(sensi_save_pth))

if __name__ == '__main__':
    import sys

    startday = sys.argv[1]
    endday = sys.argv[2]
    run_dirs_pth = sys.argv[3]
    run_name = sys.argv[4]
    sensi_save_pth = sys.argv[5]

    calc_sensi(startday, endday, run_dirs_pth, run_name, sensi_save_pth)
