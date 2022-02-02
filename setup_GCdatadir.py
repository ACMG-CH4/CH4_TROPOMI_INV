import xarray as xr
import datetime
from joblib import Parallel, delayed


def zero_pad_num_hour(n):
    nstr = str(n)
    if len(nstr) == 1:
        nstr = '0'+nstr
    return nstr


def setup_GCdatadir(startday, endday, GC_source_pth, GC_destination_pth):
    '''
    This script sets up a directory containing hourly GEOS-Chem output diagnostics
    files. The hourly files are convenient for computing virtual TROPOMI columns 
    from the GEOS-Chem simulated atmosphere (to compare with the real TROPOMI columns).

    Arguments
        startday           [str] : First day of inversion period; formatted YYYYMMDD
        endday             [str] : Last day of inversion period; formatted YYYYMMDD
        GC_source_pth      [str] : GEOS-Chem output directory
        GC_destination_pth [str] : Target GEOS-Chem data directory in inversion workspace

    '''


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

    # For each day:
    def process(d):        
        # Load the SpeciesConc and LevelEdgeDiags data
        SpeciesConc_data = xr.load_dataset(f'{GC_source_pth}/GEOSChem.SpeciesConc.{d}_0000z.nc4')
        LevelEdgeDiags_data = xr.load_dataset(f'{GC_source_pth}/GEOSChem.LevelEdgeDiags.{d}_0000z.nc4')
        
        # For each hour:
        for h in hours:
            
            # Select data for that hour
            SpeciesConc_for_hour = SpeciesConc_data.isel(time=slice(h,h+1,1))
            LevelEdgeDiags_for_hour = LevelEdgeDiags_data.isel(time=slice(h,h+1,1))
            
            # Save to new .nc4 file at destination
            SpeciesConc_save_pth = f'{GC_destination_pth}/GEOSChem.SpeciesConc.{d}_{zero_pad_num_hour(h)}00z.nc4'
            LevelEdgeDiags_save_pth = f'{GC_destination_pth}/GEOSChem.LevelEdgeDiags.{d}_{zero_pad_num_hour(h)}00z.nc4'
            SpeciesConc_for_hour.to_netcdf(SpeciesConc_save_pth)
            LevelEdgeDiags_for_hour.to_netcdf(LevelEdgeDiags_save_pth)

    results = Parallel(n_jobs=-1) (delayed(process)(day) for day in days)    
    print(f'Set up hourly data files in {GC_destination_pth}')


if __name__ == '__main__':
    import sys

    startday = sys.argv[1]
    endday = sys.argv[2]
    GC_source_pth = sys.argv[3]
    GC_destination_pth = sys.argv[4]

    setup_GCdatadir(startday, endday, GC_source_pth, GC_destination_pth)
