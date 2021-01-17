import datetime
import xarray as xr
import os

def fill_missing_hour(run_name, run_dirs_pth, prev_run_pth, start_day):
    '''
    The purpose of this script is to address the problem that output files for
    the first day of a GEOS-Chem simulation do not include data for the first 
    hour of the day; they go from 1-23h, instead of 0-23h. The solution is to
    use the final output file from each posterior simulation, which contains
    data for hour 0 of the last day, and then combine that file with the first
    output file of the most recent simulation. This needs to be done for both 
    SpeciesConc and LevelEdgeDiags, and it needs to be done for every run
    directory: i.e., for every perturbed cluster simulation, and also for every
    posterior simulation.

    Example:
    A GEOS-Chem perturbation simulations for week 1, from 20180501 to 20180508. 
    The SpeciesConc and LevelEdgeDiags files (i.e., the output files) for the
    first day, 20180501, are missing data for hour 0. We merge data from the
    final output files of my spinup simulation, which contain data for only hour
    0 of 20180501, into the latest output files for 20180501. Now the data for
    week 1 are complete, so the inversion can be run for week 1, and run the 
    posterior simulation for week 1. The final SpeciesConc and LevelEdgeDiags
    files from the posterior simulation are for 20180508, with data only for
    hour 0. We can then run the perturbation simulations for week 2, from
    20180508 to 20180515. Now the output files for 20180508 are missing data
    for hour 0. We get this data from the 20180508 output file from the last
    posterior simulation. 

    Arguments
        run_name         [str] : Name for this set of runs
        run_dirs_pth     [str] : Path to the parent folder of the GEOS-Chem run
                                 directory or directories where we want to fill
                                 missing data
        prev_run_pth     [str] : Path to the spinup or posterior run directory
        start_day        [str] : First day of simulation, for which the daily
                                 output is missing data; e.g., "20180501"
    '''

    # List run directories
    contents = os.listdir(run_dirs_pth)
    rundirs = [r for r in contents if run_name in r]

    # Process them
    for r in rundirs:
        # Load hour zero from end of spinup run or previous posterior simulation
        prev_file_SC = f'{prev_run_pth}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4'
        prev_file_LE = f'{prev_run_pth}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4'
        prev_data_SC = xr.load_dataset(prev_file_SC)
        prev_data_LE = xr.load_dataset(prev_file_LE)
        
        # Load output SpeciesConc and LevelEdgeDiags file
        output_file_SC = f'{run_dirs_pth}/{r}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4'
        output_data_SC = xr.load_dataset(output_file_SC)
        output_file_LE = f'{run_dirs_pth}/{r}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4'
        output_data_LE = xr.load_dataset(output_file_LE)
        
        # Merge output and copied datasets
        merged_data_SC = xr.merge([output_data_SC, prev_data_SC])
        merged_data_LE = xr.merge([output_data_LE, prev_data_LE])
        
        # Replace original files that were missing the first hour 
        merged_data_SC.to_netcdf(output_file_SC)
        merged_data_LE.to_netcdf(output_file_LE)


if __name__ == '__main__':
    import sys

    run_name = sys.argv[1]
    run_dirs_pth = sys.argv[2]
    prev_run_pth = sys.argv[3]
    start_day = sys.argv[4]

    fill_missing_hour(run_name, run_dirs_pth, prev_run_pth, start_day)
