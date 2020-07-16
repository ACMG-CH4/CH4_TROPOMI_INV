#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import numpy as np
from netCDF4 import Dataset
import xarray as xr
import pickle
import os

#----- define function -------
def save_obj(obj, name ):
    with open(name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open( name, 'rb') as f:
        return pickle.load(f)

def nearest_loc(loc0,table,tolerance=1):
    temp=np.abs(table-loc0)
    ind=temp.argmin()
    if temp[ind]>=tolerance:
        return np.nan
    else:
        return ind

#---- read obs error data ---
filename="/n/holyscratch01/jacob_lab/lshen/CH4/GEOS-Chem/Flexgrid/CPU_mexico_inversion/Step4_inversion_pixel/mean_error.nc"
data=xr.open_dataset(filename)
lon_GC=data['lon'].values
lat_GC=data['lat'].values
mean_error=data['error'].values
mean_error=mean_error**2
data.close()
mean_error=np.einsum('ij->ji', mean_error)

#-----------------------------
os.chdir("/n/holyscratch01/jacob_lab/lshen/CH4/GEOS-Chem/Flexgrid/CPU_mexico_inversion/data_converted/")
files=glob.glob("*.pkl")
files.sort()

all_part1=np.zeros([1207,1207],dtype=float)
all_part2=np.zeros([1207],dtype=float)
xlim=[-101.875,-86.875];ylim=[13,30]
xlim=[-107.8125,-80.9375];ylim=[10,36]
xlim=[-106.5625,-82.1875];ylim=[11,35]

for ifile in range(len(files)):
    print(files[ifile])
    met=load_obj(files[ifile])
    if met['obs_GC'].shape[0]==0:
        continue
    
    obs_GC=met['obs_GC']
    ind=np.where( (obs_GC[:,2]>=xlim[0]) & (obs_GC[:,2]<=xlim[1]) & (obs_GC[:,3]>=ylim[0]) & (obs_GC[:,3]<=ylim[1]) )
    if (len(ind[0])==0):
        continue

    obs_GC=obs_GC[ind[0],:]        
    KK=met['KK'][ind[0],:]
    NN=obs_GC.shape[0]

    #now lower the sensitivity to BC by 50%
    #KK[:,1199:]=KK[:,1199:]/2
    
    obs_error=np.zeros((NN,))    
    deltaY=obs_GC[:,0]-obs_GC[:,1]
    
    for iNN in range(NN):
        iGC=nearest_loc(obs_GC[iNN,2],lon_GC)
        jGC=nearest_loc(obs_GC[iNN,3],lat_GC)        
        obs_error[iNN]=mean_error[iGC,jGC]
            
    if (np.any(np.isnan(deltaY)) or np.any(np.isnan(KK)) or np.any(np.isnan(obs_error))):
        print('missing values', ifile, files[ifile])
        break
        
    KK_t=KK.transpose()
    
    #----------    
    KK_t2=np.zeros(KK_t.shape,dtype=float)
    for k in range(KK_t.shape[1]):
        KK_t2[:,k]=KK_t[:,k]/obs_error[k]        

    part1=KK_t2@KK
    part2=KK_t2@deltaY
    
    all_part1=all_part1+part1
    all_part2=all_part2+part2
        
        
emis_error=np.zeros(1207);emis_error.fill(0.5**2)
inv_Sa=np.diag(1/emis_error)
ratio=np.linalg.inv(all_part1+inv_Sa)@all_part2
print(ratio.max(),ratio.mean(),ratio.min())

met={}
met['all_part1']=all_part1
met['all_part2']=all_part2
met['ratio']=ratio

#======= save results ======
outputname='/n/holyscratch01/jacob_lab/lshen/CH4/GEOS-Chem/Flexgrid/CPU_mexico_inversion/Step4_inversion_pixel/inversion_result.nc'
dataset = Dataset(outputname,'w',format='NETCDF4_CLASSIC')
nvar=dataset.createDimension('nvar', 1207)

nc_all_part1 = dataset.createVariable('all_part1', np.float32,('nvar','nvar'))
nc_all_part2 = dataset.createVariable('all_part2', np.float32,('nvar'))
nc_ratio = dataset.createVariable('ratio', np.float32,('nvar'))

nc_all_part1[:,:]=all_part1
nc_all_part2[:]=all_part2
nc_ratio[:]=ratio

dataset.close()
