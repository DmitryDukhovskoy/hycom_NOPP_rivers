"""
	Elastic-anisotropic-plastic rheology (kdyn=2) uses extra prognositc variable
	In contrast to the isotropic EVP rheology, the anisotropic plastic 
	yield curve within the EAP rheology depends on the relative orientation 
	of the diamond shaped floes (unit vector r in Figure 5), with respect to 
	the principal direction of the deformation rate. 
	Local anisotropy of the sea ice cover is accounted for by an additional 
	prognostic variable, the structure tensor A 
	Add required fields: A11, A22, A12 to the restart file

	Start with copying an EVP restart file to an EAP restart file, then
	add missing fields
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
from netCDF4 import Dataset as ncFile
import numpy as np
from copy import copy
import importlib

drnm = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/restart/'
rname_out = 'cice.restart_117d-eap.nc'

print('Appending to ',rname_out)
#ncin = drnm+rname_in+'.nc'

try: ncout.close()
except: pass

ncout =ncFile(drnm+rname_out, mode='a', format='NETCDF4')
dimi = 'ni'
dimj = 'nj'

dmm = ncout.variables['stress12_1'][:]
nj, ni = dmm.shape

eps0 = 1.e-10
a05 = np.zeros((nj,ni))+0.5
a0  = np.zeros((nj,ni))+eps0

print('Adding a11_1')
A11_1 = ncout.createVariable('a11_1',np.float64,('nj','ni'))
A11_1[:,:] = a05

print('Adding a11_2')
A11_2 = ncout.createVariable('a11_2',np.float64,('nj','ni'))
A11_2[:,:] = a05

print('Adding a11_3')
A11_3 = ncout.createVariable('a11_3',np.float64,('nj','ni'))
A11_3[:,:] = a05

print('Adding a11_4')
A11_4 = ncout.createVariable('a11_4',np.float64,('nj','ni'))
A11_4[:,:] = a05 

print('Adding a12_1')
A12_1 = ncout.createVariable('a12_1',np.float64,('nj','ni'))
A12_1[:,:] = a0

print('Adding a12_2')
A12_2 = ncout.createVariable('a12_2',np.float64,('nj','ni'))
A12_2[:,:] = a0 

print('Adding a12_3')
A13_3 = ncout.createVariable('a12_3',np.float64,('nj','ni'))
A13_3[:,:] = a0 

print('Adding a12_4')
A14_4 = ncout.createVariable('a12_4',np.float64,('nj','ni'))
A14_4[:,:] = a0 

ncout.close()
print('Finished, new file: ',drnm+rname_out)


