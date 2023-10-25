# Create polar coordinate 
# for projecting HYCOM-CICE output
import numpy as np
import imp
import sub_module1
from sub_module1 import dist_sphcrd


# Define min latitude, x,y space steps:
phi_min = 50
dx = 40000   
dy = dx

#
#imp.reload(sub_module1)
#from sub_module1 import dist_sphcrd
# Determine # of grid points
dist = dist_sphcrd(phi_min,270,90,270) # half/domain width
ni0 = np.ceil(dist/dx)
# Central point - note index1 = 0
# N. Pole
i0 = ni0  
j0 = ni0
NI = np.int(2*ni0+1)
NJ = np.int(NI)
print('Grid: {0} x {1} = {2} '.format(NI,NJ,NI*NJ)) 

PHI = np.zeros([NJ,NI])
PHI[:,:] = -999
LMB = np.zeros([NJ,NI])  # longitudes
LMB[:,:] = -999
Rearth = 6371.0e3

for ii in range(NI-1):
  for jj in range(NJ-1):
    dltX = (ii-np.int(i0))*dx
    dltY = (jj-np.int(j0))*dy
    R = np.sqrt(dltX*dltX+dltY*dltY)
    phi_rad = R/Rearth  # latitude
    lmb_rad = np.arctan2(dltY,dltX)
    PHI[jj,ii] = phi_rad*180./np.pi
    LMB[jj,ii] = lmb_rad*180./np.pi 
print('Min/Max Lat: {0}, {1}'.format(np.min(PHI[PHI>-999]),np.max(PHI)))
print('Min/Max Lon: {0}, {1}'.format(np.min(LMB[LMB>-999]),np.max(LMB)))

# run ./polar_grid.py
#from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
plt.imshow(PHI, cmap='inferno'); 
plt.colorbar()
plt.show()




