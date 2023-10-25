"""
 Plot sea ice thickness fields on polar stereographic projection
"""
# To clean all variables similar to "clear all" in Matlab
# use: %run script
# Clean memory - like start a new shell, this also cleans the memory when
# script is submitted
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.colors import ListedColormap
import importlib
import mod_colormaps
#import matplotlib.gridspec as gridspec
import sys

sys.path.append('/home/ddmitry/codes/MyPython/hycom_utils')
sys.path.append('/home/ddmitry/codes/MyPython/draw_map')
sys.path.append('/home/ddmitry/codes/MyPython')


plt.close('all')

rhl = 'EVP'  # EVP or EAP
yr   = 2017
mo   = 3
mday = 15

f_hycomgrid = 0    # =1 plot field on HYCOM grid as well, for comparison

if rhl == 'EVP':
	pthnc = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_cice/'
else:
	pthnc = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_EAPcice/'

pthtopo  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/'

print("Plotting hice: {0}/{1}/{2}".format(mday,mo,yr))

def plot_hycom_grid(hice,amin,amax,flnm,btx):
	fig = plt.figure(figsize=(8,8))
	# plt.clf()
	#ax = fig.add_axes([0.1,0.1,0.8,0.8])
	#plt.imshow(hice, cmap='viridis')
	plt.imshow(hice, cmap=clrmp, 
						 norm=colors.Normalize(vmin=amin, vmax=amax))
	plt.gca().invert_yaxis()
	plt.clim(amin, amax)
	plt.xticks([])
	plt.yticks([])

	# Colorbar
	#plt.colorbar(extend='min')
	plt.colorbar(aspect=40)

	plt.title('0.04 HYCOM-CICE, hice, '+flnm)

	bottom_text(btx)
	plt.show(block=False)


from mod_utils_fig import bottom_text
import mod_read_hycom
importlib.reload(mod_read_hycom)
from mod_read_hycom import readnc_grid_topo, read_hycom
flgrid = 'depth_ARCc0.04_17DD.nc'
lonc, latc, HH = readnc_grid_topo(pthtopo,flgrid,f_lon=360)


flnm = '022_cice.{0}-{1:02d}-{2:02d}.nc'.format(yr,mo,mday)

nc = NetCDFFile(pthnc+flnm)
hice = nc.variables['hi'][:].squeeze()
#hice.shape
#hice = np.squeeze(hice)
print("hice shape=",hice.shape)


# Set up a colormap:
clrmp = copy(plt.cm.cividis)
clrmp.set_under(color=[0.9, 0.9, 0.9])
clrmp.set_bad(color=[0.5, 0.5, 0.5])

# create figure and axes instances
plt.ion()  # enables interactive mode
btx = 'plot_hice_polar.py'
amin = 0.15
amax = 1.0

# Plot on a HYCOM grid
if f_hycomgrid == 1:
  plot_hycom_grid(hice,amin,amax,flnm,btx)

#
# Close figure:
# plt.close(fig)
# plt.close('all')


#
Nval=4
importlib.reload(mod_colormaps)
#clrmpBR = mod_colormaps.clrmp_Nvalues(Nval,256)
clrmpBR = plt.cm.gist_rainbow_r
clrmpBR.set_under(color=[0.9, 0.9, 0.9])
clrmpBR.set_bad(color=[0.5, 0.5, 0.5])


fig2 = plt.figure(figsize=(8,8), constrained_layout=False)
ax1 = plt.axes([0.1, 0.1, 0.7, 0.7])

plt.sca(ax1)

m = Basemap(projection='npstere',boundinglat=50,lon_0=0,resolution='l')
# draw coastlines, state and country boundaries, edge of map.
#m.drawcoastlines(color=[0.5,0.5,0.5])

# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,360.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#m.fillcontinents(color=[0.8, 0.8, 0.8])

#hice[hice<=0.01]=np.nan
ny = hice.shape[0]; nx = hice.shape[1]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly spaced grid.
x, y = m(lons, lats) # compute map proj coordinates.
xh, yh = m(lonc,latc)  # hycom coordinates on the projections
im1 = plt.pcolormesh(xh,yh,hice,shading='flat',cmap=clrmpBR,
               norm=colors.Normalize(vmin=0.2, vmax=Nval))

plt.title('0.04 HYCOM-CICE, hice, {0}, {1}'.format(rhl,flnm))

xl1 = 1500.e3
xl2 = 7000.e3
yl1 = 2000.e3
yl2 = 7500.e3

ax1.axis('equal')
ax1.set_xlim([xl1, xl2])
ax1.set_ylim([yl1, yl2])


bottom_text(btx)


# Colorba
# Create axes for colorbar
ax2 = fig2.add_axes([ax1.get_position().x1+0.02, ax1.get_position().y0,0.02,
                    ax1.get_position().height])
plt.colorbar(im1, cax=ax2)

plt.show()











