"""
 Plot sea ice fields on polar stereographic projection
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
import matplotlib.gridspec as gridspec
import sys

sys.path.append('/home/ddmitry/codes/MyPython/hycom_utils')
sys.path.append('/home/ddmitry/codes/MyPython/draw_map')
sys.path.append('/home/ddmitry/codes/MyPython')


plt.close('all')

rhl = 'EAP'  # EVP or EAP
yr   = 2017
mo   = 8
mday = 31

f_hycomgrid = 0    # =1 plot field on HYCOM grid as well, for comparison

if rhl == 'EVP':
	pthnc = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_cice/'
#	pthnc = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_cice_old/'
else:
	pthnc = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_EAPcice/'

pthtopo  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/'

print("Plotting aice: {0}/{1}/{2}".format(mday,mo,yr))

#def get_grid():
#  pthgrid = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/'
#  flgrd = 'depth_ARCc0.04_17DD.nc'
#  nc = NetCDFFile(pthgrid+flgrd)
#  lonc = nc.variables['Longitude'][:]
#  latc = nc.variables['Latitude'][:]
#  lonc[lonc<0]=lonc[lonc<0]+360
#
#  return latc,lonc
  
def plot_hycom_grid(aice,amin,amax,flnm,btx):
	fig = plt.figure(figsize=(8,8))
	# plt.clf()
	#ax = fig.add_axes([0.1,0.1,0.8,0.8])
	#plt.imshow(aice, cmap='viridis')
	plt.imshow(aice, cmap=clrmp, 
						 norm=colors.Normalize(vmin=amin, vmax=amax))
	plt.gca().invert_yaxis()
	plt.clim(amin, amax)
	plt.xticks([])
	plt.yticks([])

	# Colorbar
	#plt.colorbar(extend='min')
	plt.colorbar(aspect=40)

	plt.title('0.04 HYCOM-CICE, aice, {0}, {1}'.format(rhl,flnm))

	bottom_text(btx)
	plt.show(block=False)


from mod_utils_fig import bottom_text
import mod_read_hycom
importlib.reload(mod_read_hycom)
from mod_read_hycom import readnc_grid_topo, read_hycom
flgrid = 'depth_ARCc0.04_17DD.nc'
lonc, latc, HH = readnc_grid_topo(pthtopo,flgrid,f_lon=360)
#lonc[lonc<0]=lonc[lonc<0]+360

#latc, lonc = get_grid()


#smo = str(mo).zfill(2)
#smday = str(mday).zfill(2)
#flnm  = '022_cice.{0}-{1}-{2}.nc'.format(yr,smo,smday)
flnm = '022_cice.{0}-{1:02d}-{2:02d}.nc'.format(yr,mo,mday)

nc = NetCDFFile(pthnc+flnm)
aice = nc.variables['aice'][:].squeeze()
#aice.shape
#aice = np.squeeze(aice)
print("aice shape=",aice.shape)


# Set up a colormap:
clrmp = copy(plt.cm.viridis)
clrmp.set_under(color=[0.9, 0.9, 0.9])
clrmp.set_bad(color=[0.5, 0.5, 0.5])

# create figure and axes instances
plt.ion()  # enables interactive mode
btx = 'plot_cice_polar.py'
amin = 0.15
amax = 1.0


if f_hycomgrid == 1:
  plot_hycom_grid(aice,amin,amax,flnm,btx)

#
# Close figure:
# plt.close(fig)
# plt.close('all')


#
#importlib.reload(mod_colormaps)
#clrmpBR = create_BlRd_clrmp(256) 
#clrmpBR = mod_colormaps.clrmp_BlRd(256)
clrmpBR = mod_colormaps.clrmp_BlGrRd(256)
clrmpBR.set_under(color=[0.9, 0.9, 0.9])
clrmpBR.set_bad(color=[0.5, 0.5, 0.5])

cmap1 = plt.cm.plasma
cmap1.set_under(color=[0.9, 0.9, 0.9])
cmap1.set_bad(color=[0.5, 0.5, 0.5])

fig2 = plt.figure(figsize=(8,8), constrained_layout=False)
#gs1 = fig2.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.9,
#                         wspace=0.05)
ax1 = plt.axes([0.1, 0.1, 0.7, 0.7])
#ax2 = plt.axes([0.85, 0.1, 0.02, 0.7])

plt.sca(ax1)

#latc = nc.variables['TLAT'][:]
#lonc = nc.variables['TLON'][:] - 1st column wrong longitudes
# if draw costline = pick resolution='h', or 'l'
# if no coastline resolution='none'
m = Basemap(projection='npstere',boundinglat=50,lon_0=0,resolution='l')
#m = Basemap(projection='npstere',llcrnrlon=-45.,
#            llcrnrlat=45., urcrnrlon=135., urcrnrlat=45.,
#            boundinglat=50.,
#            lon_0=0., lat_0=90.,rsphere=6371200.,
#            resolution='l')

# draw coastlines, state and country boundaries, edge of map.
#m.drawcoastlines(color=[0.5,0.5,0.5])

# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,360.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#m.fillcontinents(color=[0.8, 0.8, 0.8])

data = aice
ny = data.shape[0]; nx = data.shape[1]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly spaced grid.
x, y = m(lons, lats) # compute map proj coordinates.
xh, yh = m(lonc,latc)  # hycom coordinates on the projections
#im1 = plt.pcolormesh(xh,yh,aice,shading='flat',cmap=clrmpBR,
#               norm=colors.Normalize(vmin=0.15, vmax=1.0))
im1 = plt.pcolormesh(xh,yh,aice,shading='flat',cmap=cmap1,
               norm=colors.Normalize(vmin=0.15, vmax=1.0))
#plt.set_cmap('viridis_r')

plt.title('0.04 HYCOM-CICE, aice, {0}, {1}'.format(rhl,flnm))

#xl1 = 1500.
#xl2 = 8000.e3
#yl1 = 15.e3
#yl2 = 8000.e3

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











