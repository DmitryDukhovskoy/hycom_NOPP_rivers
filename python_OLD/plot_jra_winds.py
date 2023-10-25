"""
  Plot JRA-55 wind fields 
  Donwloaded from https://rda.ucar.edu/datasets/ds628.0/#!access
	6-hr surface analysis fields on model grid
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from copy import copy
#import matplotlib as mtlb
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.colors import ListedColormap
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy.interpolate import griddata
import importlib

#import sys
#sys.path.append('/usr/lib64/python3.8/site-packages/matplotlib/')


plt.close('all')

# assumed whole month data
mo  = 6
yr  = 2017
dS  = 1
dE  = 30
dth = 6       # saved time step
 
pthnc = '/nexsan/people/ddmitry/Net_ocean/JRA_55/'
uflnm = 'usurf_{0}{1:02d}{2:02}00_{0}{1:02d}{3:2d}18.nc'.format(yr,mo,dS,dE)
vflnm = 'vsurf_{0}{1:02d}{2:02}00_{0}{1:02d}{3:2d}18.nc'.format(yr,mo,dS,dE)

ncu = NetCDFFile(pthnc+uflnm)
ncv = NetCDFFile(pthnc+vflnm)

#aa = np.random.randn(7,4)
#aa[::-1] - flips array bottom to top
# Read grid
lonc = ncu.variables['g4_lon_2'][:].tolist()  # len(lonc) - size of the list
latc = ncu.variables['g4_lat_1'][::-1]  # flip: S to North

lonc.append(360.)
lonc = np.array(lonc)

nn = lonc.shape[0]
mm = latc.shape[0]

# Read U by hrs, flip S-N orientation and append extra column
# to make continuous data

ntime = np.int((dE-dS+1)*24/dth)
Uav = np.zeros((mm,nn))
Vav = np.zeros((mm,nn))

for ii in range(ntime):
  print("Reading record ",ii)
  uu = np.zeros((mm,nn))
  vv = np.zeros((mm,nn))
  uin = ncu.variables['UGRD_GDS4_HTGL'][ii,::-1,:].squeeze()
  vin = ncv.variables['VGRD_GDS4_HTGL'][ii,::-1,:].squeeze()
  uu[:,:-1] = uin
  uu[:,-1]  = uin[:,0]
  vv[:,:-1] = vin
  vv[:,-1]  = vin[:,0]

  Uav += uu
  Vav += vv

del uin, vin

Uav = Uav/np.float(ntime)
Vav = Vav/np.float(ntime)
Sav = np.sqrt(Uav**2+Vav**2)

# create figure and axes instances
plt.ion()  # enables interactive mode

# Quick check:
#plt.figure(1, figsize=(8,8))
#plt.pcolormesh(lonc,latc,Sav, cmap='viridis')
#plt.clim(0., 5.)

# Colorbar
#plt.colorbar(extend='min')
#plt.colorbar(aspect=40)

# Create 2D grid of lons/lats
lon2d_in, lat2d_in = np.meshgrid(lonc, latc)

# Create polar orthogonal grid
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
ny = Sav.shape[0]; nx = Sav.shape[1]
#lon_new, lat_new = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly spaced grid.
#xnew, ynew = m(lon_new, lat_new) # compute map proj coordinates.

# generate a grid that is equally spaced in a plot with the current pojection
lon_new, lat_new, xnew, ynew = m.makegrid(100,100, returnxy=True)
xh, yh = m(lon2d_in,lat2d_in)  # jra coordinates on the projections
xnew = xnew.astype('float64')
ynew = ynew.astype('float64')



# ============================================
# Project the data onto the new grid
points = np.column_stack((xh.flatten(),yh.flatten()))
unew = griddata(points, Uav.flatten(), (xnew, ynew), method='linear')
vnew = griddata(points, Vav.flatten(), (xnew, ynew), method='linear')
 
#
# Rotate U, V vectors to match N and E on the map projection
import mod_utils_map
#importlib.reload(mod_utils_map)
from mod_utils_map import rotate_vector_cart2polar
lonPol = lon_new
latPol = lat_new
Ur, Vr = rotate_vector_cart2polar(unew,vnew, xnew, ynew, lonPol, latPol)

I = np.argwhere(np.isnan(Ur))
nnan = np.size(I)
if nnan:
	print('Warning: rotated U,V have nans')

# ====================================
#		PLOT
# ====================================
fig1 = plt.figure(1,figsize=(8,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.75, 0.75])

plt.sca(ax1)

# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,360.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

im1 = plt.pcolormesh(xh,yh,Sav,shading='flat',cmap='Wistia')
m.fillcontinents(color=[0.8, 0.8, 0.8])
plt.clim(0., 5.)

plt.title('JRA-55, Uwind av:{0}/{2}/{3}-{1}/{2}/{3}'.format(dS,dE,mo,yr))

xl1 = 22000.
xl2 = 6000.e3
yl1 = 2200.e3
yl2 = 6000e3

plt.axis('equal')
plt.xlim([xl1, xl2])
plt.ylim([yl1, yl2])



# Colorbar
# Create axes for colorbar
ax2 = fig1.add_axes([ax1.get_position().x1+0.02, ax1.get_position().y0,0.02,
                    ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2)
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=16)

#for lbl in ticklabs:
#	print(lbl)

# Vectors:
plt.sca(ax1)
dp=2
#plt.quiver(xnew[::dp, ::dp],ynew[::dp, ::dp],Ur[::dp, ::dp],Vr[::dp, ::dp],
#               units='x', width=0.2, scale=0.4)
import mod_draw_vector
importlib.reload(mod_draw_vector)
from mod_draw_vector import draw_arrowF

nx = Ur.shape[1]
ny = Ur.shape[0]
dx = xnew[0,1]-xnew[0,0]
scl = 1.5*dx
cf_ahd = 0.3
larr_min = cf_ahd*scl*2.5
larr_max = cf_ahd*scl*6.0
vec_props={
	'v_col': [0, 0.2, 0.6],
	'lwd': 2.0,
	'cf_ahd': cf_ahd,
	'beta': 25,
	'larr_min': larr_min,
	'larr_max': larr_max
}

for jj in range(0,ny,dp):
	for ii in range(0,nx,dp):
		x1 = xnew[jj,ii]
		y1 = ynew[jj,ii]
		if x1 < xl1 or x1 > xl2 or y1 < 1600*1.e3 or y1 > 6500*1.e3:
			continue

		u1 = Ur[jj,ii]
		v1 = Vr[jj,ii]		
		x2 = x1+scl*u1
		y2 = y1+scl*v1

		draw_arrowF(x1,x2,y1,y2,**vec_props)

# Streamlines:
#m.streamplot(xnew,ynew,Ur,Vr,density=5)

plt.xlim([xl1, xl2])
plt.ylim([yl1, yl2])

btx = 'plot_jra_winds.py'
import mod_utils_fig
from mod_utils_fig import bottom_text
bottom_text(btx)


# Plot original fields, on JRA grid
f_plt0 = False
if f_plt0:
	fig2 = plt.figure(figsize=(14,9))
	plt.clf()

	ax1 = fig2.add_subplot(2,1,1)
	ax2 = fig2.add_subplot(2,1,2)

	im1 = ax1.pcolormesh(lonc,latc,Sav, shading='flat', cmap='GnBu')
	im1.set_clim(0, 5.)
	dp=3
	Q = ax1.quiver(lon2d_in[::dp, ::dp],lat2d_in[::dp, ::dp],Uav[::dp, ::dp],Vav[::dp, ::dp], 
								 units='x', width=0.2, scale=0.4)
	ax1.axis('equal')
	ax1.set_xlim([0, 180])
	ax1.set_ylim([60, 90])

	im2 = ax2.pcolormesh(lonc,latc,Sav, shading='flat', cmap='GnBu')
	im2.set_clim(0, 5.0)
	Q = ax2.quiver(lon2d_in[::dp, ::dp],lat2d_in[::dp, ::dp],Uav[::dp, ::dp],Vav[::dp, ::dp],    
								 units='x', width=0.2, scale=0.4)
	ax2.axis('equal')
	ax2.set_xlim([180, 360])
	ax2.set_ylim([60, 90])








