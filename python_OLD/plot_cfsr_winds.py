"""
	Plot CFSR winds from CICE output fields
	FCSR winds used as atm. forcing
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

plt.close('all')

# Plot original vectors on HYCOM grid for validation:
f_plt0 = False

# 1 month average
mo  = 6
yr  = 2017
dS  = 1
dE  = 30

print('Plotting av winds {0}/{1}/{2}-{3}/{4}/{5}'.format(dS,mo,yr,dE,mo,yr))
pthnc = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/{0}_cice/'.format(yr)

def get_grid():
	pthgrid = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/'
	flgrd = 'depth_ARCc0.04_17DD.nc'
	nc = NetCDFFile(pthgrid+flgrd)
	lonc = nc.variables['Longitude'][:]
	latc = nc.variables['Latitude'][:]
	lonc[lonc<0]=lonc[lonc<0]+360
	HH = nc.variables['Bathymetry'][:]

	return latc,lonc,HH


latc, lonc, HH = get_grid()
nn = lonc.shape[1]
mm = lonc.shape[0]

Uav = np.zeros((mm,nn))
Vav = np.zeros((mm,nn))

Tang = np.array([])

ndays = dE-dS+1
for ii in range(ndays):
	print("Reading record ",ii)

	mday = ii+1	
	flnm = '022_cice.{0}-{1:02d}-{2:02}.nc'.format(yr,mo,mday)
	ncuv = NetCDFFile(pthnc+flnm)
	uin = ncuv.variables['uatm'][:].squeeze()
	vin = ncuv.variables['vatm'][:].squeeze()

	Uav += uin
	Vav += vin
#
# Get angle 
# angle grid makes with lat line on T grid - center 
#	if Tang.size == 0:
#		Tang = ncuv.variables['ANGLET'][:].squeeze()
 
	ncuv.close()

del uin, vin
Uav = Uav/np.float(ndays)
Vav = Vav/np.float(ndays)
Uav[Uav>1.e20] = np.nan
Vav[Vav>1.e20] = np.nan

# Subsample smaller domain
isb1 = 300
isb2 = nn
jsb1 = 1000
jsb2 = 4000

Uav0 = np.copy(Uav)
Vav0 = np.copy(Vav)

skp = 3 
Uav = Uav0[jsb1:jsb2:skp,isb1:isb2:skp]
Vav = Vav0[jsb1:jsb2:skp,isb1:isb2:skp]
Sav = np.sqrt(Uav**2+Vav**2)
lonc = lonc[jsb1:jsb2:skp,isb1:isb2:skp]
latc = latc[jsb1:jsb2:skp,isb1:isb2:skp]
HH = HH[jsb1:jsb2:skp,isb1:isb2:skp]
#Tang = Tang[jsb1:jsb2:skp,isb1:isb2:skp]
ny = Sav.shape[0]; nx = Sav.shape[1]

Lmsk = np.zeros((ny,nx))
Lmsk[HH<0.] = 1

# Find N.Pole direction or grad(lat)
# on HYCOM/CICE grid:
import mod_utils_map
importlib.reload(mod_utils_map)
from mod_utils_map import hycom_dir_gradNorth
from mod_utils_map import dx_dy_2Dgrid
from mod_utils_map import dist_sphcrd
from mod_utils_map import dphi
import mod_draw_vector
#importlib.reload(mod_draw_vector)
from mod_draw_vector import compass


alf_north_hycom = hycom_dir_gradNorth(lonc,latc)

plt.ion()  # enables interactive mode


# Create polar orthogonal grid
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')

# generate a grid that is equally spaced in a plot with the current pojection
lon_new, lat_new, xnew, ynew = m.makegrid(200,200, returnxy=True)
xh, yh = m(lonc,latc)  # HYCOM-CICE coordinates on the projections
xnew = xnew.astype('float64')
ynew = ynew.astype('float64')

# ============================================
# Project the data onto the new grid
print("Interpolating onto new grid ...")
points = np.column_stack((xh.flatten(),yh.flatten()))
#unew = griddata(points, Uav.flatten(), (xnew, ynew), method='nearest')
#vnew = griddata(points, Vav.flatten(), (xnew, ynew), method='nearest')
unew = griddata(points, Uav.flatten(), (xnew, ynew), method='linear')
vnew = griddata(points, Vav.flatten(), (xnew, ynew), method='linear')
alf_hycom_new = griddata(points, alf_north_hycom.flatten(), 
								(xnew, ynew), method='nearest')
#
# Rotate U, V vectors to match N and E on the map projection
from mod_utils_map import rotate_vector_hycom2polar
#from mod_utils_map import dphi
#from mod_utils_map import rotate_vector

lonPol = lon_new
latPol = lat_new
Ur, Vr, plr_dirN = rotate_vector_hycom2polar(unew,vnew,alf_hycom_new,xnew,ynew,lonPol,latPol)

f_chck = False
if f_chck:
# Polar grid indices:
	ip0 = 109
	jp0 = 113
#	jh = 600 # HYCOM indices
#	ih = 550
	check_vector_rot(Uav,Vav,xh,yh,lonc,latc,unew,vnew,xnew,ynew,lonPol,latPol,ip0,jp0)

	
#

# ====================================
#   PLOT
# ====================================
print("Plotting ...")
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

plt.title('CFSv2, Winds av:{0}/{2}/{3}-{1}/{2}/{3}'.format(dS,dE,mo,yr))

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
# print(lbl)

# Vectors:
plt.sca(ax1)
dp=3
#plt.quiver(xnew[::dp, ::dp],ynew[::dp, ::dp],Ur[::dp, ::dp],Vr[::dp, ::dp],
#               units='x', width=0.2, scale=0.4)
#importlib.reload(mod_draw_vector)
from mod_draw_vector import draw_arrowF

nx = Ur.shape[1]
ny = Ur.shape[0]
dx = xnew[0,1]-xnew[0,0]
scl = 1.8*dx
cf_ahd = 0.3
larr_min = cf_ahd*scl*3.
larr_max = cf_ahd*scl*7.0
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

btx = 'plot_cfsr_winds.py'
import mod_utils_fig
from mod_utils_fig import bottom_text
bottom_text(btx)


# ==================================
# Plot original fields, on JRA grid
# ==================================
if f_plt0:
	fig2 = plt.figure(figsize=(9,9))
	plt.clf()
	ax1 = plt.axes([0.1, 0.1, 0.75, 0.75])
	plt.sca(ax1)

	clrmp2 = copy(plt.cm.Wistia)
	clrmp2.set_bad(color=[0.5, 0.5, 0.5])

	im2 = plt.pcolormesh(Sav, shading='flat', cmap=clrmp2)
	im2.set_clim(0, 5.)

	xlh1 = 0
	xlh2 = 950
	ylh1 = 0
	ylh2 = 1000
	plt.axis('equal')
	plt.xlim([xlh1, xlh2])
	plt.ylim([ylh1, ylh2])

	plt.title('CFSv2 on CICE grid, Winds av:{0}/{2}/{3}-{1}/{2}/{3}'.format(dS,dE,mo,yr))

# Colorbar
# Create axes for colorbar
	ax2 = fig2.add_axes([ax1.get_position().x1+0.02, ax1.get_position().y0,0.02,
											ax1.get_position().height])
	clb = plt.colorbar(im2, cax=ax2)
	ax2.set_yticklabels(ax2.get_yticks())
	ticklabs = clb.ax.get_yticklabels()
	clb.ax.set_yticklabels(ticklabs,fontsize=16)

# Vectors:
	plt.sca(ax1)
	nsx = Sav.shape[1]
	nsy = Sav.shape[0]
	scl = 5.0
	larr_min = cf_ahd*scl*5.
	larr_max = cf_ahd*scl*15.0
	vec_props={
		'v_col': [0, 0.2, 0.6],
		'lwd': 2.0,
		'cf_ahd': cf_ahd,
		'beta': 25,
		'larr_min': larr_min,
		'larr_max': larr_max
	}
	dph = 18

	for jj in range(0,nsy,dph):
		for ii in range(0,nsx,dph):
			x1 = ii
			y1 = jj
			if x1 < 70 or x1 > xlh2 or y1 < 100 or y1 > ylh2:
				continue

			u1 = Uav[jj,ii]
			v1 = Vav[jj,ii]
			x2 = x1+scl*u1
			y2 = y1+scl*v1

			draw_arrowF(x1,x2,y1,y2,**vec_props)

	plt.contour(latc,[50,60, 70, 80],
							colors='k',
							linewidths=1)
	plt.contour(lonc,[0, 45, 90, 135, 180, 225, 270, 315],
							colors='k',
							linewidths=1)

	bottom_text(btx)

#	ax1.axis('equal')
#	ax1.set_xlim([xlh1, xlh2])
#	ax1.set_ylim([ylh1, ylh2])


def check_vector_rot(Uav,Vav,xh,yh,lonc,latc,unew,vnew,xnew,ynew,lonPol,latPol,ip0,jp0):
	importlib.reload(mod_draw_vector)
	from mod_draw_vector import compass
	from mod_draw_vector import draw_arrowF_polar
	from mod_draw_vector import cart2polar


	xp0 = xnew[jp0,ip0]
	yp0 = ynew[jp0,ip0]
# 
# Find closest point on HYCOM grid
	D = np.sqrt((xh-xp0)**2+(yh-yp0)**2)
	a = np.where(D == np.min(D))
	jph = a[0].item()
	iph = a[1].item()

	fig4 = plt.figure(4, figsize=(8,8))
	plt.plot(xh,yh,'b.')
	plt.plot(xp0,yp0,'.',color=[0,1,0]) # Point on polar grid
	plt.plot(xh[jph,iph],yh[jph,iph],'r.')
	plt.xlim([xp0-10000., xp0+10000.])
	plt.ylim([yp0-10000., yp0+10000.])

# Original vector not interpolated on HYCOM grid
	uh = Uav[jph,iph]
	vh = Vav[jph,iph]

# Interpolated not-rotated vector on polar grid:
	up_nort = unew[jp0,ip0]
	vp_nort = vnew[jp0,ip0]

# Rotated on polar grid:
	u_rot = Ur[jp0,ip0]
	v_rot = Vr[jp0,ip0]

# angle to North on HYCOM grid:
	alfH = alf_north_hycom[jph,iph]
# HYCOM angle to North on Polar grid
	alfHnew = alf_hycom_new[jp0,ip0] # should be close to alfH
# N. Pole dir vector on HYCOM:
	grdNx = np.cos(alfHnew*np.pi/180.)
	grdNy = np.sin(alfHnew*np.pi/180.)

# angle to N. on Polar grid:
	alfPlr = plr_dirN[jp0,ip0]
	grdNPx = np.cos(alfPlr*np.pi/180.)
	grdNPy = np.sin(alfPlr*np.pi/180.)


	nf = 3
	vclr = [0,1,0]
# Vector on HYCOM grid
	ax0 = compass(uh,vh,nf,v_col=vclr)  # not-rotated velocity vector
# N. direction on HYCOM:
	ax0 = compass(grdNx,grdNy, nf, v_col=[0,0.5,0.8],  ax=ax0)

# Not-rotated vector on Polar grid - should be similar to uh,vh
	ax0 = compass(up_nort,vp_nort,nf, v_col=[1,0,0], ax=ax0)

# N. direction on polar grid:
	ax0 = compass(grdNPx,grdNPy, nf, v_col=[0.6,0.4,0], ax=ax0)

# Rotated on polar grid:
	ax0 = compass(u_rot,v_rot, nf, v_col=[0,1,0], ax=ax0)

	return




