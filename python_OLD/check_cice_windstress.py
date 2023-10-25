"""
 Check sea ice winds/stresses
"""
from IPython import get_ipython
get_ipython().magic('reset -sf') # to clean memory when script is run

from netCDF4 import Dataset as NetCDFFile
import numpy as np
from mod_utils_fig import bottom_text
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import datetime
from mpl_toolkits.basemap import Basemap, cm
from copy import copy

plt.close('all')

yr   = 2017
mo   = 7
mday = 15
hh   = 0

def get_grid():
  pthgrid = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/'
  flgrd = 'depth_ARCc0.04_17DD.nc'
  nc = NetCDFFile(pthgrid+flgrd)
  lonc = nc.variables['Longitude'][:]
  latc = nc.variables['Latitude'][:]
  lonc[lonc<0]=lonc[lonc<0]+360

  return latc,lonc

latc, lonc = get_grid()


pthnc = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_cice/'
date = datetime.datetime(yr,mo,mday,hh)

#smo = str(mo).zfill(2)
#smday = str(mday).zfill(2)
#flnm  = '022_cice.{0}-{1}-{2}.nc'.format(yr,smo,smday)
flnm = '022_cice.{0}-{1:02d}-{2:02d}.nc'.format(yr,mo,mday)

huge = 1.e20
nc = NetCDFFile(pthnc+flnm)
uatm = nc.variables['uatm'][:].squeeze()
uatm = np.array(uatm)
vatm = nc.variables['vatm'][:].squeeze()
vatm = np.array(vatm)

uatm[uatm>huge] = np.nan
vatm[vatm>huge] = np.nan

S = np.sqrt(uatm**2+vatm**2)

# Set up a colormap:
clrmp = copy(plt.cm.OrRd)
clrmp.set_under(color=[0.9, 0.9, 0.9])
clrmp.set_bad(color=[0.3, 0.3, 0.3])

# create figure and axes instances
plt.ion()  # enables interactive mode

fig1 = plt.figure(figsize=(8,8))
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,360.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

ny = S.shape[0]; nx = S.shape[1]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly spaced grid.
x, y = m(lons, lats) # compute map proj coordinates.
xh, yh = m(lonc,latc)  # hycom coordinates on the projections
#im = m.pcolormesh(xh,yh,aice,shading='flat',cmap=clrmp)
plt.pcolormesh(xh,yh,S,shading='flat',cmap=clrmp,
               norm=colors.Normalize(vmin=0., vmax=15.0))

#
# Plot wind vectors

plt.title('0.04 HYCOM-CICE, aice, '+flnm)
# Colorbar
#plt.colorbar(extend='min')
plt.colorbar(aspect=40)

bottom_text(btx)


plt.show()




