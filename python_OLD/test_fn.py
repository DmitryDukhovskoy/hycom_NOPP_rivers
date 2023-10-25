#   function dist = distance_spheric_coord(LAT1,LON1,LAT2,LON2)
#   units = m
#  this procedure calculates the great-circle distance between two
#  geographical locations on a spheriod given it
#  lat-lon coordinates with its appropiate trigonometric
#  signs. 
#  INPUT: xla1, xlo1 - first point coordinates (latitude, longitude)
#         xla2, xlo2 - second point
# all input coordinates are in DEGREES: latitude from 90 (N) to -90,
# longitudes: from -180 to 180 or 0 to 360,
# LAT2, LON2 can be either coordinates of 1 point or N points (array)
# in the latter case, distances from Pnt 1 (LAT1,LON1) to all pnts (LAT2,LON2)
# are calculated 
#  OUTPUT - distance (in m)
#  R of the earth is taken 6371.0 km
#
import numpy as np

def dist_sphcrd(xla1,xlo1,xla2,xlo2):
#  if np.absolute(xla1) > 90.0:
#    print("ERR: dist_sphcrd Lat1 > 90")
#    dist = float("nan")
#    return dist
#  if np.absolute(xla2) > 90.0:
#    print("ERR: dist_sphcrd Lat2 > 90")
#    dist = float("nan")
#    return dist
#
  dist = 0.0
  R = 6371.0e3
  cf = np.pi/180
  phi1 = xla1*cf
  phi2 = xla2*cf
  lmb1 = xlo1*cf
  lmb2 = xlo2*cf
  dphi = phi2-phi1
  dlmb = lmb2-lmb1
#
#  Central angle between 2 pnts:
#  dmm1 = (np.cos(phi1)*np.sin(dlmb))**2
#  dmm2 = (np.cos(phi2)*np.sin(phi1)-np.sin(phi2)*np.cos(phi1)*np.cos(dlmb))**2
#  dmm3 = np.absolute(np.sin(phi2)*np.sin(phi1)+np.cos(phi2)*np.cos(phi1)*np.cos(dlmb))
#  dsgm = np.arctan(np.sqrt((dmm1+dmm2)/dmm3))
#
# The great-circle distance
#  dist1 = R*dsgm
#
# Another formula:
#  dmm = np.sin(phi1)*np.sin(phi2)+np.cos(phi1)*np.cos(phi2)*np.cos(lmb2-lmb1)
#  if dmm>1.0:
#    dmm=1.0
#  dist2 = np.acos(dmm)*R
#
#  dd = np.absolute(1-dist1./dist2)
#  if dd>1.0e-3:
#    dist = dist1
#  else:
#    dist=0.5*(dist1+dist2)

  return R

