#!/bin/csh
#
set echo
#
# --- interpolate to 3-d z-levels from a single HYCOM archive file.
# --- z-levels, via linear interpolation, at Levitus depths.
#
# --- output is netCDF.
# --- this is an example, customize it for your datafile needs.
#
# --- optional title and institution.
#
setenv DT /Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid
setenv ALL /Net/ocean/ddmitry/HYCOM/hycom/ALL4

touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
  /bin/rm regional.depth.a
  /bin/ln -s $DT/depth_ARCc0.08_11.a regional.depth.a
endif
if (-z regional.depth.b) then
  /bin/rm regional.depth.b
  /bin/ln -s $DT/depth_ARCc0.08_11.b regional.depth.b
endif
#
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/ln -s $DT/regional.grid.a regional.grid.a
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -s $DT/regional.grid.b regional.grid.b
endif
#
# --- optional title and institution.
#
setenv CDF_TITLE        "0.08 HYCOM-CICE, subset NAtl no Greenland runoff"
setenv CDF_INST         "COAPS FSU"
#
#
# --- D,y,d select the archive files.
#
setenv YY 2007
setenv D /nexsan/archive/ARCc0.08_110/data/${YY}
setenv O /nexsan/archive/ARCc0.08_110/data_natl/${YY}

if (! -d ${O}) then
  mkdir -pv ${O}
endif

#foreach d ( a )
@ d1=1
while ($d1 < 367)
    setenv d `echo ${d1} | awk '{printf("%03d", $1)}'`
    setenv CDF030 ${O}/hycom008_110z_NAtl_${YY}-${d}_t.nc
    setenv CDF031 ${O}/hycom008_110z_NAtl_${YY}-${d}_s.nc
#    set ARR = (CDF033 CDF034 CDF037 CDF 038)
#    /bin/rm -f $CDF033.bin  $CDF034.bin $CDF037.bin  $CDF038.bin 
    /bin/rm -f $CDF030 $CDF031
#${ALL}/archive/src/archv2data3z <<E-o-D
${ALL}/archive/src/archv2ncdf3z <<E-o-D
${D}/110_archm.${YY}_${d}_12.a
netCDF
 000    'iexpt ' = experiment number x10 (000=from archive file)
   0    'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 1600   'idm   ' = longitudinal array size
 2520   'jdm   ' = latitudinal  array size
  41    'kdm   ' = number of layers
  34.0  'thbase' = reference density (sigma units)
   0    'smooth' = smooth the layered fields (0=F,1=T)
 400    'iorign' = i-origin of plotted subregion
 300    'jorign' = j-origin of plotted subregion
 650    'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
 800    'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
   1    'itype ' = interpolation type (0=sample,1=linear)
   55    'kz    ' = number of depths to sample
    0.0  'z     ' = sample depth  1
    1.0  'z     ' = sample depth  2
    2.0  'z     ' = sample depth  3
    4.0  'z     ' = sample depth  4
    6.0  'z     ' = sample depth  5
    8.0  'z     ' = sample depth  6
   10.0  'z     ' = sample depth  7
   12.0  'z     ' = sample depth  8
   14.0  'z     ' = sample depth  9
   16.0  'z     ' = sample depth 10
   18.0  'z     ' = sample depth 11
   20.0  'z     ' = sample depth 12
   25.0  'z     ' = sample depth 13
   30.0  'z     ' = sample depth 14
   35.0  'z     ' = sample depth 15
   40.0  'z     ' = sample depth 16
   45.0  'z     ' = sample depth 17
   50.0  'z     ' = sample depth 18
   60.0  'z     ' = sample depth 19
   70.0  'z     ' = sample depth 20
   80.0  'z     ' = sample depth 21
   90.0  'z     ' = sample depth 22
  100.0  'z     ' = sample depth 23
  125.0  'z     ' = sample depth 24
  150.0  'z     ' = sample depth 25
  175.0  'z     ' = sample depth 26
  200.0  'z     ' = sample depth 27
  225.0  'z     ' = sample depth 28
  250.0  'z     ' = sample depth 29
  275.0  'z     ' = sample depth 30
  300.0  'z     ' = sample depth 31
  325.0  'z     ' = sample depth 32
  350.0  'z     ' = sample depth 33
  400.0  'z     ' = sample depth 34
  500.0  'z     ' = sample depth 35
  600.0  'z     ' = sample depth 36
  700.0  'z     ' = sample depth 37
  800.0  'z     ' = sample depth 38
  900.0  'z     ' = sample depth 39
 1000.0  'z     ' = sample depth 40
 1250.0  'z     ' = sample depth 41
 1500.0  'z     ' = sample depth 42
 1750.0  'z     ' = sample depth 43
 2000.0  'z     ' = sample depth 44
 2250.0  'z     ' = sample depth 45
 2500.0  'z     ' = sample depth 46
 2750.0  'z     ' = sample depth 47
 3000.0  'z     ' = sample depth 48
 3250.0  'z     ' = sample depth 49
 3500.0  'z     ' = sample depth 50
 3750.0  'z     ' = sample depth 51
 4000.0  'z     ' = sample depth 52
 4500.0  'z     ' = sample depth 53
 5000.0  'z     ' = sample depth 54
 6000.0  'z     ' = sample depth 55
  0     'botio ' = bathymetry  I/O unit (0 no I/O)
  0     'mltio ' = mix.l.thk.  I/O unit (0 no I/O)
  0     'tempml' = temperature jump across mixed-layer (degC,  0 no I/O)
  0     'densml' =   density jump across mixed-layer (kg/m3, 0 no I/O)
  0     'infio ' = intf. depth I/O unit (0 no I/O, <0 label with layer #)
  0     'wvlio ' = w-velocity  I/O unit (0 no I/O)
  0    'uvlio ' = u-velocity  I/O unit (0 no I/O)
  0    'vvlio ' = v-velocity  I/O unit (0 no I/O)
  0     'splio ' = speed       I/O unit (0 no I/O)
  30    'temio ' = temperature I/O unit (0 no I/O)
  31    'salio ' = salinity    I/O unit (0 no I/O)
  0     'tthio ' = density     I/O unit (0 no I/O)
  0     'keio  ' = kinetic egy I/O unit (0 no I/O)
E-o-D

@ d1++
echo ${d1}
# Move output files:

end
chmod -R 755 ${O}
