#!/bin/csh
## to make the CICE grid from regional grid.[ab]
setenv R ARCc0.04
setenv TV 17DD
setenv T /Net/ocean/ddmitry/HYCOM/ARCc/${R}/topo
touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
  /bin/rm regional.depth.a
  /bin/ln -s ${T}/depth_${R}_${TV}.a regional.depth.a
endif
if (-z regional.depth.b) then
  /bin/rm regional.depth.b
  /bin/ln -s ${T}/depth_${R}_${TV}.b regional.depth.b
endif
#
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/ln -s ${T}/regional.grid.a regional.grid.a
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -s ${T}/regional.grid.b regional.grid.b
endif
touch      fort.30 fort.31 fort.32 
/bin/rm -f fort.30 fort.31* fort.32*
#
# --- D,y,d select the archive files.
/Net/ocean/ddmitry/HYCOM/hycom/ALL4/cice/src/grid2cice <<E-o-D
   1    'i1st ' = 1st hycom i-point on cice grid
   1    'j1st ' = 1st hycom j-point on cice grid
3200    'imt   ' = 1st cice global array dimension
5040    'jmt   ' = 2nd cice global array dimension
E-o-D
 
setenv FLI regional.cice.T${TV}.r
touch regional.cice.r
if ( -z regional.cice.r ) then
  echo "File regional.cice.r was not created "
  exit 1
else
  echo "moving regional.cice.r --> ${T}/${FLI}"
  mv regional.cice.r ${T}/${FLI}
endif

rm -f regional.cice.txt

exit 0
## 
 
