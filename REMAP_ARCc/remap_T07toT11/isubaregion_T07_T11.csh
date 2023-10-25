Use isubregion_*.csh instead of this script
in isubregion need to subset the region into ARCc from
GLBb first (use remap_GLB2ARC.F90 for this)

THis script:
Could not make this scrsipt work
in theory it should work for:
take GLBb archive file
provide with gmap file (remapping indices to
remap GLBb -> ARCc)
provide topo files for GLBb and ARCc 
and the code will subset the region and 
interpolate  it into new topo/bathymetry
Gives an error in reading topo file: max value from 
*a file does not match *b

#!/bin/csh -x
# 
# Remap archive file from
# GLBb to ARCc using
# map file
# and from
# Topo T07 -> T11
# 
setenv DF /Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files
setenv DG /Net/kronos/ddmitry/hycom/GLBb0.08/expt_19.0/output
setenv DT /Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid

touch regional.grid.a
/bin/rm -f regional.grid.[ab]
/bin/ln -s ${DT}/regional.grid.a # ARCc0.08
/bin/ln -s ${DT}/regional.grid.b
#
# See isubaregion.f:
#c --- 'flnm_reg'  = target sub-region grid       filename
#c --- 'flnm_map'  = target sub-region grid map   filename
#c --- 'flnm_top'  = target bathymetry filename, or 'NONE'
#c --- 'flnm_tin'  = input  bathymetry filename, or 'NONE'
#c --- 'flnm_in'   = input  archive    filename
#c --- 'flnm_out'  = output archive    filename
#c --- 'cline_out' = output title line (replaces preambl(5))
#/Net/ocean/ddmitry/HYCOM/hycom/ALL/subregion/src/isubaregion <<E-o-D
/Net/ocean/ddmitry/HYCOM/hycom/ALL4/subregion/src/isubaregion <<E-o-D
${DT}/regional.grid.a
${DT}/regional.gmapi_GLBb0.08.a
${DT}/depth_ARCc0.08_11.a
${DT}/depth_GLBb0.08_07.a
${DF}/archv_arcT07.1993_001_00.a
${DG}/190_archv.1993_001_00.a
GLBb0.08 T07 to ARCc0.08 T11 grid
  1600	     'idm   ' = longitudinal array size
  2520	     'jdm   ' = latitudinal  array size
    0        'iceflg' = ice in output archive flag (0=none,1=energy loan model)
    0        'smooth' = smooth interface depths (0=F,1=T)
E-o-D

exit 0
