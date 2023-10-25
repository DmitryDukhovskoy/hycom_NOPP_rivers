#!/bin/csh -x
# 
# Remap nest archive files
# from GLBb0.08 GOFS3.1 reanalysis
# GOFS3.1 uses T11 (41 L)
# In Step 1: GLBb GOFS3.1 -> ARCc0.08
# see: ../remap_gridGLB2ARC_archm/
# Both domain have same T11
# however this topo does not
# match my T11 for ARCc0.08 - 
# the Pacific OB is closed
# Also converted fields are not "closed" at the Pacific OB
#
# Step2: Thus need to reinterploate into
# my T11 topo
# 
# Note that the input
# archive file has to be on the same grid/domain
# as the output grid, the output grid 
# can be a finer-grid but must be an integer
# multiple of the original grid
#
setenv DF /Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files
setenv DT /Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid

touch regional.grid.a
/bin/rm -f regional.grid.[ab]
/bin/ln -s ${DT}/regional.grid.a # ARCc0.08
/bin/ln -s ${DT}/regional.grid.b

setenv flnm_in  ${DF}/archv_arcT11L41_stp1.2015_130_00
setenv flnm_tin ${DT}/depth_ARCc0.08_GOFS3.1_11
setenv flnm_out ${DF}/archv_arcT11L41_stp2.2015_130_00
setenv flnm_top ${DT}/depth_ARCc0.08_11

# Check existing files, the code cannot overwrite it
/bin/rm -f ${flnm_out}.[ab]

# Input:
#c --- 'flnm_in'   = input  filename
#c --- 'flnm_tin'  = input  bathymetry filename
#c --- 'flnm_out'  = output filename
#c --- 'flnm_top'  = output bathymetry filename
#c --- 'cline_out' = output title line (replaces preambl(5))
#/Net/ocean/ddmitry/HYCOM/hycom/ALL/subregion/src/isubregion <<E-o-D
/Net/ocean/ddmitry/HYCOM/hycom/ALL4/subregion/src/isubregion <<E-o-D
${flnm_in}.a
${flnm_tin}.a
${flnm_out}.a
${flnm_top}.a
From GLBb0.08_19.0 T07 -> ARCc0.08 T07 -> ARCc0.08 T11 grid
  1600	     'idm   ' = longitudinal array size
  2520	     'jdm   ' = latitudinal  array size
    1        'irefi ' = longitudinal input  reference location
    1        'jrefi ' = latitudinal  input  reference location
    1        'irefo ' = longitudinal output  reference location
    1        'jrefo ' = latitudinal  output  reference location
    1        'ijgrd ' = integer scale factor between input and output grids
    1        'iceflg' = ice in output archive flag (0=none,1=energy loan model)
    0        'smooth' = smooth interface depths (0=F,1=T)
E-o-D

exit 0
