Preapre hycom ARCc0.08 for NOPP experiments
HYCOM: 2.2.99Pi-900_relo_cice_v.40
17 term eq. of state, sigma 2

Vertical layers: 41
Topography T11
Note that after remapping Global T11 -> ARCc T11
Pacific OB remains open, whereas it should be "closed", i.e.
there should be wall (depths >0), see:
hycom_NOPP_rivers/topo/correct_ARCc_T11.m


Use restart: from GLBb0.08 expt_19.0: 1993 01 01
Model is forced by CFSR and CFSv2
Ocean: GLBb0.08: 19.0, 19.1
GLBb0.08: Topo T11, 32 layers, 17-term sigma-2 eq. of state

Global grid is remaped into ARCc
vertical layers are remaped into 41 layers

==================================================
newton 13> pwd
/u/home/wallcraf/hycom/ARCb0.08/topo
newton 14> ll
total 626704
drwxr-xr-x   2 wallcraf 0375G018    4096 Apr  1 20:52 ./
drwxr-xr-x   3 wallcraf 0375G018    4096 Apr  1 20:51 ../
-rw-r--r--   1 wallcraf 0375G018 16138240 Apr  1 20:51 depth_ARCb0.08_11.a
-rw-r--r--   1 wallcraf 0375G018     441 Apr  1 20:51 depth_ARCb0.08_11.b
-rw-r--r--   1 wallcraf 0375G018 258048000 Apr  1 20:51 regional.cice_11.r
-rw-r--r--   1 wallcraf 0375G018 32276480 Apr  1 20:52 regional.gmapi_GLBb0.08.a
-rw-r--r--   1 wallcraf 0375G018     127 Apr  1 20:52 regional.gmapi_GLBb0.08.b

The gmapi file is a copy of GLBa0.08 to ARCc0.08, since the actual grid 
does not change between GLBa0.08 and GLBb0.08 (and ARCc0.08 and ARCb0.08).

The reanalysis uses the 17-term sigma2 equation of state, but is 32 layers. 
Since the extra 9 layers are fixed depth and near the surface 
it is easy to convert an archive file from 32 to 41 layers.

Alan.

On 04/04/2016 10:28 AM, Dukhovskoy, Dmitriy wrote:
could you point me to the blkdat file with 41 layers?

On shepard, see:

shepard03 168> pwd
/p/home/wallcraf/hycom/GLBb0.08/expt_69.1
shepard03 169> ll blk*t ice_in
-rw-r----- 1 wallcraf 0375G018 18273 Aug 24  2015 blkdat.input
-rw-r----- 1 wallcraf 0375G018  4791 Aug 24  2015 ice_in

Note that of you run with hourly atmospheric 
forcing (CFSR) you need to specify cfsr in the ice_in file.

Alan.
=================================

Prepare T11 topography:
use ARCb0.08_11 frmo Alan, the 1st and lsat rows are not walls
due to big mismatch in depth between T07 and T11 at last row
wierd density fields are generated during interpolatino of 
input fields from T07 -> T11
first, need to close off the N/S boundaries at j 
see ../topo/correct_ARCc_T11.m -> *ARCc0.08_11*[ab] and nc files


Need to prepare restart files on ARCc0.08 grid, T11, 41 layers
steps:
1) Remap global archive files to ARCc grid
in directory:
anls_mtlb_utils/hycom_arc08/NOPP_rivers/REMAP_ARCc/remap_gridGLB2ARC_archv
   Use 190_archv.1993_001_00.[ab] from GLBb0.08 expt_19.0 (32 lr, T07)
   -> archv_arcT07.1993_001_00 Topo T07, 32 layers

  Use Fortran code: remapGLB.x (remap_GLB2ARC.F90) for remapping



2) Convert arc08_T07_32l_archv.1993_001_00 -> arc08_T11_32l_archv.1993_001_00
  ARCc0.08 T07 -> ARCc0.08 T11 32 layers
run script isubregion_T07_T11.csh
(output: archv_arcT11.1993_001_00.a)

May try to interpolate these files to 41 layers assuming 
added layers are fixed-depth 
see code:
hycom_arc08/NOPP_rivers/REMAP_ARCc/remap32to41lrs_fixedZ.m

read README.remap in that directory  (REMAP_ARCc/):
   - interpolate 32 -> 41 layers
   and prepare nest files and restart HYCOM files



Alternatively (continue from step 2):
3) Interpolate into 41 layers 
   arc08_T11_32l_archv.1993_001_00->arc08_T11_41l_archv.1993_001_00 

4) Run HYCOM in new configuration starting from rest with no atm / river forcing
   to generate a template restart file (see expt_10.1)

5) use this restart template and interpolated archive (41 layers, T11)
   to create a new restart:
/home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/REMAP_ARCc/...
 archv2restart_41lr_T11.com


