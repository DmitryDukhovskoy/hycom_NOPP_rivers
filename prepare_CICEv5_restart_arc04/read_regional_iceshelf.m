% CICE uses HYCOM's "ice shelf" bathymetry (land/sea boundary)
% gordon02 123> pwd
%/app/projects/hycom/GLBc0.04/topo
%gordon02 124> ll *11is27*
%-rw-r----- 1 shriver 0375G018 253984768 May 11  2018 depth_GLBc0.04_11is27.a
%-rw-r----- 1 shriver 0375G018       434 May 11  2018 depth_GLBc0.04_11is27.b
%-rw-r----- 1 shriver 0375G018 761940356 Apr 14 14:39 hycom_latlonmask_11is27.nc
%lrwxrwxrwx 1 shriver 0375G018        22 Mar 22  2017 regional.cice_11is27.r -> regional.cice_11is17.r 
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

pthice = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo/';


TV = '17DD';  % topo name
ntopo = 17;
PTH.data = '/Net/mars/ddmitry/hycom/ARCc0.04/force/rivers/';
PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
PTH.dataout = '/Net/mars/ddmitry/hycom/ARCc0.04/force/riversNCAR/';
%flt=[pth,'regional.grid.b'];


fltopo = sprintf('%sdepth_ARCc0.04_%s.nc',PTH.topo,TV);
fprintf('Year %i, Mo=%i\n',YearGr,MM);

% Get topo and grid:
HH   = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn]= size(HH);
IDM = nn;
JDM = mm;


% read bathymetry from cice.regional.r
%
% in grid2cice.f (located at:
% /pacific/abozec/HYCOM/hycom/ALL2/cice/src
% ice grid is written as direct access array:
% each record is IDMo x JDMo double precision (real*8)
%c        kmt    land mask array (0,1)
%c        ulati  latitude  of u-cell centers (radians)
%c        uloni  longitude of u-cell centers (radians)
%c        htn    length of northern edge of t-cell (m)
%c        hte    length of eastern  edge of t-cell (m)
%c        anglet conversion on t-cell between cice and lat-long grids (radians)
%c        tlati  latitude  of t-cell centers (radians)
%c        tloni  longitude of t-cell centers (radians)
%
fgrd = sprintf('%sregional.cice.T17DD.r',pthice);
fgri = fopen(fgrd,'r');
frewind(fgri);
dmm=fread(fgri,1,'float64','ieee-be');
kmt=fread(fgri,[IDM,JDM],'float64','ieee-be');
kmt=kmt';
dmm=fread(fgri,1,'float64','ieee-be');
fclose(fgri);


% Plot sea ice land mask:
figure(1); clf;
pcolor(kmt); shading flat;






