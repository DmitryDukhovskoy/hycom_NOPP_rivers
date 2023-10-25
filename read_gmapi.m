% Some file Alan W. Sent:
% "The gmapi file is a copy of GLBa0.08 to ARCc0.08, 
% since the actual grid does not change between GLBa0.08 
% and GLBb0.08 (and ARCc0.08 and ARCb0.08)."
% These are remapping indices GLBa-> ARCc
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

PTH.topo='/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
fltopo=sprintf('%sdepth_ARCc0.08_11.nc',PTH.topo);
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LAT  = alat;
LON  = elon;
[m,n]= size(HH);
[mm,nn]= size(HH);

IDM=nn;
JDM=mm;
IJDM = IDM*JDM;
npad=4096-mod(IJDM,4096);

Lrec=(IJDM+npad)*4; % 1 record length, bytes

fina = sprintf('%sregional.gmapi_GLBb0.08.a',PTH.topo);
fid  = fopen(fina,'r','ieee-be');

xmap = fread(fid,IJDM,'float32');
xmap = reshape(xmap,[IDM JDM])';
fseek(fid,4*(npad+IJDM),-1);
ymap = fread(fid,IJDM,'float32');
ymap = reshape(ymap,[IDM JDM])';

fclose(fid);


