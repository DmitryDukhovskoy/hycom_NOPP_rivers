% Read tbaric fields
% thermobaric reference state
% default 2, Arctic = 1 transition
% at GIN sill
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.08';

pthin   = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/relax/110/41layers_T11/',R);
%pthout  = '/Net/mars/ddmitry/hycom/ARCc0.04/relax/010/';
pthtopo = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);

T  = 11;
TV = sprintf('%2.2i',T);

% Get new topo and grid:
fltopo_new=sprintf('%sdepth_%s_%s.nc',pthtopo,R,TV);
HH   = nc_varget(fltopo_new,'Bathymetry');
alat = nc_varget(fltopo_new,'Latitude');
elon = nc_varget(fltopo_new,'Longitude');
LAT  = alat;
LON  = elon;
[mm,nn]= size(HH);

IDM=nn; % ARCc0.08 (old) grid
JDM=mm; % ARCc0.08
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);



knew=41;

% Read in sigma layers:
fina = sprintf('%stbaric_T%s.a',pthin,TV);
finb = sprintf('%stbaric_T%s.b',pthin,TV);

fid = fopen(fina,'r');
dmm = fread(fid,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
A   = reshape(dmm,IDM,JDM)';
A(A>1e10)=nan;
fclose(fid);

