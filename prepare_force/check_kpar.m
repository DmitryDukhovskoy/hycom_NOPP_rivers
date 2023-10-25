% Check SeaWifs ocean color
% kpar.a to ARCc0.04
% also check range in created fields by using Alan's code:
% /nexsan/people/ddmitry/Net_ocean/HYCOM/hycom/ALL4/bin/hycom_range kpar_arc04.a 3200 5040
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

pthin   = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/force/seawifs/';
pthout  = '/Net/mars/ddmitry/hycom/ARCc0.04/force/';
pthtopo2= '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
fltopo2 = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo2); % 

fprintf('Reading topo %s\n',fltopo2);
HH2  = nc_varget(fltopo2,'Bathymetry');
LAT2 = nc_varget(fltopo2,'Latitude');
LON2 = nc_varget(fltopo2,'Longitude');
[m2,n2]= size(HH2);
IDM2=n2;
JDM2=m2;
IJDM2=IDM2*JDM2;
npad2=4096-mod(IJDM2,4096);
toto=ones(npad2,1);
toto2 = toto;

fouta = sprintf('%skpar_arc04.a',pthout);
foutb = sprintf('%skpar_arc04.b',pthout);
fida = fopen(fouta,'r');
%fidb = fopen(foutb,'r');

fprintf('Reading fields: %s\n',fouta);

imp = 12; % month to plot


frewind(fida);
stat = fseek(fida,(imp-1)*(IJDM2+npad2)*4,-1);
dmm=fread(fida,IJDM2,'float32','ieee-be');  % read 2D field (1 layer)
dm1=fread(fida,npad2,'float32','ieee-be');  % read npad 

F = reshape(dmm,IDM2,JDM2)';

