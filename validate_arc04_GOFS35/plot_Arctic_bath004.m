% Plot tracer flux to Greenland
% across specified contour -
% isobath around Greenland
% see trFlux_contour_greenl.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.04';
expt = 22; % GOFS3.5
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
btx = 'plot_Arctic_bath004.m';

ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

xyl(1,1) = 1;
xyl(1,2) = 1;
xyl(2,1) = nn;
xyl(2,2) = mm;

fn = 1;
sub_plot_bath3(HH,LON,LAT,fn,xyl);
bottom_text(btx,'pwd',1);

