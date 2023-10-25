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

YR1=2000;
YR2=2016;

f_map=0; % fluxes color-coded along contour, map, UV

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthfig=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_green_xsct/',...
		  regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';
btx = 'plot_TrFlx_contour.m';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

xyl(1,1) = 250;
xyl(1,2) = 900;
xyl(2,1) = 1600;
xyl(2,2) = 2000;

fn = 1;
sub_plot_bath3(HH,LON,LAT,fn,xyl);
btx = 'plot_Arctic_domain.m';
bottom_text(btx,'pwd',1);

