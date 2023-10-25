% Plot results 
% Heat/Vol fluxes 
% HtVlm_flux_contour_mean_season.m

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=1993;
YR2=2009;

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat2/';
btx='plot_HtVlm_flux_contour.m';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section

finp=sprintf('%sarc08_%3.3i_Greenl_VolHeat_seasons_%i-%i.mat',pthmat,expt,YR1,YR2);
fprintf('Loading %s\n',finp);
load(finp);

cff = 1e-6;
%im  = 1;
fn  = 1;
hf1 = VHFLX.Hflx_tot*cff;
xl1 = 450;
xl2 = 1050;
yl1 = 380;
yl2 = 1100;
c1  = -100;
c2  = 100;
stl = sprintf('Mean HeatFlx, W/m*%3.1d, %i-%i',cff,YR1,YR2);
sub_plot_hflx1d(HH,LON,LAT,HFLX,hf1,cff,fn,stl,xl1,xl2,yl1,yl2,c1,c2);
bottom_text(btx,'pwd',1);
drawnow




