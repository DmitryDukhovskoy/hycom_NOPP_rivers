% Plot sea ice fields from 0.04 HYCOM-CICE
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt=012;
yr=2005;
mo=7;
dm=3;

%pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthbin  = '/nexsan/people/ddmitry/hycom/ARCc0.08/110/cice/';
pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i_cice/',expt,yr);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
btx='plot_ice_arc04.m';

ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);

Lmsk=HH*0;
Lmsk(HH<0)=1;

fin = sprintf('%s%3.3i_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',...
	      pthbin,expt,yr,mo,dm);

%A=squeeze(nc_varget(fin,'divu'));
%u=squeeze(nc_varget(fin,'uocn'));
%v=squeeze(nc_varget(fin,'vocn'));
%u=squeeze(nc_varget(fin,'uvel'));
%v=squeeze(nc_varget(fin,'vvel'));
%S=sqrt(u.^2+v.^2);

%sairx=squeeze(nc_varget(fin,'strairx')); % N/m2
%sairy=squeeze(nc_varget(fin,'strairy'));
%S=sqrt(sairx.^2+sairy.^2);

%pcolor(S); shading flat

%hs = squeeze(nc_varget(fin,'hs')); % snow
hc = squeeze(nc_varget(fin,'aice')); % 
hi = squeeze(nc_varget(fin,'hi'));


%i
stl=sprintf('0.04 HYCOM-%3.3i, %4.4i/%2.2i/%2.2i, Ice Thickness, m',expt,yr,mo,dm);
xl1=200;
xl2=3200;
yl1=600;
yl2=4000;
c1=0;
c2=5;
fn=15;
sub_plot_hice(hi,Lmsk,xl1,xl2,yl1,yl2,c1,c2,fn);
contour(hc,[0.15 0.15],'r-');
title(stl);

bottom_text(btx,'pwd',1);

%
% Histogram of ice thicknesses
nff=2;
hbin=[0.5:10];
sub_hist_hice(hi,hbin,nff);
title(stl);

bottom_text(btx,'pwd',1,'position',[0.08 0.45 0.3 0.04]);





