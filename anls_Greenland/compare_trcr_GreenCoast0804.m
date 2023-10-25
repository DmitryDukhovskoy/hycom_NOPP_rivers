% Compare tracer concentration
% along Greenland coast
% data extracted in trcr_GreenCoast*.m
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.04';
pthfig = sprintf('/Net/mars/ddmitry/hycom/ARCc0.04/011/fig_trac/');
pthmat = '/Net/ocean/ddmitry/HYCOM/ARCc/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';

plr=1;
yr=2005;

fmt04 = sprintf('%sARCc0.04_trcr_NLr%2.2i_GreenCoast_%i.mat',pthmat,plr,yr);
fprintf('Loading %s\n',fmt04);
load(fmt04);
TRC4=TRC;

fmt08 = sprintf('%sARCc0.08_trcr_NLr%2.2i_GreenCoast_%i.mat',pthmat,plr,yr);
fprintf('Loading %s\n',fmt08);
load(fmt08);
TRC8=TRC;

ftopo = sprintf('%s/depth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

%fmt08 = sprintf('%sarc08_green_segm.mat',pthmat);
%fprintf('Loading %s\n',fmt08);
%load(fmt08); % SGM

ip = 120;
D4=TRC4.Dpth;
dst4=TRC4.Dist_m(ip,:)*1e-3;
tr4 =TRC4.Tr(ip,:);
mtr4=TRC4.mass_Tr(ip,:);

dst8=TRC8.Dist_m(ip,:)*1e-3; 
tr8 =TRC8.Tr(ip,:);
mtr8=TRC8.mass_Tr(ip,:);

figure(1); clf;

axes('position',[0.08 0.7 0.8 0.25]);
contour(HH,[0 0],'k');
hold on;
contour(HH,[-500 -500],'Color',[0.5 0.5 0.5]);
Ix = TRC8.Indx(ip,:);
[jj,ii] = ind2sub(size(HH),Ix);
plot(ii,jj,'r.');
axis('equal');
set(gca,'xtick',[ ],...
	'ytick',[ ],...
	'xlim',[450 1000],...
	'ylim',[370 950]);


axes('position',[0.08 0.4 0.8 0.25]);
plot(dst4,tr4);
hold on;
plot(dst8,tr8,'r-');
title('Tracer conc., kg/m3, surf. layer, 008-red, 004-blue');


axes('position',[0.08 0.08 0.8 0.25]);
plot(dst4,mtr4);
hold on;
plot(dst8,mtr8,'r-');
title('Depth-intgr (mass/grid cell), kg, 008-red, 004-blue');



