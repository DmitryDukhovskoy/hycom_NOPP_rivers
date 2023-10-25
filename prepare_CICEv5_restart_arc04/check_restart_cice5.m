% CICE5 is mapd from GLBc0.04 GOFS3.5 CICE5
% using D. Hebert code
% on Gordon
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear


pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthdat = '/Net/kronos/ddmitry/hycom/ARCc0.04/cice/';

frst = sprintf('%scice.restart.2017010109_interp.nc',pthdat);


ftopo = sprintf('%s/depth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

%
%uu = nc_varget(frst,'uvel');

% Plot conc:
ncat = 5;
dmm = nc_varget(frst,'aicen');
aice = squeeze(sum(dmm,1));
clear dmm

figure(1); clf;
contour(HH,[0 0],'k');
hold on;

aice(aice==0)=nan;
pcolor(aice); shading flat;
caxis([0.2 1]);

axis('equal');
set(gca,'xlim',[0 nn],...
        'ylim',[0 mm]);

clb=colorbar;
stl = sprintf('%s',frst);
title(frst,'Interpreter','none');
btx = 'check_restart_cice5.m';

bottom_text(btx,'pwd',1);


% Ice thickness:
% volume per 1 m2:
dmm = nc_varget(frst,'vicen');
hice = squeeze(sum(dmm,1));
clear dmm

%cmp = colormap('jet');

figure(2); clf;
contour(HH,[0 0],'k');
hold on;

hice(hice==0)=nan;
pcolor(hice); shading flat;
caxis([1 5]);

axis('equal');
set(gca,'xlim',[0 nn],...
        'ylim',[0 mm]);

clb=colorbar;
stl = sprintf('Hice, %s',frst(40:end));
title(frst,'Interpreter','none');

bottom_text(btx,'pwd',1);




