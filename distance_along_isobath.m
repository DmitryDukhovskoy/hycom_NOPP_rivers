% Estimate approximate distance along isobath in the Arctic Ocean
% use [x,y] = getpts 
% to interactively get pnts along the contours:
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps/
startup

format long g
clear all
close

yr   = 2005;
iday = 60; % check winter day to avoid warm surface layer
regn = 'ARCc0.08';
expt = 110;

pfld  = 'temp';
%f_extr = 1;  % =0 - load in extracted depth of Atl. Water, =1 -extract
s_fig  = 0;

sfig=0;

rg = 9806;

pthfig  = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/%3.3i/fig_AtlLayer/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);

% Get topo:
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn] = size(LON);

figure(1); clf;
contour(HH,[0 0],'k');
hold on;
contour(HH,[-500 -500],'b');
contour(HH,[-2000 -2000],'g');

% 500m isobath around Arctic Ocean
IJ= [        1055         959
        1096        1034
        1129        1080
        1143        1163
        1178        1206
        1209        1242
        1172        1289
        1174        1366
        1193        1430
        1199        1510
        1168        1569
        1086        1569
        1028        1551
         973        1567
         911        1608
         876        1641
         833        1680
         776        1680
         739        1639
         718        1631
         694        1670
         634        1685
         564        1695
         503        1670
         443        1621
         433        1544
         462        1482
         517        1438
         552        1387
         587        1325
         601        1266
         630        1253
         667        1217
         706        1191
         772        1137
         823        1096
         885        1065
         942        1008
         985         929
         996         864];

plot(IJ(:,1),IJ(:,2),'r.');
np = length(IJ);
clear X Y D
for ik=1:np
  i0=IJ(ik,1);
  j0=IJ(ik,2);
  x0=LON(j0,i0);
  y0=LAT(j0,i0);
  X(ik,1)=x0;
  Y(ik,1)=y0;
end

for k=1:np-1
  x1=X(k);
  y1=Y(k);
  x2=X(k+1);
  y2=Y(k+1);
  D(k,1)=distance_spheric_coord(y1,x1,y2,x2); % m
end
D(k+1,1)=D(k);
sD = sum(D); % total distance, m

uatl=0.02; % average speed of Atl. water 2 cm/s
T = sD/uatl;
T = T/(3600*24); % days to travel along the isobath

