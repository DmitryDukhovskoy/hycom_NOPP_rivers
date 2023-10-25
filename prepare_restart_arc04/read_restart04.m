% Read & plot hycom restart fields
% for ARCc0.04
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

%pthbin  = '/Net/mars/ddmitry/hycom/ARCc0.04/restart/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
%fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);
pthbin  = '/Net/mars/ddmitry/hycom/ARCc0.04/restart/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);

fprintf('Reading topo: %s\n',fltopo);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn]= size(HH);
[m,n]= size(HH);
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

fina = sprintf('%srestart_105a.a',pthbin);
finb = sprintf('%srestart_105a.b',pthbin);

%fld = 'u';
%fld = 'temp';
%fld = 'dpmxl';
%F = read_hycom_restart(fina,finb,fld,IDM,JDM);

%fld = 'u'; % make nTr=1 to read u, v, etc.
%fld = 'v';
%fld = 'dp';
%fld = 'temp';
%fld = 'saln';
fld = 'tracer';
lp = 1;
nTr = 5; 
fprintf('Reading restart %s, #tr %i, lev %i: %s\n',fld,nTr,lp,fina);
F = read_hycom_restart(fina,finb,fld,IDM,JDM,'r_tracer',nTr,'r_layer',lp);

if ndims(F)==3
  A = squeeze(F(lp,:,:));
else
  A = F;
end

A(A>1e20) = nan;
%c1 = min(min(A));
%c2 = max(max(A));

c1=0;
c2=1;

fprintf('Plotting ...\n');
figure(1); clf;
pcolor(A); shading flat;
hold on;
contour(HH,[0 0],'k');
caxis([c1 c2]);
colorbar

stt = sprintf('%s, Fld = %s',fina,fld);
title(stt,'Interpreter','none','Fontsize',12);









