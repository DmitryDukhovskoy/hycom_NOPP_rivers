% Read hycom restart fields
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

%pthbin  = '/Net/mars/ddmitry/hycom/ARCc0.04/restart/';
%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);
pthbin  = '/Net/mars/ddmitry/hycom/ARCc0.08/restart/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
fltopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo);

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

fina = sprintf('%srestart_arc008_105a.a',pthbin);
finb = sprintf('%srestart_arc008_105a.b',pthbin);

%fld = 'u';
%fld = 'temp';
fld = 'tracer';
fprintf('Reading restart %s: %s\n',fld,fina);
lp = 1;
nTr = 5; 
F = read_hycom_restart(fina,finb,fld,IDM,JDM,'r_tracer',nTr,'r_layer',lp);
%F = read_hycom_restart(fina,finb,fld,IDM,JDM);

if ndims(F)==3
  A = squeeze(F(lp,:,:));
else
  A = F;
end

A(A>1e20) = nan;
%c1 = min(min(A));
%c2 = max(max(A));

c1=0;
c2=2;

figure(1); clf;
pcolor(A); shading flat;
hold on;
contour(HH,[0 0],'k');
caxis([c1 c2]);
colorbar

stt = sprintf('%s, Fld = %s',fina,fld);
title(stt,'Interpreter','none','Fontsize',12);









