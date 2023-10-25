% Plot hycom restart fields
% for ARCb0.08 GOFS3.5-like run by NRL
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

hg = 1e20;

%pthbin  = '/Net/mars/ddmitry/hycom/ARCc0.04/restart/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
%fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);
pthbin  = '/Net/kronos/ddmitry/hycom/GOFS3.5/';
%pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.04/restart_022/';  % expt 02.2 - HYCOM-CICE5, GOFS3.5
pthtopo = '/Net/kronos/ddmitry/hycom/GOFS3.5/';
fltopa = sprintf('%sdepth_ARCb0.08_NRL_11.a',pthtopo);
fltopb = sprintf('%sdepth_ARCb0.08_NRL_11.b',pthtopo);

fprintf('Reading topo: %s\n',fltopa);

IDM = 1600;
JDM = 2520;
IJDM = IDM*JDM;

depth_fid=fopen(fltopa,'r');
[h,count]=fread(depth_fid,IJDM,'float32','ieee-be');
HH=reshape(h,IDM,JDM)';
fclose(depth_fid);

I=find(HH>1e20);
HH=-HH;
HH(I)=100;
[mm,nn]= size(HH);
[m,n]= size(HH);
Ib = find(HH<0);

fina = sprintf('%srestart_NRL_arc08_gofs35.a',pthbin);
finb = sprintf('%srestart_NRL_arc08_gofs35.b',pthbin);


%fld = 'u';
fld = 'spd';
%fld = 'dp';
%fld = 'temp';
rlr = 1;
if ~strncmp(fld,'spd',3)
  F = read_hycom_restart(fina,finb,fld,IDM,JDM,'r_layer',rlr);
else
  U = read_hycom_restart(fina,finb,'u',IDM,JDM,'r_layer',rlr);
  V = read_hycom_restart(fina,finb,'v',IDM,JDM,'r_layer',rlr);
  F = sqrt(U.^2+V.^2);
  clear U V
end

  % Check land point:
IL=find(HH<0);
JL=find(F(IL)>1e20);
if ~isempty(JL),
  fprintf('Field huge values are in the ocean points !!!! \n');
  [Jn,In]=ind2sub(size(HH),JL);
end;

I=find(F>1e20);
F(I)=nan;
%F=F./(1000*9.8);

figure(1); clf;
pcolor(F); shading flat;

hold on;
contour(HH,[0 0],'k');
axis('equal');
set(gca,'xlim',[1 nn],...
        'ylim',[1 mm]);


title(sprintf('%s, %s',fina, fld),'Interpreter','none');
colorbar







