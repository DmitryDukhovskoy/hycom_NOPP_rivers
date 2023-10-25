% Plot streamlines prepared
% in streamline_meanUV.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;

regn = 'ARCc0.08';
expt = 110;
yrF  = 2005

pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmat = sprintf('%sNAtlGreenl_particles_%i.mat',pthmat,yrF);

fprintf('Loading %s\n',fmat);
load(fmat);

%YRPLT = [2011,2012,2013,2014,2015];
% Years to average or get monthly field:
YRS = yrF; % more than 1 year - mean over these years
np  = length(YRS);
imo = 13; % >12 - annual mean
if np>1, imo=13; end;

ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
hmsk=HH;
hmsk(hmsk>=0)=0;
hmsk(HH<0)=1;

[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);
xlim1 = 20;
xlim2 = nn-1;
ylim1 = 100;
ylim2 = mm-100;

nd = length(PRTCL.TRACK);
for id=1:nd
  I = PRTCL.TRACK(id).I';
  J = PRTCL.TRACK(id).J';
  IP(id,:)=I;
  JP(id,:)=J;
end

np=size(IP,2);

GC  = sub_greenl_isobath(HH,LON,LAT);
IIg = GC.cntr_Iindx;
JJg = GC.cntr_Jindx;

lcmp = [0 0 0; 1 1 1];
figure(1); clf;
pcolor(hmsk); shading flat;
colormap(lcmp); 
freezeColors;


hold on;
contour(HH,[-5000:1000:-10],'Color',[0.8 0.8 0.8]);

cpp=0;
XP=[];
YP=[];
arrow_dst = 1.5;
dlim = 0.25; % min dist between isolines
for ip=1:2:np
  X0=IP(:,ip);
  Y0=JP(:,ip);
  ncc=length(X0);
  nc2=round(ncc/2);
  dx1=((X0(1)-X0(end)).^2+(Y0(1)-Y0(end)).^2);
  dx2=((X0(1)-X0(nc2)).^2+(Y0(1)-Y0(nc2)).^2);
  if dx1<0.1 & dx2<0.1, continue; end;
% Check if the streamlines are too crowded:
  xmm=IP;
  ymm=JP;
  xmm(:,ip)=nan;
  ymm(:,ip)=nan;
  for ik=1:ncc
    x0=X0(ik);
    y0=Y0(ik);
    D=((xmm-x0).^2+(ymm-y0).^2);
    dmn=nanmin(nanmin(D));
    if dmn<dlim, X0(ik)=nan; Y0(ik)=nan; end;
  end
  
  plot(X0,Y0,'-','Color',[0 0.6 0.8]);
% Plot arowhead   
  x1=X0(end-1);
  x2=X0(end);
  y1=Y0(end-1);
  y2=Y0(end);
  if isempty(XP), XP=X0'*1e9; YP=Y0'*1e9; end;
  D=((X0(end)-XP).^2+(Y0(end)-YP).^2);
  dmn=min(D);
  if dmn>arrow_dst
    cpp=cpp+1;
    XP(cpp,:)=X0';
    YP(cpp,:)=Y0';
    cf=1;
    beta=15;
    v_col = [0 0.4 0.6];
    lwd=1;
    draw_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
  end
  
end
%contour(HH,[-800 -800],'r','Linewidth',2);
plot(IIg,JJg,'r-','Linewidth',2);

axis('equal');
%set(gca,'xlim',[400 1050],...
%	'ylim',[140 1100]);
set(gca,'xlim',[400 1250],...
	'ylim',[130 1200]);
stl=sprintf('ARCc0.08-110, Mean U, %i',YRS);
title(stl);
btx = 'plot_streamlines_meanUV.m';
bottom_text(btx,'pwd',1);



