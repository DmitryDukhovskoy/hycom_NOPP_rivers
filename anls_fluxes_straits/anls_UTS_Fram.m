% Plot Fram xsection U,T,S
% calc EOF
% analyze variability
% Calc fluxes
% Extracted 2D arrays in extract_daily_UTS_Fram.m
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt = 112;
TV   = 11;  % topo version
sctnm='Fram';

s_mat=1;
YR1=2005;
YR2=YR1;

rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
btx='anls_UTS_Fram.m';


fmat = sprintf('%sarc08-%3.3i_UTS_daily_%s_%4.4i.mat',...
	       pthmat,expt,sctnm,YR1);

if s_mat==0,
  fprintf('No mat files saved ... \n');
else
  fprintf('Data will be saved -> %s\n',fmat);
end


ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

fprintf('Loading %s\n',fmat);
load(fmat);

TM=SCT.TM;
nrc=length(TM);
IS=SCT.I;
JS=SCT.J;
X=SCT.X(1:end-1);
Y=SCT.Y(1:end-1);
%X=SCT.SGM_Lon;
%Y=SCT.SGM_Lat;


f_sct=0;
if f_sct==1
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:1000:-10],'Color',[0.7 0.7 0.7]);
  contour(LAT,[79 79],'b');
  axis('equal');
  set(gca,'xlim',[800 1200],...
	  'ylim',[800 1100]);
  plot(IS,JS,'r.');
  title('Section, 79N');
  
end


f_chck=0;
if f_chck==1
  nfg=2;
  it=20;
  dnmb=TM(it);
  dv=datevec(dnmb);
  Hb=SCT.SGM_Hb;
  T=squeeze(SCT.Temp(it,:,:));
  S=squeeze(SCT.Saln(it,:,:));
  V=squeeze(SCT.Unrm(it,:,:));
  dH=squeeze(SCT.dH_thkn(it,:,:));
  [ma,na]=size(dH);
  ZZ=zeros(1,210);
  a=-cumsum(dH,1);
  ZZ=[ZZ;a];
  
  stl=sprintf('arc08-%3.3i %s, T, %4.4i/%2.2i/%2.2i',...
	      expt,sctnm,dv(1:3));
  cntr=[-2:0.5:10];
  cc0=0;
  sub_plot_Txsct(nfg,Hb,T,ZZ,X,dnmb,stl,cntr,cc0);
  bottom_text(btx,'pwd',1);
  
  
  stl=sprintf('arc08-%3.3i %s, S, %4.4i/%2.2i/%2.2i',...
	      expt,sctnm,dv(1:3));
  cntr=[31:0.5:35];
  cc0=34.8;
  nfg=3;
  sub_plot_Sxsct(nfg,Hb,S,ZZ,X,dnmb,stl,cntr,cc0);
  bottom_text(btx,'pwd',1);
  
% V component
% Project on a 79 Lat to avoid "steps" discontinuities
% due to zig-zaging segments
  dI=diff(SCT.I);
  I0=find(dI>0);
  V0=V(:,I0);
  Hb0=Hb(I0);
  ZZ0=ZZ(:,I0);
  X0=X(I0);
  
  stl=sprintf('arc08-%3.3i %s, V, %4.4i/%2.2i/%2.2i',...
	      expt,sctnm,dv(1:3));
  cntr=[-0.5:0.1:0.5];
  cc0=0;
  nfg=3;
  sub_plot_Vxsct(nfg,Hb0,V0,ZZ0,X0,dnmb,stl,cntr,cc0);
  bottom_text(btx,'pwd',1);


end




