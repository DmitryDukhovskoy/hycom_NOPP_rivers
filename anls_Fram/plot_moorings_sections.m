% Calculate volume transport across straits/sections
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;

regn = 'ARCc0.08';
expt = '110';
TV   = 11;  % topo version
%segm = 'BeringS';
 segm = 'FramS';
%segm = 'BarOp';
%segm = 'DavisS';
YR1  = 1993;
YR2  = 2016;

%Zmn = 100; % average over the top Zmn m
rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%s/fig_sections/';

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

% SBE AWI mooring locations
SBE = Fram_moorAWI;
nf=length(SBE);

% lon/lat -> index space

cmp=flipud(colormap_cold);
HHb=HH;
HH(HH>0)=nan;
HH(1000:end,:)=nan;
HH(1:600,:)=nan;
HH(:,1:600)=nan;


f_add=1;
if f_add==1
  sctnm='Fram';
  yr=2006;
  pthfram  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
  fmat_out = sprintf('%sarc08-112_UTSonZ_daily_%s_%4.4i.mat',...
                 pthfram,sctnm,yr);
  fprintf('Loading %s\n',fmat_out);
  load(fmat_out);
  Isct=SCTZ.I;
  Jsct=SCTZ.J;
end  



figure(1); clf;
%axes('Position',[0.05 0.2 0.9 0.72]);
pcolor(HH); shading interp;
hold on;
colormap(cmp);
caxis([-5000 0]);
contour(HH,[-5000:500:-10],'Color',[0.7 0.7 0.7]);



xl1=900;
xl2=1150;
yl1=830;
yl2=980;
axis('equal');
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'Color',[0 0 0]);


lon=LON;
lon(lon>60)=nan;
lon(lon<-160)=nan;
contour(LAT,[60:2:80],'Color',[0.5 0.5 0.5]);
contour(lon,[-24:4:24],'Color',[0.5 0.5 0.5]);
contour(lon,[0 0],'Color',[0.5 0.5 0.5],'Linewidth',1.5);

if f_add==1
% Plot HYCOM section   
  plot(Isct,Jsct,'-','Linewidth',2,'Color',[1 0.7 0]);  
end


for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  XY=[x,y];
  IJ = sub_XY2indx(XY,LON,LAT);
  plot(IJ(1),IJ(2),'r.','Markersize',28);
  nm=SBE(ik).Name;
  if ik==1
    text(IJ(1)+1,IJ(2)+1,nm,'Fontsize',16);
  elseif ik==2
    text(IJ(1)+2,IJ(2)-1,nm,'Fontsize',16);
  elseif ik==3
    text(IJ(1)-3,IJ(2)+2,nm,'Fontsize',16);
  elseif ik==4
    text(IJ(1),IJ(2)+2,nm,'Fontsize',16);
  else
    text(IJ(1)-3,IJ(2)+2,nm,'Fontsize',16);
  end    
end
%
% F17 on shelf mooring - added to SBE
%x17=-8;
%y17=79.029;
%XY=[x17,y17];
%IJ = sub_XY2indx(XY,LON,LAT);
%plot(IJ(1),IJ(2),'r.','Markersize',16);
%nm='F17';
%text(IJ(1),IJ(2),nm,'Fontsize',12);


stl = sprintf('Moorings AWI and NPI, Section: 0.08 HYCOM ');
title(stl);

hcb=colorbar;
set(hcb,'Position',[0.91 0.15 0.01 0.73],...
	'TickLength', 0.02, ...
	'Fontsize',14);

btx='plot_moorings_sections.m';
bottom_text(btx,'pwd',1,'Position',[0.05 0.05 0.4 0.02]);



