% Plot/analyze fluxes across Greenland contour
% output from vhflx_greenl_cntr008.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2005;
YR2 = 2009;
Cp = 4200; % J/kg K
%Tref= -273.15; % Ref T to calc. H flux
Tref1= -1.8; % Ref T to calc. H flux
Tref2=0;

expt = 110;

rg=9806;  % convert pressure to depth, m
hgg=1e20; 

plr=0; % highlight this interface
btx = 'vhflx_greenl_cntr008.m';

fprintf('Oceanic Heat Flux, Greenland Section, %i-%i\n',YR1,YR2);

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthfig=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/%3.3i/fig_green_xsct/',...
		  regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DXh,DYh]=sub_dx_dy(LON,LAT);
Acell=DXh.*DYh;

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section
nh = length(Hs);
dL = diff(GC.Distance_m);
dL(nh)=dL(nh-1); % over land anyway, dL is unimportant
GC.dL = dL;

% ================================================
% Plot Greenland map and the contour
% ================================================
f_pltgr=0;
if f_pltgr==1
  fn=10;
  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  
  bottom_text(btx,'pwd',1);
end

BOX=[];
VF=[];
HF1=[];
HF2=[];
TM=[];

cc=0;
for iyr=YR1:YR2
  YR=iyr;
  fmat = sprintf('%s%3.3i_GreenlCntr_HVflx_daily_%i.mat',...
		   pthmat,expt,YR);

  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  cc=cc+1;
  if cc==1
    for isc=1:6
      BOX(isc).Name=SCT(isc).Name;
      BOX(isc).T=SCT(isc).Tbox;
      BOX(isc).S=SCT(isc).Sbox;
      
      dmm=SCT(isc).TM;
      SCT(isc).TM=dmm(1:365);
      dmm=SCT(isc).Tbox;
      SCT(isc).Tbox=dmm(1:365);
      dmm=SCT(isc).Sbox;
      SCT(isc).Sbox=dmm(1:365);
      
    end
  else
    for isc=1:6
      dmm=BOX(isc).T;
      BOX(isc).T=[dmm,SCT(isc).Tbox];
      dmm=BOX(isc).S;
      BOX(isc).S=[dmm,SCT(isc).Sbox];
    end
  end
  
    
  Vf=HFLX.VolFlx;
  vf=nansum(Vf,2); % integrate along the contour, m3/s
  Hf=HFLX.HFlx1;
  hf1=nansum(Hf,2);
  Hf=HFLX.HFlx2;
  hf2=nansum(Hf,2);
  
  TM=[TM,HFLX.TM];
  VF=[VF;vf];
  HF1=[HF1;hf1];
  HF2=[HF2;hf2];
  
end

YY=[];
DV=datevec(TM);
for iyr=YR1:YR2
  I=find(DV(:,1)==iyr);
  dmm=[0:length(I)-1]./length(I)+iyr;
  dmm=dmm(:);
  YY=[YY;dmm];
end

mvf=mean(VF);
m10=prctile(VF,10);
m90=prctile(VF,90);

cff=1e-12;
HF1=HF1*cff;
HF2=HF2*cff;

mhf1=mean(HF1);
h110=prctile(HF1,10);
h190=prctile(HF1,90);

mhf2=mean(HF2);
h210=prctile(HF2,10);
h290=prctile(HF2,90);

figure(1); clf;
axes('Position',[0.1 0.7 0.8 0.2]);
plot(YY,VF);
hold
plot([YY(1) YY(end)],[m10 m10],'r--');
plot([YY(1) YY(end)],[m90 m90],'r--');
set(gca,'tickdir','out',...
	'xlim',[YY(1) YY(end)],...
	'xtick',[2000:0.5:2020],...
	'ylim',[1.01*min(VF) 1.01*max(VF)],...
	'xgrid','on',...
	'ygrid','on');

stl=sprintf('0.08 HYCOM-110, GrContour Daily Vol Flux, %8.1f m3/s, %i-%i',...
	    mvf,YR1,YR2);
title(stl);

% Hflux Tref=-1.8
axes('Position',[0.1 0.4 0.8 0.2]);
plot(YY,HF1);
hold
plot([YY(1) YY(end)],[h110 h110],'r--');
plot([YY(1) YY(end)],[h190 h190],'r--');
set(gca,'tickdir','out',...
	'xlim',[YY(1) YY(end)],...
	'xtick',[2000:0.5:2020],...
	'ylim',[1.01*min(HF1) 1.01*max(HF1)],...
	'xgrid','on',...
	'ygrid','on');

stl=sprintf('0.08 HYCOM-110, GrContour Day HeatFlux, %4.1f*1e12 W, Trf=-1.8C, %i-%i',...
	    mhf1,YR1,YR2);
title(stl);

% Hflux Tref=-1.8
axes('Position',[0.1 0.1 0.8 0.2]);
plot(YY,HF2);
hold
plot([YY(1) YY(end)],[h210 h210],'r--');
plot([YY(1) YY(end)],[h290 h290],'r--');
set(gca,'tickdir','out',...
	'xlim',[YY(1) YY(end)],...
	'xtick',[2000:0.5:2020],...
	'ylim',[1.01*min(HF2) 1.01*max(HF2)],...
	'xgrid','on',...
	'ygrid','on');

stl=sprintf('0.08 HYCOM-110, GrContour Day HeatFlux, %4.1f*1e12 W, Trf=0C, %i-%i',...
	    mhf2,YR1,YR2);
title(stl);

btx='plot_vhflx_greenl_cntr008.m';
bottom_text(btx,'pwd',1);

% Plot boxes
CLR=[0.4 0.4 0.4;...
     1 0.5 0;...
     0 0.5 1;...
     0.6 0.2 1;...
     1 0. 0;...
     0 .6 0.3];

figure(2);clf;
axes('Position',[0.08 0.6 0.7 0.3]);
hold on;
for ik=1:6
  clr=CLR(ik,:);
  t=BOX(ik).T;
  plot(YY,t,'Linewidth',2,'Color',clr);
end
set(gca,'tickdir','out',...
	'xlim',[YY(1) YY(end)],...
	'xtick',[2000:0.5:2020],...
	'ylim',[-1.2 4],...
	'xgrid','on',...
	'ygrid','on');

stl=sprintf('0.08 HYCOM-110, T Boxes, %i-%i',...
	    YR1,YR2);
title(stl);


axes('Position',[0.08 0.15 0.7 0.3]);
hold on;
for ik=1:6
  clr=CLR(ik,:);
  s=BOX(ik).S;
  plot(YY,s,'Linewidth',2,'Color',clr);
end
set(gca,'tickdir','out',...
	'xlim',[YY(1) YY(end)],...
	'xtick',[2000:0.5:2020],...
	'ylim',[32.9 34.4],...
	'xgrid','on',...
	'ygrid','on');

stl=sprintf('0.08 HYCOM-110, S Boxes, %i-%i',...
	    YR1,YR2);
title(stl);

axes('Position',[0.8 0.5 0.13 0.4]);
hold on;
for ik=1:6
  clr=CLR(ik,:);
  y1=1-ik/10;
  plot([0.1 0.5],[y1 y1],'-','Linewidth',2,'Color',clr);
  spp=sprintf('%i %s',ik,BOX(ik).Name);
  text(0.6,y1,spp);
end
set(gca,'xlim',[0 0.8],...
	'ylim',[0.3 1],...
	'xtick',[],...
	'ytick',[]);

bottom_text(btx,'pwd',1);
  


  
