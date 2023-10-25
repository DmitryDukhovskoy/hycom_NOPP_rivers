% Plot tracer flux from Greenland
% across specified contour -
% isobath around Greenland
% see trFlux_contour_greenl.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=2000;
YR2=2016;

f_map=0; % fluxes color-coded along contour, map, UV

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthfig=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_green_xsct/',...
		  regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';
btx = 'plot_TrFlx_contour.m';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section

% ================================================
% Plot Greenland map and the contour
% ================================================
f_pltgr=0;
if f_pltgr==1
  fn=10;
  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  bottom_text(btx,'pwd',1);
end



TflxD = [];
Tflxz = [];
HF1D = [];
cc = 0;
for yr=YR1:YR2
  fmat = sprintf('%s%3.3i_Greenl_TrFlx_%i.mat',...
		   pthmat,expt,yr);

  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  for im=1:12
    cc = cc+1;
    TM(cc) = FWFLX(im).TM;
    ZZ = FWFLX(im).ZZ;
    iz = find(ZZ<-50);
    dmm = FWFLX(im).TrFlux_kg_s;
    dall= nansum(dmm);
    dmm(iz)=0;
    dmz = nansum(dmm); 
   
    TflxD = [TflxD;dall]; % Tr flux whole depth along cntr
    Tflxz = [Tflxz;dmz];  % Tr flux upper 50m, cntr
    
  end
  
end

%A=B
dst = FWFLX(1).DistCntr; % m
dltX = diff(dst);
dltX(end+1) = dltX(end);
dltX0=dltX;
Ir = find(dltX<100); % segment points repeated - bug in contour subroutine
dltX(Ir) = nan;
dX=nanmean(dltX);

% Mean Tr Flux, kg/s m
aa=nanmean(TflxD)./dltX;
aa(aa==0)=nan;
fprintf('Mean tr flux, kg/s m = %8.6d\n',nanmean(aa));

% Plot Fluxes: whole depth
% kg/s per 1 m
% Filter spatially along the contour
% to get rid of highly oscillatory flux
Wn = 1/40;
[Bf,Af] = butter(9,Wn,'low');
[Bfh,Afh] = butter(9,1/20,'low');

% Convert all data to kg/s per 1 m
TflxD=TflxD./dX;
Tflxz=Tflxz./dX;

% Overall mean:
mTflxD = nanmean(TflxD); % whole depth
dmm    = mTflxD; % kg/s -> kg/s per 1m
yy     = filtfilt(Bf,Af,dmm); % W/m
yy(find(Hs>=0))=nan;
mTflxD=yy;

mTflxz = nanmean(Tflxz); % whole depth
dmm    = mTflxz; % kg/s -> kg/s per 1m
yy     = filtfilt(Bf,Af,dmm); % W/m
yy(find(Hs>=0))=nan;
mTflxz=yy;



nrc=length(TM);
for ik=1:nrc
  dmm    = TflxD(ik,:); % kg/s -> kg/s per 1m
  yy     = filtfilt(Bf,Af,dmm); % W/m
  yy(Hs>=0)=nan;
  fTflxD(ik,:)=yy;
  
  dmm    = Tflxz(ik,:); % kg/s -> kg/s per 1m
  yy     = filtfilt(Bf,Af,dmm); % W/m
  yy(Hs>=0)=nan;
  fTflxz(ik,:)=yy;
end
p10 = prctile(fTflxD,10,1);
p10(isnan(p10))=0;
yy  = filtfilt(Bfh,Afh,p10);
p10 = yy;
p10(Hs>=0)=nan;
p90 = prctile(fTflxD,90,1);
p90(isnan(p90))=0;
yy  = filtfilt(Bfh,Afh,p90);
p90 = yy;
p90(Hs>=0)=nan;

p10z = prctile(fTflxz,10,1);
p10z(isnan(p10z))=0;
yy  = filtfilt(Bfh,Afh,p10z);
p10z = yy;
p10z(Hs>=0)=nan;
p90z = prctile(fTflxz,90,1);
p90z(isnan(p90z))=0;
yy  = filtfilt(Bfh,Afh,p90z);
p90z = yy;
p90z(Hs>=0)=nan;

% Fraction of basin-shelf flux 
% in Lab Sea:
TFlx = nansum(mTflxD);
TLb  = nansum(mTflxD(610:1133));
rF=TLb/TFlx;

% Whole depth, mean fluxes
xsct = dst*1e-3; % km
yl1  = 1.02*min(p10);
yl2  = 1.02*max(p90);
yt1  = -30;
yt2  = 10;
dy   = 5;
fnmb = 1;
sub_plot_flxsct1d(mTflxD, xsct,fnmb,yl1,yl2,yt1, yt2, dy);
plot(xsct,p10,'Color',[0.6 0.6 0.6],'linewidth',1.8);
plot(xsct,p90,'Color',[0.6 0.6 0.6],'linewidth',1.8);
stl = sprintf('%s-%i, Mean 10/90prct TrFlux Btm-Srf, kg/s per 1m, %i-%i',...
	      regn,expt,YR1,YR2);
title(stl);
bottom_text(btx,'pwd',1,'position',[0.02 0.3 0.4 0.05]);

% Plot Flux in the upper 50m:
xsct = dst*1e-3; % km
yl1  = 1.02*min(p10z);
yl2  = 1.02*max(p90z);
yt1  = -15;
yt2  = 15;
dy   = 1;
fnmb = 3;
sub_plot_flxsct1d(mTflxz, xsct,fnmb,yl1,yl2,yt1, yt2, dy);
plot(xsct,p10z,'Color',[0.6 0.6 0.6],'linewidth',1.6);
plot(xsct,p90z,'Color',[0.6 0.6 0.6],'linewidth',1.6);
stl = sprintf('%s-%i, Mean 10/90prct TrFlux 50m, kg/s per 1m, %i-%i',...
	      regn,expt,YR1,YR2);
title(stl);
bottom_text(btx,'pwd',1,'position',[0.02 0.3 0.4 0.05]);

%keyboard

% Analyze seasonality:
npp = length(dltX);
cc=0;
TrD_mo=zeros(12,npp);
Trz_mo=zeros(12,npp);
for ik=1:nrc
  im=mod(ik,12);
  if im==0, im=12; end;
  dmm = TflxD(im,:);
  TrD_mo(im,:)=TrD_mo(im,:)+dmm;
  dmm = Tflxz(im,:);
  Trz_mo(im,:)=Trz_mo(im,:)+dmm;
end
cc=nrc/12;
TrD_mo=TrD_mo/cc;
Trz_mo=Trz_mo/cc;

for im=1:12
  dmm = TrD_mo(im,:);
  dmm(isnan(dmm))=0;
  yy  = filtfilt(Bf,Af,dmm);
  yy(Hs>=0)=nan;
  TrD_mo(im,:)=yy;

  dmm = Trz_mo(im,:);
  dmm(isnan(dmm))=0;
  yy  = filtfilt(Bf,Af,dmm);
  yy(Hs>=0)=nan;
  Trz_mo(im,:)=yy;
end

% Plot monthly fluxes:
dmo=1;
fnmb=4;
yt1=-30;
yt2=20;
dy=5;
stl = sprintf('%s-%i, Monthly TrFlux WhDpth, kg/s per 1m, %i-%i',...
	      regn,expt,YR1,YR2);
sub_plot_Month_flxsct1d(TrD_mo,dmo,xsct,fnmb,yt1,yt2,dy,stl);
bottom_text(btx,'pwd',1,'position',[0.4 0.3 0.5 0.1]);

fnmb=5;
yt1=-20;
yt2=20;
dy=2;
stl = sprintf('%s-%i, Monthly TrFlux 50m, kg/s per 1m, %i-%i',...
	      regn,expt,YR1,YR2);
sub_plot_Month_flxsct1d(Trz_mo,dmo,xsct,fnmb,yt1,yt2,dy,stl);

bottom_text(btx,'pwd',1,'position',[0.4 0.3 0.5 0.1]);


% --------------------------
% Annual Mean
% --------------------------
nyr = nrc/12;
cff=1e-9;
A=nansum(TflxD,2).*dX*3600*24*365*1e-3; % kg/(s*m) * m*s=> ton per year
dmm=reshape(A,[12,nyr]);
Mdy=mean(dmm)*cff;

A=nansum(Tflxz,2).*dX*3600*24*365*1e-3; % kg/(s*m) * m*s=> ton per year
dmm=reshape(A,[12,nyr]);
Mzy=mean(dmm)*cff;


yrx = [YR1:YR2];
figure(6); clf;
axes('Position',[0.09 0.55 0.85 0.35]);
bar(yrx,Mdy);
set(gca,'tickdir','out',...
	'xlim',[YR1-0.5 YR2+0.5]);
stl=sprintf('%s-%i, Annual TrFlux WhDpth, ton/yr*%3.1d',regn,expt,cff);
title(stl);

axes('Position',[0.09 0.08 0.85 0.35]);
bar(yrx,Mzy);
set(gca,'tickdir','out',...
	'xlim',[YR1-0.5 YR2+0.5]);
stl=sprintf('%s-%i, Annual TrFlux 50m, ton/yr*%3.1d',regn,expt,cff);
title(stl);
bottom_text(btx,'pwd',1);


if f_map==1
  fnmb=8;
  cff = 1;
  %im  = 1;
  fn  = im;
  hf1 = mTflxD;
  xl1 = 450;
  xl2 = 1050;
  yl1 = 380;
  yl2 = 1100;
  c1  = -10;
  c2  = 10;
  stl = sprintf('%s-%3.3i,Mean TrFlx WhDpth, kg/s m, %i-%i',regn,expt,YR1,YR2);
  sub_plot_trflx_map(HH,LON,LAT,GC,hf1,...
		     cff,fnmb,stl,xl1,xl2,yl1,yl2,c1,c2);
  bottom_text(btx,'pwd',1);
%  drawnow
  
end



