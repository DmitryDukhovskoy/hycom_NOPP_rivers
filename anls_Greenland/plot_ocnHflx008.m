% Calculate ocean heat flux to Greenland
% across specified contour -
% isobath around Greenland
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 1993;
YR2 = 2016;

expt = 110;

s_mat = 0; % =0 - extract data, no save; =1 extract & save, =2 - load saved
s_fig = 0;


rg=9806;  % convert pressure to depth, m
hgg=1e20; 

plr=0; % highlight this interface
btx = 'plot_ocnHflx008.m';


regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthfig=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_green_xsct/',...
		  regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

fmat = sprintf('%s%3.3i_Greenl_HVflx_%i.mat',...
		   pthmat,expt,YR1);
load(fmat);

f_pltgr=0;
if f_pltgr==1
  Hg = HH;

  Hg(1:380,:)=nan;
  Hg(1100:end,:)=nan;
  Hg(:,1050:end)=nan;
  Hg(:,1:450)=nan;

  fprintf('Plotting Greenland heat contour ...\n');
  fn=2;
  domname = '0';
  sub_plot_bath(Hg,LON,LAT,fn,domname);
  contour(Hg,[-100:10:-5],'Color',[0.85 0.85 0.85]);
  contour(Hg,[-3500:100:-10],'Color',[0.6 0.6 0.6]);
  contour(Hg,[-4000:500:-50],'Color',[0.25 0.25 0.25]);
  IIs = HFLX(1).GrCntr_II;
  JJs = HFLX(1).GrCntr_JJ;
  x   = HFLX(1).DistCntr*1e-3; % m->km

  plot(IIs,JJs,'b-','Linewidth',2);
  for km=0:500:max(x)
    d=abs(x-km);
    i0=find(d==min(d));
    if km==0
       plot(IIs(i0),JJs(i0),'r.','Markersize',14);
       plot(IIs(i0),JJs(i0),'rd','Markersize',6);
    else
      plot(IIs(i0),JJs(i0),'r.','Markersize',11);
    end
  end

  set(gca,'xlim',[450 1050],...
	  'ylim',[380 1100],...
	  'xtick',[],...
	  'ytick',[]);
  
  title('Contour (~800m) for heat flux calculation');

  bottom_text(btx,'pwd',1);
end

Vflx = [];
Hflx = [];
HF1D = [];
cc = 0;
for yr=YR1:YR2
  fmat = sprintf('%s%3.3i_Greenl_HVflx_%i.mat',...
		   pthmat,expt,yr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  for im=1:12
    cc = cc+1;
    TM(cc) = HFLX(im).TM;
    dmm = HFLX(im).Vol_flux_m3s;
    dmm = nansum(dmm);
    Vflx = [Vflx;dmm];
    dmm = HFLX(im).Hflux_W;
    HF1D(cc,:) = nansum(dmm);  % flux along contour;
    dmm = nansum(nansum(dmm));
    Hflx = [Hflx;dmm];         % total flux
  end
  
end

% Save monthly mean Flux:
%fmm = sprintf('%sARCc08_%3.3i_Greenl_moHflx_%i-%i.mat',pthmat,expt,YR1,YR2);
%save(fmm,'Hflx','TM');


% Annual mean flux:
nyr=length(TM)/12;
mHflx=reshape(Hflx,12,nyr);
sgmHF=std(mHflx,1);
mHflx=mean(mHflx);
yrs=[YR1:YR2];
figure(15); clf;
axes('Position',[0.08 0.5 0.85 0.45]);
hb=bar(yrs,mHflx,0.95);
hold on;
for ir=1:length(yrs)
  lmt1=mHflx(ir)-sgmHF(ir);
  lmt2=mHflx(ir)+sgmHF(ir);
  yr=yrs(ir);
  plot([yr yr],[lmt1 lmt2],'k-','Linewidth',2);
  plot([yr-0.05 yr+0.05],[lmt1 lmt1],'k-','Linewidth',2);
  plot([yr-0.05 yr+0.05],[lmt2 lmt2],'k-','Linewidth',2);
end

set(hb,'Facecolor',[0.8 0.4 0]);
set(gca,'tickdir','out',...
	'xlim',[YR1-0.5 YR2+0.5],...
	'xtick',[1990:2:2020],...
	'xgrid','on',...
	'ygrid','on');

stl=sprintf('ARCc0.08-110, Annual heat flux Greenland contour, %i-%i',YR1,YR2);
title(stl);
bottom_text(btx,'pwd',1,'position',[0.08 0.4 0.7 0.05]);

Hs = HFLX.Hbottom;

% ===============================
% Plot Monthly mean heat fluxes
% along the contour, convert to
% W/m (of the contour)
%
% Calculate monthly climatology:
nm=12;  % now 11 months
[nrc,npp] = size(HF1D);
nyr = nrc/nm; 
A=zeros(nm,npp);
imo=0;
for ii=1:nrc
  dmm=HF1D(ii,:);
  imo=imo+1;
  if (imo==nm+1),
    imo=1;
  end
  A(imo,:)=A(imo,:)+dmm;
end
Hf1d=A/nyr;

dst = HFLX(1).DistCntr; % m
dltX = diff(dst);
dltX(end+1) = dltX(end);
dltX0=dltX;
Ir = find(dltX<100); % segment points repeated - bug in contour subroutine
dltX(Ir) = nan;
dX=nanmean(dltX);

% Filter spatially along the contour
% to get rid of highly oscillatory flux
Wn = 1/40;
[Bf,Af] = butter(9,Wn,'low');

for im=1:nm
  dmm = Hf1d(im,:);
  yy = filtfilt(Bf,Af,dmm)./dX; % W/m
  yy(find(Hs>=0))=nan;
  Hf1d(im,:) = yy;
end

% Plot monthly heat fluxes along the contour:
for im=1:6:12
  cff = 1e-6;
  %im  = 1;
  fn  = im;
  hf1 = Hf1d(im,:)*cff;
  xl1 = 450;
  xl2 = 1050;
  yl1 = 380;
  yl2 = 1100;
  c1  = -100;
  c2  = 100;
  stl = sprintf('Mean HeatFlx, W/m*%3.1d, Mo=%2.2i, %i-%i',cff,im,YR1,YR2);
  sub_plot_hflx1d(HH,LON,LAT,HFLX,hf1,cff,fn,stl,xl1,xl2,yl1,yl2,c1,c2);
  bottom_text(btx,'pwd',1);
  drawnow
  
end

stl = sprintf('%s-%i, Monthly HtFlux, W/m*%3.1d, %2.2i/%i-%i',...
	    regn,expt,cff,YR1,YR2);
dmo = 1;
fnmb=20;
yt1=-700;
yt2=700;
dy=100;
sub_plot_Month_flxsct1d(Hf1d*cff,dmo,dst*1e-3,fnmb,yt1,yt2,dy,stl);
bottom_text(btx,'pwd',1,'position',[0.4 0.3 0.5 0.1]);





% Plot 2D section along the contour:
fnmb = 20;
nlr = 41;
im=6;
dx = HFLX(1).DistCntr*1e-3; % m->km
Hs = HFLX(1).Hbottom;
ZZ = HFLX(im).ZZ;
HFz= HFLX(im).Hflux_W;
% plot(dx,Hs); % plot Bottom profile along contour
% set(gca,'xlim',[0 max(dx)],'xtick',[0:500:max(dx)])
% title('Bottom profile along the heat contour');
% Convert W -> W/m2
ar = ZZ*0;
[Dst,dmb] = meshgrid(dx,[1:nlr]);
cff = 1e-9; 
%      tstr=sprintf('HFlx, W/m2*%3.1d, %s',cff,datestr(dnmb));
Wmp = HFz*cff; %

% Filter along the contour to 
% get rid of high-freq. noise:
Inan=find(isnan(Wmp));
Wmp(Inan)=0;
dZ = diff(ZZ);
dZ=[ZZ(1,:);dZ];
dZ=abs(dZ);
for k=1:nlr
  dh = dZ(k,:);
  gArea = dltX0.*dh;
  gArea(gArea==0)=nan;
  dmm = Wmp(k,:)./gArea; % W/m2
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); % 
%  yy(find(Hs>=0))=nan;
  Wmp(k,:) = yy;
end  

Wmp(Inan)=nan;
% For plotting add surface layer      
ZZp = [ZZ(1,:);ZZ];
ZZp(1,:)=0;
Wmp  = [Wmp(1,:); Wmp]; % for plotting
Dstp = [Dst(1,:);Dst];
c1=-5;
c2=5;
dnmb = datenum(yr,im,15);
tstr=sprintf('HFlx, W/m2*%3.1d, %s',cff,datestr(dnmb));
sub_plot_hflxZcntr(Wmp,ZZp,Dstp,fnmb,tstr,Hs,c1,c2);

bottom_text(btx,'pwd',1);

% --------------------------
% Plot Hovm diagram of heat fluxes
% time vs distance along the contour
cl2 = colormap_red(100);
cl1 = colormap_blue(100);
for ik=1:2;
  cl2(ik,:) = [1 1 1];
  cl1(ik,:) = [1 1 1];
end
cl1 = flipud(cl1);
cmp = [cl1;cl2];
cmp = smooth_colormap(cmp,10);
ncnt = length(cmp);
%c1  = -2; % W/m
%c2  = 2;  
cntr= (c1:(c2-c1)/ncnt:c2);


DV=datevec(TM);
cyr=0;
yrold=0;
cc=0;
for ik=1:length(TM);
  yr=DV(ik,1);
  if (yr~=yrold),
    cc=0;
    yrold=yr;
    nds=length(find(DV(:,1)==yr));
  end
  cc=cc+1;
  YRSd(ik,1)=yr+cc/nds;
end

  
figure(16); clf;
pcolor(dst/1e3,YRSd,HF1D); shading flat;
colormap(cmp); 
caxis([-1e13 1e13]);
set(gca,'tickdir','out',...
	'xlim',[0 9300],...
	'ylim',[1993 2017],...
	'ytick',[1993:2016]);
colorbar;
stl=sprintf('Depth-integr heat flux vs Gr. controu distance, 0.08 HYCOM-%3.3i',expt);
title(stl);

bottom_text(btx,'pwd',1);
