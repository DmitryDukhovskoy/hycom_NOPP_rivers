% Extract U,V and S to calculate
% across Lom. Ridge
% Do for 0-50m 
% and whole water column
% Collocate U,V with T,S, Tr points
% 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 1;

YR1 = 2003;
YR2 = 2007;

expt = 110;
Sref = 34.8;
Zlvl = -50;  % integrate down to 50m and bottom
ntr  = 4;
hgg  = 1e20;
rg   = 9806;


pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthfig = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/%3.3i/fig_lomonosov/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx = 'plot_Flx_Lomonosov.m';

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
%[DX,DY]=sub_dx_dy(LON,LAT);

SCT = sub_set_xsct(HH,LON,LAT);
Hs  = SCT.Hbottom;
dst = SCT.Distance_m;
dfd = diff(dst);
dXmn = mean(dfd);

Wn = 1/40;
[Bf,Af] = butter(9,Wn,'low');

cyy=0;
clear FW FWd FTR
for yr = 1993:2016
%  if yr==2012, continue; end
  
  fprintf('Reading %i\n',yr);
  fmat = sprintf('%s110_STrFlx_lomonosov_%i.mat',pthmat,yr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  if ~exist('IIs','var');
    IIs = TFLX(1).Indx;
    JJs = TFLX(1).Jndx;
    z0  = TFLX(1).Z0_m;
    npp = length(IIs);
  end
  
% Calculate annual mean fluxes for each tracer
  FTr = zeros(4,npp); % tracer fluxes
  FTrD = zeros(4,npp); % tracer fluxes
  FFw = zeros(1,npp); % FW flx in the upper 50m, m3/s
  FFwD= zeros(1,npp); % over whole water column
  cc = 0;
  for imo = 1:12
    dmm = TFLX(imo).FWflxZ_m3s;
    FFw = FFw+dmm;
    dmm = TFLX(imo).FWflx_m3s;
    FFwD= FFwD+dmm;
    for nTr=1:4
      dmm = TFLX(imo).TrFlxZ_kgs(nTr,:);
      FTr(nTr,:) = FTr(nTr,:)+dmm;
      dmm = TFLX(imo).TrFlx_kgs(nTr,:);
      FTrD(nTr,:) = FTrD(nTr,:)+dmm;
    end
    cc=cc+1;
  end
  FTr = FTr/cc;
  FTrD = FTrD/cc;
  FFw = FFw/cc;
  FFwD= FFwD/cc;
  
  dmm = FFw;
  FFw = filtfilt(Bf,Af,dmm); % m3/s
  FFw(Hs>=0) = nan;
  dmm = FFwD;
  FFwD= filtfilt(Bf,Af,dmm); % m3/s
  FFwD(Hs>=0) = nan;
  for nTr=1:4
    dmm = FTr(nTr,:);
    fmm = filtfilt(Bf,Af,dmm);
    FTr(nTr,:) = fmm;
  end
  
  cyy = cyy+1;
  FW(cyy,:)  = FFw;
  FWd(cyy,:) = FFwD;
  for nTr = 1:4
    FTR(nTr).Flx_Tr(cyy,:) = FTr(nTr,:);
  end
  
  
end

Fmn  = nanmean(FW);
Fmin = min(FW);
Fmax = max(FW);
FDmn = nanmean(FWd);
FDmin = min(FWd);
FDmax = max(FWd);



% =============================================
% Plot FW flux + towards Canada along the section
% =============================================
dp = 10;
ngrd = 150;
cff = max(abs(Fmax)/ngrd);
iFmn  = sub_indx_flx(Fmn,dp,cff,IIs,JJs,SCT);
iFmax = sub_indx_flx(Fmax,dp,cff,IIs,JJs,SCT);
iFmin = sub_indx_flx(Fmin,dp,cff,IIs,JJs,SCT);
%iFDmn = sub_indx_flx(FDmn,dp,cff,IIs,JJs,SCT);
%iFDmax= sub_indx_flx(FDmax,dp,cff,IIs,JJs,SCT);
%iFDmin= sub_indx_flx(FDmin,dp,cff,IIs,JJs,SCT);

fn=1;
sub_plot_flx_lomonos(HH,SCT,iFmn,iFmax,iFmin,fn);

stl = sprintf('Mean FWFlux, 1993-2016, upper50m, Srf=34.8, max(Fmean)=%5.1f m2/s',...
	      max(Fmn)/dXmn);
title(stl,'Fontsize',10);

bottom_text(btx,'pwd',1,'Fotsize',10);

if s_fig==1
  fgnm=sprintf('%s%3.3i_mean_FWFlx_Lomonos',pthfig,expt);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r250',fgnm);
end

% =============================================
% Plot Tracer flux + towards Canada along the section
% =============================================
POS(1,:) = [0.08 0.52 0.46 0.46];
POS(2,:) = [0.48 0.52 0.46 0.46];
POS(3,:) = [0.08 0.03 0.46 0.46];
POS(4,:) = [0.48 0.03 0.46 0.46];

% Here, Tr = 1 - Mackenzie
%          = 2 - East Euras.R.
%          = 3 - West Euras.R.
%          = 4 - Pacific W.

trnm{1}='Mackeznie';
trnm{2}='E.Euras.R.';
trnm{3}='W.Euras.R.';
trnm{4}='Pacific.W.';

for nTr=1:4
  Tr  = FTR(nTr).Flx_Tr;
  Tmn = nanmean(Tr);
  Tmax= max(Tr);
  Tmin= min(Tr);
  dp = 10; % averaging over grid points to plot 1 flux value
  ngrd = 150;
  cff = max(abs(Tmax)/ngrd);
  iFmn  = sub_indx_flx(Tmn,dp,cff,IIs,JJs,SCT);
  iFmax  = sub_indx_flx(Tmax,dp,cff,IIs,JJs,SCT);
  iFmin = sub_indx_flx(Tmin,dp,cff,IIs,JJs,SCT);
  
  fn=nTr+1;
  sub_plot_flx_lomonos(HH,SCT,iFmn,iFmax,iFmin,fn);

  stl=sprintf('Mean TrFlux %s, 1993-2016, 50m, max(FTr mean)=%5.1f kg/s*1m',...
		trnm{nTr},max(Tmn)/dXmn);
  title(stl,'Fontsize',10);

  bottom_text(btx,'pwd',1,'Fotsize',10);

  if s_fig==1
    fgnm=sprintf('%s%3.3i_mean_FlxTr%i_Lomonos',pthfig,expt,nTr);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end
end  


% -----------------------------------------------
% Plot time series of integrated fluxes by years
% -----------------------------------------------
xyr = [1993:2016];
CLR = [0 0.5 1; ...
       1 0.5 0; ...
       0. 1 0.5; ...
       0.2 0.2 0.2];

for nTr=1:4
  Tr  = FTR(nTr).Flx_Tr;
  TrTot = nansum(Tr.*dXmn,2);
  
  figure(10+nTr); clf;
  axes('Position',[0.08 0.42 0.85 0.4]);
%  hb = plot(xyr,TrTot,'Linewidth',2);
  hb = bar(xyr,TrTot,0.96);
  set(hb,'FaceColor',[0.6 0.6 0.6]);
  set(gca,'xlim',[xyr(1) ceil(xyr(end))+0.5],...
	'tickdir','out',...
	'ylim',[1.1*min(TrTot) 1.1*max(TrTot)],...
	'xtick',[xyr(1):1:2020],...
	'xgrid','on',...
        'ygrid','on');

  stl=sprintf('Annual TrFlux Lom.Sect. %s, %i, 50m, kg/s',...
		trnm{nTr});
  title(stl,'Fontsize',12);

  bottom_text(btx,'pwd',1,'Fotsize',10,'Position',[0.02 0.2 0.8 0.1]);

  if s_fig==1
    fgnm=sprintf('%s%3.3i_TotalFlxTr%i_Lomonos',pthfig,expt,nTr);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end
end  
