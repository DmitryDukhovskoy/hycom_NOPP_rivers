% Plot Tracer fluxes through Fram Strait
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2000;
YR2 = 2016;

expt = 110;
Zlvl = -50;  % integrate down to 50m and bottom
ntr  = 4;
hgg  = 1e20;
rg   = 9806;

s_fig = 0;

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthfig = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/%3.3i/fig_fram/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx = 'plot_TrFlx_Fram.m';
cBr = sub_BeringTrConc(0);% coefficient for Bering Str., tracer, kg/m3

fprintf('Calculating years: %i - %i\n',YR1,YR2);

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

SCT = sub_set_Fram_xsct(HH,LON,LAT);

clear Fmn Fmin Fmax
cc = 0;
cmn = 0;
for yr=YR1:YR2
  fmat = sprintf('%s%3.3i_TrFlx_Fram_%i.mat',...
		   pthmat,expt,yr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  for imo=1:12
    Tflx  = TFLX(imo).TrFlx_kgs;
    Tflxz = TFLX(imo).TrFlxZ_kgs;
  
    [ntr,pnt] = size(Tflx);
    
% Long-term mean  
    if ~exist('Fmn','var')
      Fmn  = zeros(ntr,pnt);
      Fmin = zeros(ntr,pnt);
      Fmax = zeros(ntr,pnt)-1e9;
      Fmnz = zeros(ntr,pnt);
      Fminz= zeros(ntr,pnt);
      Fmaxz= zeros(ntr,pnt)-1e9;
      Ftot = zeros(ntr,1);
      Ftotz= zeros(ntr,1);
    end

    if yr>2005 % before 2005 - no flux
    for nTr=1:4
      I = find(Tflx<Fmin);
      Fmin(I) = Tflx(I);
      I = find(Tflxz<Fminz);
      Fminz(I) = Tflxz(I);
      I = find(Tflx>Fmax);
      Fmax(I) = Tflx(I);
      I = find(Tflxz>Fmaxz);
      Fmaxz(I) = Tflxz(I);
    end
    end

    cc = cc+1;
    ftot = sum(Tflx,2);
    ftotz= sum(Tflxz,2);
    Ftot(:,cc) = ftot;
    Ftotz(:,cc)= ftotz;
    
    if yr>2005
      cmn=cmn+1;
      Fmn  = Fmn+Tflx;
      Fmnz = Fmnz+Tflxz;
    end
    
    
  end
  
end

Fmn = Fmn/cmn;
Fmnz= Fmnz/cmn;

%nTr = 1;
%fmn = Fmn(nTr,:);
%fmin= Fmin(nTr,:);
%fmax= Fmax(nTr,:);

%fmnz = Fmnz(nTr,:);
%fminz= Fminz(nTr,:);
%fmaxz= Fmaxz(nTr,:);
%
%Wn = 1/10;
%[Bf,Af] = butter(9,Wn,'low');
%dmm = fmnz;
%fmnz = filtfilt(Bf,Af,dmm);
%dmm = fminz;
%fminz = filtfilt(Bf,Af,dmm);
%dmm = fmaxz;
%fmaxz = filtfilt(Bf,Af,dmm);

% =============================================
% Plot Mean/max/min Tracer flux + towards Canada along the section
% Upper 50 m
% =============================================
trnm{1}='Mackeznie';
trnm{2}='E.Euras.R.';
trnm{3}='W.Euras.R.';
trnm{4}='Pacific.W.';
xl1 = 850;
xl2 = 1150;
yl1 = 700;
yl2 = 1150;
IIs = TFLX(1).Indx;
JJs = TFLX(1).Jndx;
z0  = TFLX(1).Z0_m;
npp = length(IIs);
Hs  = SCT.Hbottom;
dst = SCT.Distance_m;
dfd = diff(dst);
dXmn = mean(dfd);


for nTr=1:4
  iTr = nTr+1; % actual tracer #
  fprintf('Plotting Tr %s\n',trnm{nTr});
  Tmn = Fmnz(nTr,:);
  Tmax= Fmaxz(nTr,:);
  Tmin= Fminz(nTr,:);
  dp = 5; % averaging over grid points to plot 1 flux value
  ngrd = 200;
  cff = max(max([abs(Tmax), abs(Tmin)])/ngrd);
  iFmn  = sub_indx_flx(Tmn,dp,cff,IIs,JJs,SCT);
  iFmax  = sub_indx_flx(Tmax,dp,cff,IIs,JJs,SCT);
  iFmin = sub_indx_flx(Tmin,dp,cff,IIs,JJs,SCT);
  
  fn=nTr+1;
  sub_plot_flx_fram(HH,SCT,iFmn,iFmax,iFmin,fn,xl1,xl2,yl1,yl2);

  stl=sprintf('Mean TrFlux %s, 1993-2016, 50m, min(Fmin)=%5.1f kg/s*1m',...
		trnm{nTr},min(Tmin)/dXmn);
  title(stl,'Fontsize',10);

  bottom_text(btx,'pwd',1,'Fotsize',10);

  if s_fig==1
    fgnm=sprintf('%s%3.3i_mean_FlxTr%i_Fram',pthfig,expt,iTr);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r300',fgnm);
  end
end  

% -----------------------------------------------
% Plot time series of integrated fluxes by years
% -----------------------------------------------
xyr = [YR1:1/12:2016.99];
for nTr=1:4
  TrTot = Ftotz(nTr,:);
  iTr = nTr+1;
  
  for iyr=1993:2016
    I=find(xyr<iyr+1 & xyr>=iyr);
    TrY(I)=nanmean(TrTot(I));
  end
  
    
  figure(10+nTr); clf;
  axes('Position',[0.08 0.5 0.85 0.35]);
  plot(xyr,TrTot,'Linewidth',2);
  hold on
  plot(xyr,TrY,'r-');
  set(gca,'xlim',[2000 ceil(xyr(end))],...
	'tickdir','out',...
	'ylim',[1.05*min(TrTot) 1.05*max(TrTot)],...
	'xtick',[1990:1:2020],...
	'Fontsize',16,...
	'xgrid','on',...
	'ygrid','on');


  stl=sprintf('Monthly TrFlux Fram Str %s, 50m, kg/s',...
		trnm{nTr});
  title(stl,'Fontsize',12);

  set(gcf,'Position',[944 508 1596 831]);
  bottom_text(btx,'pwd',1,'position',[0.05 0.3 0.8 0.1],'Fotsize',10);

  if s_fig==1
    fgnm=sprintf('%s%3.3i_mnthTotalFlxTr%i_Fram',pthfig,expt,iTr);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end
end  
