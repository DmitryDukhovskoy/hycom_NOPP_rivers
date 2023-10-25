% Extract U,V and S to calculate
% across Lom. Ridge
% Do for 0-50m 
% and whole water column
% Collocate U,V with T,S, Tr points
% For selected years
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 1;

pyr = [2007,2012,2013,2014];
nyr = length(pyr);

expt = 110;
Sref = 34.8;
Zlvl = -50;  % integrate down to 50m and bottom
ntr  = 4;
hgg  = 1e20;
rg   = 9806;


pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthfig = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/%3.3i/fig_lomonosov/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx = 'plot_Flx_slctYrs_Lomonosov.m';


CLR=[0   0   1;
     0  1   0.5;
     1  0   0;
     0 0.5   1;
     1  0  1;
     0 0.5  0.5;
     0.5 0.5 0;
     0  1    0;
     0.5 0.5 0.5;
     1 0.5  0;
     1  0 0.5;
     0.5 0 1;
     0. 0  0;
     0  0.5  0;
     0.2 1 0.5];


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
clear FW FWd FTR YRS
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
    end
    cc=cc+1;
  end
  FTr = FTr/cc;
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
  
  YRS(cyy)=yr;
end

Fmn  = nanmean(FW);
Fmin = min(FW);
Fmax = max(FW);
FDmn = nanmean(FWd);
FDmin = min(FWd);
FDmax = max(FWd);


% =============================================
% Plot Tracer flux + towards Canada along the section
% =============================================
trnm{1}='Mackeznie';
trnm{2}='E.Euras.R.';
trnm{3}='W.Euras.R.';
trnm{4}='Pacific.W.';

hmsk=HH;
hmsk(HH<0)=nan;
cll = colormap([0.65 0.65 0.65]);

xl1 = 600;
xl2 = 1200;
yl1 = 1000;
yl2 = 1850;

for nTr=1:4
  
  figure(nTr); clf;

% Plot Mean FWflux , upper 50m
  pcolor(hmsk); shading flat;
  colormap(cll);
  freezeColors;
%contour(HH,[0 0],'k');
  hold on;

  contour(HH,[-5000:500:-10],'Color',[0.75 0.75 0.75]);
  Ept=SCT.End_Pnts;
  for ig=1:2
    plot([Ept(ig,1) Ept(ig+1,1)],[Ept(ig,2) Ept(ig+1,2)],...
       '-','Color',[0.2 0.2 0.2],...
       'Linewidth',1.6);
  end
  
  axis('equal');
  set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[],...
	'ytick',[]);

  Tr  = FTR(nTr).Flx_Tr;
  Tmax = max(max(abs(Tr)));
  dp = 10; % averaging over grid points to plot 1 flux value
  ngrd = 180;
  cff = max(abs(Tmax)/ngrd);
  clear cyr
  
  tmax = 0;
  for iir = 1:nyr
    ipl = find(YRS==pyr(iir));
    Tmn = Tr(ipl,:);
    iFmn  = sub_indx_flx(Tmn,dp,cff,IIs,JJs,SCT);
    
    tmax=max([tmax,max(abs(Tmn))]);
    clr=CLR(iir,:);
  
%  sub_plot_flx_lomonos(HH,SCT,iFmn,iFmax,iFmin,fn);
    plot(iFmn(:,1),iFmn(:,2),'-','Color',clr,'Linewidth',2);
    cyr{iir}=sprintf('%i',pyr(iir));
  
  end

  stl=sprintf('TrFlux %s, 50m, max(abs(allFTr))=%5.1f kg/s*1m',...
		trnm{nTr},tmax/dXmn);
  title(stl,'Fontsize',10);

  axes('Position',[0.05 0.5 0.05 0.08]);
  hold on;
  for iir = 1:nyr
    clr=CLR(iir,:);
    plot([0 1],[iir iir],'-','Color',clr);
    text(1.2,iir,cyr(iir));
  end
  set(gca,'visible','off');
  
  
  bottom_text(btx,'pwd',1,'Fotsize',10);

  if s_fig==1
    fgnm=sprintf('%s%3.3i_FlxTr%i_Yrs%i_lomonos',pthfig,expt,nTr,pyr(1));
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end
end  
     
