% Monthly mean FW flux in Fram strait 
% from NPI mooring observations
% http://www.mosj.no/en/climate/ocean/freshwater-flux-fram-strait.html
%
% Norwegian Polar Institute (2018). 
% Freshwater flux in the Fram Strait. 
% Environmental monitoring of Svalbard and Jan
% Mayen (MOSJ). 
% URL: http://www.mosj.no/en/climate/ocean/freshwater-flux-fram-strait.html
%
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

s_fig = 0;
regn = 'ARCc0.08';
expt = 110;  
YR0 = 2005; % for mean UV - year when particles initialized
YR1 = YR0;
YR2 = 2010;

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_Gr_prt/';


shelfM=1; % use time ser with shelf mooring
FWT=sub_read_fram(shelfM);
Fm=FWT.Annual_km3_yr;
TM=FWT.TM;
YRS=FWT.Years;
%i1=find(YRS==2002);
i2=find(YRS==2009);
DV=datevec(TM);
mFm=nanmean(Fm(1:i2));
dFm=-(Fm-mFm); %+ anomaly = + FW flux to Greenland now

Fm(1)=Fm(2);
dFdT=diff(Fm); % Note sign convention is same as FWT, "-" southward
Fm(1)=nan;

yr1=DV(1,1)
yr2=DV(end,1);
mo1=DV(1,2);
mo2=DV(end,2);
yrr=[yr1+mo1/12:1/12:yr2+mo2/12];

FYR = [yr1:yr2];

[PYR,dP,mnP] = sub_prcp_NCEPR;
%
% Greenland anomaly
frv = sprintf('%sGreenland_annAnom.mat',pthmat);
f_annGr=0;
if f_annGr==1
  [TM,Fgr] = sub_read_Greenland_v3; 
  DVg = datevec(TM);
  ngr = length(TM);
  nyr =ngr/12;
  dmm = reshape(Fgr,[12,nyr]);

  Fyr = sum(dmm); % km3/yr
  Ygr = [DVg(1,1):DVg(end,1)];
  ii=find(Ygr==1992);
  GRmn = mean(Fyr(1:ii));
  dGR  = Fyr-GRmn;
  
  save(frv,'GRmn','dGR','Ygr');
else
  fprintf('Loading %s\n',frv);
  load(frv);
end  

figure(1); clf;
axes('Position',[0.1 0.68 0.86 0.25]);
clr=[0. 0.5 0.6];
hold on
%plot([1990 2017],[0 0],'Color',[0.7 0.7 0.7]);
p1=plot(FYR+0.5,dFm,'-','Color',[0. 0.5 0.6],'Linewidth',2);
plot(FYR+0.5,dFm,'.','Color',[0. 0.7 0.8],'MarkerSize',18);
yl1=-1600;
yl2=1600;
nyy=8;
dyy=abs(yl2-yl1)/nyy;

set(gca,'tickdir','out',...
	'xlim',[1990 2017],...
	'ylim',[yl1 yl2],...
	'xtick',[1990:2017],...
	'ytick',[yl1:dyy:yl2],...
	'xgrid','on',...
	'ygrid','on',...
	'ycolor',clr,...
	'Fontsize',14);
title('Anomalies Fram FWFlux km3/yr');

axes('Position',[0.1 0.37 0.86 0.25]);
clr=[0.5 0. 0.4];
hold on
%plot([1990 2017],[0 0],'k-','Color',[0.7 0.7 0.7]);
p2=plot(PYR+0.5,dP,'-','Color',clr,'Linewidth',2);
plot(PYR+0.5,dP,'.','Color',[0.7 0. 0.4],'MarkerSize',18);
yl1=-1200;
yl2=1200;
nyy=8;
dyy=abs(yl2-yl1)/nyy;
set(gca,'tickdir','out',...
	'xlim',[1990 2017],...
	'ylim',[yl1 yl2],...
	'xtick',[1990:2017],...
	'ytick',[yl1:dyy:yl2],...
	'yaxislocation','right',...
	'ycolor',clr,...
	'Fontsize',14);
title('Anomalies Precip km3/yr');


axes('Position',[0.1 0.07 0.86 0.25]);
hold on
plot([1990 2017],[0 0],'k-','Color',[0.7 0.7 0.7]);
p3=plot(Ygr+0.5,dGR,'-','Color',[0.7 0.3 0.1],'Linewidth',2);
plot(Ygr+0.5,dGR,'.','Color',[0.9 0.4 0.],'MarkerSize',18);
%hb=bar(FYR,dFm);
%hc=bar(PYR,dP);

set(gca,'tickdir','out',...
	'xlim',[1990 2017],...
	'ylim',[-200 600],...
	'xtick',[1990:2017],...
	'ytick',[-1600:200:1600],...
	'Fontsize',14);
title('Greenland FWFlux anom, wrt 1957-1992');

spp = legend([p1,p2,p3],{'Fram','Prcp','GrFWF'});
set(spp,'Position',[0.8 0.85 0.1 0.08],...
	'FontSize',14);

%title('Anomalies Components FWFlx, SPNA, km3/yr');

sp1=sprintf('Fram, mean(%i-%i)=%6.1f km3/yr',YRS(1),YRS(i2),abs(mFm));
sp2=sprintf('Prcp NCEPR II, mean(1990-2016)=%6.1f km3/yr',mnP);
sp3=sprintf('GFWF, mean(1957-1992)=%6.1f km3/yr',GRmn);
SPP={sp1;sp2;sp3};

stt = text(2000,-10,SPP);
set(stt,'Fontsize',14);

btx = 'plot_FramFW_precip_v2.m';
bottom_text(btx,'pwd',1,'position',[0.01 0.01 0.4 0.05]);

set(gcf,'Position',[603 31 1482 1311]);
%drawnow

if s_fig==1
  fgnm = sprintf('%sFram_prcip_GrFW_SPNA',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end

% Plot time series of FW fluxes:
figure(3); clf;
axes('Position',[0.1 0.4 0.86 0.5]);
hold on
%plot([1990 2017],[0 0],'k-');
p1=plot(FYR+0.5,abs(Fm),'-','Color',[0. 0.5 0.6],'Linewidth',2.5);
plot(FYR+0.5,abs(Fm),'.','Color',[0. 0.7 0.8],'MarkerSize',18);
plot([2004 2016],[abs(mFm) abs(mFm)],...
     'k--','Linewidth',2,...
     'Color',[0. 0.8 0.9]);
set(gca,'tickdir','out',...
	'xlim',[2003.5 2016.5],...
	'ylim',[1400 3500],...
	'xtick',[1990:2017],...
	'ytick',[0:200:4000],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
stl=sprintf('Fram Strait FWT (km3/yr), mean=%5.1f km3/yr, Sref=34.9',abs(mFm));
title(stl);

bottom_text(btx,'pwd',1,'position',[0.01 0.01 0.4 0.05]);



% Plot time series of FW/dt (km3/yr2) fluxes:
figure(4); clf;
axes('Position',[0.1 0.4 0.86 0.5]);
hold on
%plot([1990 2017],[0 0],'k-');
dFdT(length(FYR))=nan;
dFdT(1)=nan;
p1=plot(FYR+1,-dFdT,'-','Color',[0. 0.5 0.6],'Linewidth',2.5);
plot(FYR+1,-dFdT,'.','Color',[0. 0.7 0.8],'MarkerSize',18);
plot([2000 2017],[0 0],...
     'k--','Linewidth',1);

set(gca,'tickdir','out',...
	'xlim',[2003.5 2016.5],...
	'ylim',[-1500 1500],...
	'xtick',[1990:2017],...
	'ytick',[-2500:250:2000],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
stl=sprintf('Rate of change of Fram FWT (km3/yr2), Sref=34.9',abs(mFm));
title(stl);

bottom_text(btx,'pwd',1,'position',[0.01 0.01 0.4 0.05]);

