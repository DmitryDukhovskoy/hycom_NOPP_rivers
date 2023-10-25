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

ufld = 'mean'; % mean - annual mean fields, daily - daily U
regn = 'ARCc0.08';
expt = 110;  
YR0 = 2005; % for mean UV - year when particles initialized
YR1 = YR0;
YR2 = 2010;

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_Gr_prt/';

dd1 = datenum(1997,9,15);
dd2 = datenum(2009,08,15);
cc=0;
for iyr=1997:2009
  for im=1:12
    dnmb=datenum(iyr,im,15);
    if dnmb<dd1, continue; end;
    if dnmb>dd2, break; end;
    cc=cc+1;
    TM(cc,1)=dnmb;
  end
end

DV=datevec(TM);
yr1=DV(1,1);
yr2=DV(end,1);
F=load('Fram_FW_NPI.txt'); % km3/yr, monthly

cc=0;
for iyr=yr1:yr2
  I=find(DV(:,1)==iyr);
  dmm=mean(F(I));
  cc=cc+1;
  Fm(cc,1)=dmm;
end

mo1=DV(1,2);
mo2=DV(end,2);
yrr=[yr1+mo1/12:1/12:yr2+mo2/12];

figure(1); clf;
axes('Position',[0.12 0.4 0.86 0.5]);
plot(yrr,F,'k-','Color',[0.7 0.7 0.7],'Linewidth',2);
hold on;
cc=0;
for iyr=yr1:yr2
  cc=cc+1;
  plot([iyr iyr+1],[Fm(cc) Fm(cc)],'k','Linewidth',2);
end

set(gca,'tickdir','out',...
	'xlim',[yrr(1) yrr(end)],...
	'ylim',[-3500 0],...
	'xtick',[yr1:yr2],...
	'ytick',[-3500:500:0],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);

title('FW Flux Fram Strait, NPI, km3/yr');
btx = 'plot_Fram_FW_NPI.m';
bottom_text(btx,'pwd',1,'position',[0.03 0.25 0.4 0.05]);

% Estimate FW flux anomlay during 2007-2008
mnF=mean(Fm(2:10)); % mean 1998-2006
dF = Fm-mnF;
