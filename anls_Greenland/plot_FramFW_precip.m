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

%dd1 = datenum(1997,9,15);
%dd2 = datenum(2009,08,15);
%cc=0;
%for iyr=1997:2009
%  for im=1:12
%    dnmb=datenum(iyr,im,15);
%    if dnmb<dd1, continue; end;
%    if dnmb>dd2, break; end;
%    cc=cc+1;
%    TM(cc,1)=dnmb;
%  end
%end
%
%DV=datevec(TM);
%yr1=DV(1,1);
%yr2=DV(end,1);

shelfM=1; % use time ser with shelf mooring
FWT=sub_read_fram(shelfM);
Fm=FWT.Annual_km3_yr;
TM=FWT.TM;
YRS=FWT.Years;
%i1=find(YRS==2002);
i2=find(YRS==2009);
DV=datevec(TM);

%F=load('Fram_FW_NPI.txt'); % km3/yr, monthly

%cc=0;
%for iyr=yr1:yr2
%  I=find(DV(:,1)==iyr);
%  dmm=mean(F(I));
%  cc=cc+1;
%  Fm(cc,1)=dmm;
%end

% Reference mean: 1998-2002
% and 2003-2009
%mF1=nanmean(Fm(1:i1));
%mF2=nanmean(Fm(i1+1:i2));
%dFm(1:i1)=Fm(1:i1)-mF1;
%dFm(i1+1:end)=Fm(i1+1:end)-mF2;
%dFm=-dFm; % + anomaly = + FW flux to Greenland now
mFm=nanmean(Fm(1:i2));
dFm=-(Fm-mFm); %+ anomaly = + FW flux to Greenland now

% Rate of change of FWT (km3/yr2):
Fm(1)=Fm(2);
dFdT=-diff(Fm); % Note sign convention is same as FWT, changed to + southward
Fm(1)=nan;


yr1=DV(1,1)
yr2=DV(end,1);
mo1=DV(1,2);
mo2=DV(end,2);
yrr=[yr1+mo1/12:1/12:yr2+mo2/12];

FYR = [yr1:yr2];

dmm=dFm;
dmm(isnan(dmm))=0;
%smFm=cumsum(dmm); % cumulative FWT anomaly, Fram
smFm=cumsum(dFdT); % Note sign, where + direction is 
smFm=[nan;smFm];
dFdT(length(Fm))=nan; % make same length
dFdT(1)=nan;

%
% Precipitation
%
[PYR,dP,mnP] = sub_prcp_NCEPR;
smPr=cumsum(dP);
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
smGR=cumsum(dGR);


% =================================
% Beaufort Gyre observations
% BGOS
% =================================
arBG = 1.023e12; % area of BG region in OBGS, m2
YRB=[2003:2017];
% FWC x1e3 km3
fwc=[16.9; 17.2; 18.2; 18.7; 19.8;,...
     21.8; 21.8; 21.9; 21.8; 22.1; ...
     20.6; 21.4; 22.4; 23.4; 23.5]*1e3;

mfwc = fwc(1);
dfw = diff(fwc); % change, km3/yr
dBG = [0;dfw];
smBG=fwc-mfwc;  % FWC anomaly BG, integrated in time, km3



figure(1); clf;
axes('Position',[0.1 0.4 0.86 0.5]);
hold on
%plot([1990 2017],[0 0],'k-');
p1=plot(FYR+0.5,dFm,'-','Color',[0. 0.5 0.6],'Linewidth',2.5);
plot(FYR+0.5,dFm,'.','Color',[0. 0.7 0.8],'MarkerSize',18);

p2=plot(PYR+0.5,dP,'-','Color',[0.5 0. 0.4],'Linewidth',2.5);
plot(PYR+0.5,dP,'.','Color',[0.7 0. 0.4],'MarkerSize',18);

p3=plot(Ygr+0.5,dGR,'-','Color',[0.7 0.3 0.1],'Linewidth',2.5);
plot(Ygr+0.5,dGR,'.','Color',[0.9 0.4 0.],'MarkerSize',18);
%hb=bar(FYR,dFm);
%hc=bar(PYR,dP);

p4=plot(YRB+0.5,dBG,'-','Color',[0.7 0.3 0.7],'Linewidth',2.5);
plot(YRB+0.5,dBG,'.','Color',[0.9 0.5 .9],'Markersize',18);

set(gca,'tickdir','out',...
	'xlim',[1990 2018],...
	'ylim',[-1550 2050],...
	'xtick',[1990:2017],...
	'ytick',[-1500:250:2000],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);

spp = legend([p1,p2,p3,p4],{'Fram','Prcp','GrFWF','BG'});
set(spp,'Position',[0.7 0.15 0.11 0.18],...
	'FontSize',14);

title('Rate of change Components FWFlx, SPNA, km3/yr');

btx = 'plot_FramFW_precip.m';
bottom_text(btx,'pwd',1,'position',[0.03 0.2 0.4 0.05]);

set(gcf,'Position',[856 483 1700 852]);
if s_fig==1
  fgnm = sprintf('%sFram_prcip_GrFW_SPNA',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end

% ---------------------------------
% Plot Time-integrated FW anomalies
% ---------------------------------
figure(2); clf;
axes('Position',[0.1 0.4 0.86 0.5]);
hold on
%plot([1990 2017],[0 0],'k-');
p1=plot(FYR+0.5,smFm,'-','Color',[0. 0.5 0.6],'Linewidth',2.5);
plot(FYR+0.5,smFm,'.','Color',[0. 0.7 0.8],'MarkerSize',18);

p2=plot(PYR+0.5,smPr,'-','Color',[0.5 0. 0.4],'Linewidth',2.5);
plot(PYR+0.5,smPr,'.','Color',[0.7 0. 0.4],'MarkerSize',18);

p3=plot(Ygr+0.5,smGR,'-','Color',[0.7 0.3 0.1],'Linewidth',2.5);
plot(Ygr+0.5,smGR,'.','Color',[0.9 0.4 0.],'MarkerSize',18);
%hb=bar(FYR,dFm);
%hc=bar(PYR,dP);

p4=plot(YRB+0.5,smBG,'-','Color',[0.7 0.3 0.7],'Linewidth',2.5);
plot(YRB+0.5,smBG,'.','Color',[0.9 0.5 .9],'Markersize',18);

set(gca,'tickdir','out',...
	'xlim',[1990 2018],...
	'ylim',[-3000 8000],...
	'xtick',[1990:2017],...
	'ytick',[-4000:1000:8000],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);

spp = legend([p1,p2,p3,p4],{'Fram','Prcp','GrFWF','BG'});
set(spp,'Position',[0.7 0.15 0.11 0.18],...
	'FontSize',14);

title('TimeIntgr Anomalies Components FWFlx, SPNA, km3');

sp1=sprintf('Fram, mean(%i-%i)=%6.1f km3/yr',YRS(1),YRS(i2),abs(mFm));
sp2=sprintf('Prcp NCEPR II, mean(1990-2009)=%6.1f km3/yr',mnP);
sp3=sprintf('GFWF, mean(1957-1992)=%6.1f km3/yr',GRmn);
SPP={sp1;sp2;sp3};

stt = text(1990,-4500,SPP);
set(stt,'Fontsize',14);

set(gcf,'Position',[856 483 1700 852]);

btx = 'plot_FramFW_precip.m';
bottom_text(btx,'pwd',1,'position',[0.03 0.2 0.4 0.05]);

if s_fig==1
  fgnm = sprintf('%sFram_tintgr_prcip_GrFW_SPNA',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end

