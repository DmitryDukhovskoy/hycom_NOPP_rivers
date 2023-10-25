% Fit some analytical function to
% Bamber river data 
% km3/mo
%
% Note logi function works well
% g(x)=exp(x)/(1+exp(x))
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup

close all
clear

s_fig=0;

hg=2^100; 

pthG   = '/nexsan/people/ddmitry/Net_ocean/arctic_AOregimes/data/GreenlandRunoffv3/';
pthfig = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/fig_Greenland_river/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat/';
btmtxt = 'anls_GreenlandRunoff_v3.m';

YearGr=2004;
mday=[31;28;31;30;31;30;31;31;30;31;30;31];
if mod(YearGr,4)==0,
  mday(2)=29;
end


% New Bamber's data 
% for comparison:
fmat = sprintf('%sGreenlFWFv3.mat',pthmat);
fprintf('Loading %s\n',fmat);
load(fmat)

imn = 36;         % # of yrs, calculated reference mean: 1958-1993
YRn = GrF.Years;
FWFn = GrF.FWFtotal_km3yr;
FWF  = GrF.TotalFWF_km3mo;
fmn93 = mean(FWFn(1:imn)); % mean flux 1958-1993
%% Means  for components:
dmm       = GrF.Tundra_Runoff_km3yr;
tundra_mn = mean(dmm(1:imn)); % ref mean tundra
dmm       = GrF.GrMeltWater_km3yr;
meltw_mn  = mean(dmm(1:imn)); % ref mean tundra
dmm       = GrF.IceDisch_km3yr;
idsch_mn  = mean(dmm(1:imn)); % ref mean tundra


% Anomlay 1993 - 2016
yr1 = YRn(1);
yr2 = YRn(end);
  
% Group by months/years
% km3/mo
nrc = length(FWF);
nyr = nrc/12;
FTMn = reshape(FWF,[12, nyr]);
Fmin = min(FTMn);
Fmax = max(FTMn);

% Get data from Bamer's file by components
% for different FWFlux components
fnm = sprintf('%sFWF17.v3.nc',pthG);
tmm = double(nc_varget(fnm,'TIME'));
Xgr = nc_varget(fnm,'lon');
Ygr = nc_varget(fnm,'lat');
Rt  = double(nc_varget(fnm,'runoff_tundra')); % tundra runoff, km3/mo
Rg  = double(nc_varget(fnm,'runoff_ice')); % GrIS runoff - meltwater
D   = nc_varget(fnm,'solid_ice'); % solid ice
LGr = nc_varget(fnm,'LSMGr'); % Greenland mask
TM  = datenum(1958,1,1)+tmm; 
DV  = datevec(TM);
nrc = length(TM);

% Total FW flux = D+Rg+Rt (no CAA, Svalbard here);
clear FWF iD iR*
for it=1:nrc
  dv = datevec(TM(it));
  fprintf('Greenland Runoff, %i/%2.2i/%2.2i\n',dv(1:3));
  dmm = squeeze(D(it,:,:));
  Di = dmm.*LGr; % Greenland mask
  dmm = squeeze(Rt(it,:,:));
  Rti = dmm.*LGr; % Gr tundra
  dmm = squeeze(Rg(it,:,:));
  Rgi = dmm.*LGr; % Greenland meltwater
  FWF(it,1) = nansum(nansum(Di))+...
      nansum(nansum(Rti))+...
      nansum(nansum(Rgi));
  iD(it,1) = nansum(nansum(Di)); % total ice discharge, km3/mo
  iRt(it,1)= nansum(nansum(Rti));
  iRg(it,1)= nansum(nansum(Rgi));
end  

TY = [0:nrc-1]/12+DV(1,1) ;
% Do yearly fluxes:
YR = [DV(1,1):DV(end,1)];
nyr = nrc/12;
A = reshape(iD,[12 nyr]);
sD = sum(A);  % annual discharge - same as in saved fmat
A = reshape(iRt,[12 nyr]);
sRt = sum(A); % annual tundra - same as in saved fmat
A = reshape(iRg,[12 nyr]);
sRg = sum(A); % annual Greenland meltwater

FG = [sD;sRt;sRg]'; % all components together: solid disch, tundra, meltw.




YRs=[1993:2016];
FWFs=FWFn(36:end);

xx=[-8:15];
y0=mean(FWFs(1:5));
dY=mean(FWFs(end-8:end))-mean(FWFs(1:5));
yy=y0+dY/2+(dY/2)*tanh(xx);
%yy=dY/2*tanh(xx);

% Fit anomaly using logit function:
Fmn=mean(FWFn(1:imn));
dFW=FWFs-Fmn;
dF0=mean(dFW(end-8:end));
tm=[1:length(FWFs)];
m=1/2;
b=-4;
y=dF0*exp(m*tm+b)./(1+exp(m*tm+b));

% Fit linear regression to dFW to get 
% annual increase rate:
% Assuming GFW flux accelerating, want to find 
% acceleration rate
YY=[ones(length(FWFs),1),tm'];
dFW=dFW(:);
[BB,bint,RR,rint,STT]=regress(dFW,YY);
rg=YY*BB;

figure(4); clf;
axes('Position',[0.1 0.5 0.85 0.42]);
hold
plot(YRs,dFW,'k-','Linewidth',2);
plot(YRs,rg,'k-','Color',[0.6 0.6 0.6],'Linewidth',2);
set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'xtick',[1951:1:2016],...
	'ytick',[0:100:1500],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
title(sprintf('GFWA flux change, km3/yr^2, slope=%6.2f',BB(2)));
btx='fitFn_GreenlandRunoff.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.35 0.05]);

% Fit linear regression to dFW 
% for fraction of the total GFWA
% - part that is fluxed into Subpolar Gyre
% estimated based on Gr runoff (Labr + Irm segment)
% + net flux through Davis Str + Denmark Str.
%
%to get 
% annual increase rate:
% Assuming GFW flux accelerating, want to find 
% acceleration rate
% Based on Gr runoff and fluxes through straits, estimated
% Labr-Irm FW flux = 70 km3/yr + 62 km3/yr (Denmark) + 35 (Davis)=167 km3/yr
% ~ 80%  of total (208.6)
fr=0.8;
dFWfr=fr*dFW(:);
[BBf,bint2,RR2,rint2,STT2]=regress(dFWfr,YY);
rgf=YY*BBf;

keyboard


figure(5); clf;
axes('Position',[0.1 0.5 0.85 0.42]);
hold
plot(YRs,dFWfr,'k-','Linewidth',2);
plot(YRs,rgf,'k-','Color',[0.6 0.6 0.6],'Linewidth',2);
set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'xtick',[1951:1:2016],...
	'ytick',[0:100:1500],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
title(sprintf('SPG fraction, GFWA flux change, km3/yr^2, slope=%6.2f',BBf(2)));
bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.35 0.05]);

% Plot GFWA and GFWA-SPG fraction together
figure(6); clf;
axes('Position',[0.1 0.5 0.85 0.42]);
hold
% Total GFWA
plot(YRs,dFW,'k-','Linewidth',2.2);
plot(YRs,rg,'k--','Linewidth',1.6);

% Fraction of GFWA
plot(YRs,dFWfr,'k-','Color',[0.7 0.7 0.7],'Linewidth',2.2);
plot(YRs,rgf,'k--','Color',[0.7 0.7 0.7],'Linewidth',1.6);
set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'xtick',[1951:1:2016],...
	'ylim',[-50 550],...
	'ytick',[0:100:1500],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
title('GFWA and fraction fluxed into SPG, km3/yr');

axes('Position',[0.1 0.2 0.6 0.2]);
txt1=sprintf('GFWA:      %5.2f + %5.2f*t',BB(1),BB(2));
txt2=sprintf('GFWA SPG:  %5.2f + %5.2f*t',BBf(1),BBf(2));
text(0.1,0.5,txt1,'Fontsize',14);
text(0.1,0.3,txt2,'Fontsize',14);
set(gca,'xlim',[0.1 0.9],...
	'ylim',[0.2 0.7],...
	'visible','off');

bottom_text(btx,'pwd',1,'Position',[0.08 0.15 0.35 0.05]);



f_tan=0;
if f_tan==1
figure(1); clf;
axes('Position',[0.1 0.55 0.85 0.38]);
plot(dFW);
hold
plot(y);
title('Analytical Function Fit to GFWF anomaly: dF0*e^{m*t+b}/(1+e^{m*t+b})');

axes('Position',[0.1 0.09 0.85 0.38]);
%plot(YRn,FWFn,'k');
hb = bar(YRs,FWFs,0.9);
hold on
plot(YRs,yy,'r-');
set(hb,'edgecolor','none',...
       'FaceColor',[0.7 0.7 0.7]);

set(gca,'tickdir','out',...
	'xlim',[1993 yr2+0.5],...
	'xtick',[1950:2:2016],...
	'ytick',[0:250:1500],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',12);
title('tanh(t) fit to Data v3, Total Greenl. FW Flux, km3/yr');


% Plot solution
% Constant F=F0;
t=[0:24];
k=1/2; % tau=1/k - lifetime of Greenland Fr Water
F0=208.6; % mean GFWA flux, km3/yr
V=F0/k*(1-exp(-k*t));


btx = 'fitFn_GreenlandRunoff_v3.m';
bottom_text(btx,'pwd',1,'position',[0.05 0.04 0.4 0.08]);

end


