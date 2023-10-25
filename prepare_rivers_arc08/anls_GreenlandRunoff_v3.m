% Plot Bamber river data - updated, version 3
% km3/mo
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

% Calculate standard error of the mean
% = std(x)/sqrt(n)
sgm=std(FWFn(1:imn));
sem=sgm/sqrt(imn); 


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

Ftot_mo=[]; % keep total FWFlux by months
TY = [0:nrc-1]/12+DV(1,1) ;
% Do yearly fluxes:
YR = [DV(1,1):DV(end,1)];
nyr = nrc/12;
A = reshape(iD,[12 nyr]);
Ftot_mo=A;
sD = sum(A);  % annual discharge - same as in saved fmat
A = reshape(iRt,[12 nyr]);
Ftot_mo=Ftot_mo+A;
sRt = sum(A); % annual tundra - same as in saved fmat
A = reshape(iRg,[12 nyr]);
Ftot_mo=Ftot_mo+A;
sRg = sum(A); % annual Greenland meltwater

FG = [sD;sRt;sRg]'; % all components together: solid disch, tundra, meltw.





figure(1); clf;
axes('Position',[0.1 0.55 0.85 0.35]);
%plot(YRn,FWFn,'k');
hb = bar(YRn,FWFn,0.9);
set(hb,'edgecolor','none',...
       'FaceColor',[0.7 0.7 0.7]);
hold on
plot([YRn(1)-1 YRn(end)+1],[fmn93 fmn93],'k','Color',[0.3 0.3 0.3]);

% This is km3/mo!
%for kk=1:nyr
%  y0=YRn(kk);
%  plot([y0 y0],[Fmin(kk) Fmax(kk)],'k');
%  plot([y0-0.1 y0+0.1],[Fmin(kk) Fmin(kk)],'k');
%  plot([y0-0.1 y0+0.1],[Fmax(kk) Fmax(kk)],'k');
%end

set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+0.5],...
	'xtick',[1950:5:2016],...
	'ytick',[0:250:1500],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
title('Data v3, Total Greenl. FW Flux, km3/yr');

% Plot cumulative anomaly
dF = FWFn-fmn93;
cdF = cumsum(dF);
dF93= dF(imn:end); % anomaly 1993-end
% 

axes('Position',[0.1 0.09 0.85 0.35]);
%plot(YRn,FWFn,'k');
hb = bar(YRn,cdF,0.9);
set(hb,'edgecolor','none',...
       'FaceColor',[0.7 0.7 0.7]);
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+0.5],...
	'xtick',[1950:5:2016],...
	'ytick',[-1000:1000:6000],...
	'yminortick','on',...
	'ylim',[-500 6000],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
stl=sprintf('cum dlt(FWFlux), wrt to 1958-1993 mean %6.1f km3/yr',fmn93);
title(stl);

btx = 'anls_GreenlandRunoff_v3.m';
bottom_text(btx,'pwd',1);
drawnow

% Plot seasonal signal, km3/mo
figure(2); clf;
axes('Position',[0.1 0.35 0.7 0.55]);
hbb=boxplot(FTMn');
set(hbb,{'linew'},{1.6});

set(gca,'tickdir','out',...
	'xlim',[0 13],...
	'xtick',[1:12],...
	'ytick',[0:50:500],...
	'xgrid','on',...
	'ygrid','on',...
	'ylim',[0 455],...
	'Fontsize',16);
title('Monthly Total FWFlux, km3/mo');
btx = 'anls_Greenland_runoff.m';
bottom_text(btx,'pwd',1,'Position',[-0.05 0.21 0.9 0.1]);
drawnow

if s_fig==1
  fgnm=sprintf('%sGreenland_FWF_Bamb_box',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end

%
%
% Plot bar stacked diagram
% total and cumulative anomalies:
% Anomalies by components
dD  = sD-idsch_mn;
dRt = sRt-tundra_mn;
dRg = sRg-meltw_mn;

cdD  = cumsum(dD);  % ice disch
cdRt = cumsum(dRt); % tundra runoff
cdRg = cumsum(dRg); % meltwater
cAnm = [cdD', cdRt', cdRg'];


cmp = [0, 0, 0; 0.5, 0.5, 0.5; 0.85, 0.85, 0.85];
figure(3); clf;
axes('Position',[0.1 0.55 0.85 0.35]);
hb2 = bar(YRn,FG,'stacked');
set(hb2,'edgecolor','none');
for ik=1:3
  hb2(ik).FaceColor = cmp(ik,:);
end
hold on
plot([YRn(1)-1 YRn(end)+1],[fmn93 fmn93],'k--','Color',[0.2 0.2 0.2]);
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+0.5],...
	'xtick',[1950:5:2016],...
	'ytick',[0:250:1500],...
	'xgrid','on',...
	'ygrid','on');
stl=sprintf('Greenl FWFlux Comp, 1958-1993 mean %6.1f km3/yr',fmn93);
title(stl);

lg = legend('Disch.', 'Tundra', 'Meltw');
set(lg,'Position',[0.084 0.91 0.086 0.073]);


% Plot cumulative anomaly
dF = FWFn-fmn93;
cdF = cumsum(dF);
axes('Position',[0.1 0.09 0.85 0.35]);
%plot(YRn,FWFn,'k');
%hb = bar(YRn,cdF,0.9); % total cum anomaly
%set(hb,'edgecolor','none',...
%       'FaceColor',[0.7 0.7 0.7]);
hb3 = bar(YRn,cAnm,'stacked');
set(hb3,'edgecolor','none');
for ik=1:3
  hb3(ik).FaceColor = cmp(ik,:);
end

set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+0.5],...
	'xtick',[1950:5:2016],...
	'ytick',[-1000:1000:6000],...
	'yminortick','on',...
	'ylim',[-500 6000],...
	'xgrid','on',...
	'ygrid','on');
stl=sprintf('cum dlt(FWFlux), wrt to 1958-1993 mean %6.1f km3/yr',fmn93);
title(stl);

btx = 'anls_GreenlandRunoff_v3.m';
bottom_text(btx,'pwd',1);

% Plot cumulative anomaly for 1993-2016
% Fraction of meltwater wrt to total Gr FWF anom:
frmlt=cAnm(:,end)./cdF';
frmlt(frmlt>1)=1;

figure(4); clf; 
axes('Position',[0.1 0.1 0.87 0.35]);
%plot(YRn,FWFn,'k');
%hb = bar(YRn,cdF,0.9); % total cum anomaly
%set(hb,'edgecolor','none',...
%       'FaceColor',[0.7 0.7 0.7]);
hb3 = bar(YRn,cAnm,'stacked');
set(hb3,'edgecolor','none');
for ik=1:3
  hb3(ik).FaceColor = cmp(ik,:);
end

% fraction of meltwater:
for ipp=imn:length(frmlt);
  stt=sprintf('%3.2f',frmlt(ipp));
  text(YRn(ipp)-0.4,cdF(ipp)+300,stt,'Fontsize',12);
end
text(1993,1200,'Fraction of meltwater','Fontsize',12);

set(gca,'tickdir','out',...
	'xlim',[1992.5 yr2+0.5],...
	'xtick',[1990:2016],...
	'ytick',[-1000:1000:6000],...
	'yminortick','on',...
	'ylim',[-500 6000],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
stl=sprintf('cum dlt(FWFlux), wrt to 1958-1993 mean %6.1f km3/yr',fmn93);
title(stl);

lg = legend('Disch.', 'Tundra', 'Meltw');
set(lg,'Position',[0.084 0.5 0.086 0.073]);

btx = 'anls_GreenlandRunoff_v3.m';
bottom_text(btx,'pwd',1);







