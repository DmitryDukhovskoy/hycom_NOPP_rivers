% Plot time series of annual
% tracer content converted 
% to the FWCanomaly= Greenland FWFlux anomaly
% accumulated in the Subpolar Gyre
%
% Tracer content calculated in vol_intgr_regn_trcr008.m
%
% Note: updated regions, now all regions
% are adjacent to each other
% So that Subpolar Gyre
% combines reg #1 (Labr)+Reg#2(IcelSea)+Reg#7 (Centr NAtl)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig=0;
regn = 'ARCc0.08';
expt = 110;  

s_mat = 0; % =0 - load saved time series of FWC anomaly
nvl   = 41;  % # of v. layers in the model
nTr   = 1;   % tracer to plot

% Experiments:
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);
btx = 'tracer_sbpGyre.m';

%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


% Define Regions:
%f_pltbox=0;
%BX = sub_define_boxes(HH,LON,LAT,f_pltbox);
%nbx = 6; % only N. Atl
% create mask for each box
%[XX,YY] = meshgrid((1:nn),(1:mm));
%for ib=1:nbx
%  iBG = BX(ib).IJ;
%  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
%  IN = find(INP==1 & HH<h0);
%  BX(ib).IN = IN;
%  fprintf(' Tracer integrated for region: %i %s\n',ib,BX(ib).Name);
%end

% Greenland FW F anomaly
frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat);
fprintf('Loading %s\n',frv);
load(frv); % annual cumulative Greenland FW anomaly cFWF


fmatgr=sprintf('%strcr_sbpGr_fwc.mat',pthmat);
if s_mat==1
  rbx=0;
  % Combine all years:
  cc=0;
  Vfw=0;
  for iyr=1993:2016
    yr   = iyr;
    fmat = sprintf('%strcr%2.2i_regn_VrtLrMass_%i.mat',pthmat,nTr,iyr);
    fprintf('Loading %s\n',fmat);
    load(fmat);
    nbx=length(VTRCR);

    cc=cc+1;
    for ibx=1:nbx
      nm=VTRCR(ibx).Name;
      dzm=VTRCR(ibx).DZM;
      mtrlr=VTRCR(ibx).TR;  % mass tracer in layers, in the box
      mTRd=nansum(mtrlr,1); % integrated over the whole depth & Box, kg

      ism = find(Ygr==iyr);
      fwf0 = cFWF(ism); % km3
      VFW(ibx).Name=nm;
      VFW(ibx).Year(cc)=yr;
      
  % Convert to the FWC anomaly  
      vfwl_mo=zeros(41,12);
      vfw_mo=0;
      frc_mo=vfwl_mo;
      for imo=1:12
        fmatT = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,iyr,imo);
        fprintf(' :: Loading %s\n',fmatT);
        load(fmatT);
        ibtm=5; % over whole depth to bottom
        Tr_dom=squeeze(TRCR(ibtm).MassTr_kg); % tracer integrated over whole water depth
        MTr_dom = nansum(nansum(Tr_dom)); % whole domain
        rr=mTRd(imo)/MTr_dom; 
%        rrl=mtrlr(:,imo)./MTr_dom; % by layers
        vfw_mo(imo)=fwf0*rr; % km3 
%        vfwl_mo(:,imo)=fwf0*rrl; 
%      frc_mo(:,imo)=rrl;
      end
  % Estimate volume of Greenland surplus FW in grid cell
  % = ratio of FW tracer*total FWF (integrated over time 1990-date requested)
%    MtarL(:,cc)=nanmean(mtrlr,2);  % mass Tracer, annual by layers
      VFW(ibx).Vfw(cc) = mean(vfw_mo); % km3 Gr FWF anomaly integrated in time in the box
    end

  end % year
  
  fprintf('Saving %s\n', fmatgr);
  save(fmatgr,'VFW');
else
  fprintf('Loading %s\n', fmatgr);
  load(fmatgr);
end

% Time series of Greenland FWF cumulative
ism = find(Ygr==1993);
fwf0 = cFWF(ism:end); % km3


yr   = 2016;
fmat = sprintf('%strcr%2.2i_regn_VrtLrMass_%i.mat',pthmat,nTr,yr);
fprintf('Loading %s\n',fmat);
load(fmat);
nbx=length(VTRCR);


figure(1); clf;
for ibx=1:nbx
  nm=VFW(ibx).Name;
  YR=VFW(ibx).Year;
  Vfw=VFW(ibx).Vfw;

  subplot(4,2,ibx);
  plot(YR,Vfw);
  stl=sprintf('%s\n',nm);
  title(stl);
end

% Fraction of total FW
cll=0;
if cll==1
CLR=[0 0.6 0.9; ...
     0 0.2 0.6; ...
     0 0.9 0.5; ...
     0.3 0.7 0.1;...
     0.9 0.8 0;...
     0.9 0.4 0;...
     0.9 0 0.9; ...
     0.6 0.2 0.1; ...
     0.7 0.7 0.7; ...
     0 0 0];
end
CLR=[0.47 0. 0.8;
     0.3 0.7 1;
     .72 0.15 0.15;
     0.64 0.49 0.;
     .75 0. 0.91;
     0.04 0.47 0;
     0.12 1 0.137;
     0 0 0;
     0.85 1 0.9;
     0.95 0.7 1];

% Plot FWC fraction wrt to 
% cumulative GFWA
figure(2); clf;
Vtot=0;
Vend=0;
for ibx=1:nbx
  nm=VFW(ibx).Name;
  YR=VFW(ibx).Year;
  Vfw=VFW(ibx).Vfw; % km3

  subplot(4,2,ibx);
  clr=CLR(ibx,:);
  plot(YR,Vfw./fwf0,'Color',clr,'linewidth',2);
  stl=sprintf('FW%%, %s\n',nm);
  title(stl);
  Vtot=Vtot+Vfw;
  yl2=ceil(max(Vfw./fwf0)*1e3)/1e3;
  dyy=yl2/4;
  set(gca,'tickdir','out',...
	  'xlim',[1993 2016],...
	  'ylim',[0 yl2],...
	  'ytick',[0:dyy:yl2],...
	  'xtick',[1990:2:2020],...
	  'xgrid','on',...
	  'ygrid','on',...
	  'Color','none',...
	  'YColor',clr,...
	  'Box','off');
  if mod(ibx,2)==0
    set(gca,'YaxisLocation','Right');
  end

  Vend(ibx)=Vfw(end)./fwf0(end); % fraction GFWA end of simulation
end
Vend(ibx+1)=Vtot(end)/fwf0(end);

% In the whole SPNA region:
ibx=nbx+1;
subplot(4,2,ibx);
clr=CLR(8,:);
plot(YR,Vtot./fwf0,'Color',clr,'linewidth',2);
yl2=ceil(max(Vtot./fwf0)*1e3)/1e3;
dyy=yl2/4;
set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'ylim',[0 yl2],...
	'ytick',[0:dyy:yl2],...
	'xtick',[1990:2:2020],...
	'xgrid','on',...
	'ygrid','on',...
	'Color','none',...
	'YColor',clr,...
	'Box','off');
if mod(ibx,2)==0
  set(gca,'YaxisLocation','Right');
end
stl=sprintf('FW%%, SPNA');
title(stl);
bottom_text(btx,'pwd',1);  

if s_fig==1
  fgnm=sprintf('%sarc08_%3.3i_fwc_fraction_basins.eps',pthfig,expt);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end

% FWC in terms of m of fresh water
%
figure(3); clf;
Mtot=0;
Atot=0;
Argn=0;
Mfw_end=0;
for ibx=1:nbx
  nm=VFW(ibx).Name;
  YR=VFW(ibx).Year;
  Vfw=VFW(ibx).Vfw*1e9; % m3
  IN=VTRCR(ibx).IN;
  Acl=sum(Acell(IN)); % m2
  Argn(ibx)=Acl*1e-6; % km2
  mfw=Vfw./Acl;
  Mtot=Mtot+Vfw;
  Atot=Atot+Acl;
  Mfw_end(ibx)=mfw(end);
  Bname{ibx}=nm(1:4);
  
  subplot(4,2,ibx);
  clr=CLR(ibx,:);
  plot(YR,mfw,'Color',clr,'linewidth',2);
  stl=sprintf('FWC, m, %s\n',nm);
  title(stl);
  Vtot=Vtot+Vfw;
  set(gca,'tickdir','out',...
	  'xlim',[1993 2016],...
	  'ylim',[0 1.05*max(mfw)],...
	  'xtick',[1990:2:2020],...
	  'xgrid','on',...
	  'ygrid','on',...
	  'Color','none',...
	  'Box','off');
  if mod(ibx,2)==0
    set(gca,'YaxisLocation','Right');
  end
  
end
Mtot=Mtot/Atot;
Mfw_end(ibx+1)=Mtot(end);
Bname{ibx+1}='SPNA';

% In the whole SPNA region:
ibx=nbx+1;
subplot(4,2,nbx+1);
clr=CLR(8,:);
plot(YR,Mtot,'Color',clr,'linewidth',2);
set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'ylim',[0 1.05*max(Mtot)],...
	'xtick',[1990:2:2020],...
	  'xgrid','on',...
	  'ygrid','on',...
	'Color','none',...
	'Box','off');
  if mod(ibx,2)==0
    set(gca,'YaxisLocation','Right');
  end
stl=sprintf('FWC, m, SPNA');
title(stl);
bottom_text(btx,'pwd',1);  

if s_fig==1
  fgnm=sprintf('%sarc08_%3.3i_fwc_basins.eps',pthfig,expt);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end

fprintf('============================\n');
fprintf('Region   Area(km2)  dFWC(m)  FWC%%\n');
fprintf('============================\n');
for ibx=1:nbx
  nm=VFW(ibx).Name;
  fprintf('%s %5.2d %4.2f %6.3f%%\n',nm(1:5),Argn(ibx),Mfw_end(ibx),Vend(ibx)*100);
end
fprintf('%s %5.2d %4.2f %6.3f%%\n','SPNA',Atot*1e-6,Mfw_end(ibx+1),Vend(ibx+1)*100);
fprintf('============================\n');

figure(4); clf;
axes('Position',[0.08 0.65 0.5 0.3]);
hold on;
%hb=bar(Mfw_end);
di=0.45;
for ibx=1:nbx+1
  i1=ibx-di;
  i2=ibx+di;
  j1=0;
  j2=Mfw_end(ibx);
  clr=CLR(ibx,:);
  hb=patch([i1 i1 i2 i2],[j1 j2 j2 j1],clr);
  set(hb,'Edgecolor','none');
end

set(gca,'tickdir','out',...
	'xtick',[1:8],...
	'xticklabel',Bname,...
	'xlim',[0.5 8.5],...
	'ylim',[0 0.5],...
	'ytick',[0:0.1:0.5],...
	'Fontsize',14);
title('dlt FWC from GFWA, m')

axes('Position',[0.08 0.12 0.5 0.3]);
hold on;
%hb=bar(Mfw_end);
di=0.45;
for ibx=1:nbx+1
  i1=ibx-di;
  i2=ibx+di;
  j1=0;
  j2=Vend(ibx);
  clr=CLR(ibx,:);
  hb=patch([i1 i1 i2 i2],[j1 j2 j2 j1],clr);
  set(hb,'Edgecolor','none');
end

set(gca,'tickdir','out',...
	'xtick',[1:8],...
	'xticklabel',Bname,...
	'xlim',[0.5 8.5],...
	'ylim',[0 0.5],...
	'ytick',[0:0.1:0.5],...
	'Fontsize',14);
title('Fraction of GFWA')

bottom_text(btx,'pwd',1);  

if s_fig==1
  fgnm=sprintf('%sarc08_%3.3i_GFWA_bars_basins.eps',pthfig,expt);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end






a=0;
if a==1
% Plot FWC for constant runoff:
Mup=sum(MtrL(1:25,:));
Mtot=sum(MtrL,1);

figure(10); clf;
plot(Mup/max(Mup));
hold on

% analytical solution assuming constant FWFlux=F0
dF0=280; % mean GFWF anom, km3/yr
t=[0:24];
k=1/5;
V=dF0/k*(1-exp(-k*t));
plot(V/max(V));


x=[0:24];
k=1/2;
b=-6;
y=dF0*exp(k*x+b)./(1+exp(k*x+b));
plot(y);

c0=exp(-12)/(k+2);
V2=exp(-12)*exp(2*t)/(k+2)+exp(-k*t)*c0;

end



