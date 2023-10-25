% Plot Extracted time series of sea ice area & volume
% from test experiments
% extracted in extr_ice_area_vol.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

pthout   = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat   = pthout;

fmatout = sprintf('%sseaice_VolArea_tests.mat',pthmat);
fprintf('Loading %s\n',fmatout);
load(fmatout);

CLR = [0.9 0.5 0; ...
       0.8 1 0.6; ...
       0   0.6 0.9; ...
       0.6 0   0.9; ...
       1  0.  0.4; ...
       0  0.8   0.2];

%
% Monthly mean from NSIDC NOAA - extent (including holes inside the contour)
% Note that calculated CICE area is ice area (excluding holes)
Ainoaa=[8.6e6; ...
        5.6e6; ...
        3.6e6; ...
        3.3e6; ...
        5.2e6; ...
        8.0e6; ...
        10.3e6]; 

for ii=1:7
  mo=5+ii;
  Tnoaa(ii,1) = datenum(2017,mo,1);
  Tnoaa(ii,2) = datenum(2017,mo,30);
end;
          

nexpt = length(ICE);
figure(1); clf;
set(gcf,'Position',[931         783        1109         556]);
axes('Position',[0.09 0.5 0.85 0.4]);
hold on;

amx = 0;
amn =1e10;
for ixx = 1:nexpt
  TM = ICE(ixx).TM;
  Ai = ICE(ixx).Area_km2;
  Vi = ICE(ixx).Vol_km3;

  DV  = datevec(TM);
  YR  = DV(1,1);
  dJ1 = datenum(YR,1,1);
  jD  = TM-dJ1+1;

  clr = CLR(ixx,:);
  plot(jD,Ai,'-','Color',clr,'Linewidth',2);

  amx = max([amx,max(Ai)]);
  amn = min([amn,min(Ai)]);


end;

for ii=1:length(Tnoaa)
		d1 = Tnoaa(ii,1)-dJ1+1;
		d2 = Tnoaa(ii,2)-dJ1+1;
		ai = Ainoaa(ii);

		plot([d1 d2],[ai ai],'k-','Linewidth',2);
end
   
amx = max([amx,max(Ainoaa)]);
amn = min([amn,min(Ainoaa)]);

icc=0; 
for imm=6:12
  dnmb = datenum(2017,imm,1);
  jd = dnmb-dJ1+1;
  plot([jd jd],[0 1.2*amx],'--','Color',[0.8 0.8 0.8]); 
  icc=icc+1;
  xtck(icc)=jd;
end

dnmb = datenum(2017,12,31);
jd = dnmb-dJ1+1;
plot([jd jd],[0 1.2*amx],'--','Color',[0.8 0.8 0.8]);

set(gca,'tickdir','out',...
        'xtick',xtck,...
        'xlim',[xtck(1) 365],...
        'ylim',[0.9*amn 1.1*amx],...
        'ytick',[0:2e6:1.1*amx],...
        'ygrid','on',...
        'fontsize',12);
title('Sea Ice Area / Extent, km^2');

% Legend
axes('Position',[0.1 0.1 0.6 0.3]);
hold on;
yy0=length(ICE)+1;
xx1=0.05;
xx2=xx1+0.2;
xx3=xx2+0.1;

for ixx=1:length(ICE);
  nm = ICE(ixx).Name;
  clr = CLR(ixx,:);
  iy = yy0-ixx+1;
  plot([xx1 xx2],[iy iy],'-','Color',clr,'Linewidth',2);
  text(xx3,iy,nm,'Fontsize',12,'Interpreter','none');
end
iy = iy-1;
plot([xx1 xx2],[iy iy],'-','Color',[0 0 0],'Linewidth',2);
text(xx3,iy,'NOAA IceExt','Fontsize',12);

set(gca,'xlim',[0 5],...
        'ylim',[1 yy0+1],...
        'visible','off')

btx = 'plot_ice_area_vol.m';
bottom_text(btx,'pwd',1);


   
figure(2); clf;
set(gcf,'Position',[925         775        1109         556]);
axes('Position',[0.09 0.5 0.85 0.4]);
hold on;

vmx = 0;
vmn =1e10;
for ixx = 1:nexpt
  TM = ICE(ixx).TM;
  Vi = ICE(ixx).Vol_km3;

  DV  = datevec(TM);
  YR  = DV(1,1);
  dJ1 = datenum(YR,1,1);
  jD  = TM-dJ1+1;

  clr = CLR(ixx,:);
  plot(jD,Vi,'-','Color',clr,'Linewidth',2);

  vmx = max([vmx,max(Vi)]);
  vmn = min([vmn,min(Vi)]);

end;

icc=0;
for imm=6:12
  dnmb = datenum(2017,imm,1);
  jd = dnmb-dJ1+1;
  plot([jd jd],[0 1.2*vmx],'--','Color',[0.8 0.8 0.8]);
  icc=icc+1;
  xtck(icc)=jd;
end

dnmb = datenum(2017,12,31);
jd = dnmb-dJ1+1;
plot([jd jd],[0 1.2*vmx],'--','Color',[0.8 0.8 0.8]);

set(gca,'tickdir','out',...
        'xtick',xtck,...
        'xlim',[xtck(1) 365],...
        'ylim',[0 1.1*vmx],...
        'ytick',[0:5000:1.2*vmx],...
        'ygrid','on',...
        'fontsize',12);
title('Ice Volume, km^3');

% Legend
axes('Position',[0.1 0.1 0.6 0.3]);
hold on;
yy0=length(ICE)+1;
xx1=0.05;
xx2=xx1+0.2;
xx3=xx2+0.1;

for ixx=1:length(ICE);
  nm = ICE(ixx).Name;
  clr = CLR(ixx,:);
  iy = yy0-ixx+1;
  plot([xx1 xx2],[iy iy],'-','Color',clr,'Linewidth',2);
  text(xx3,iy,nm,'Fontsize',12,'Interpreter','none');
end
iy = iy-1
plot([xx1 xx2],[iy iy],'-','Color',[0 0 0],'Linewidth',2);
text(xx3,iy,'NOAA IceExt','Fontsize',12);

set(gca,'xlim',[0 5],...
        'ylim',[1 yy0+1],...
        'visible','off')

btx = 'plot_ice_area_vol.m';
bottom_text(btx,'pwd',1);





