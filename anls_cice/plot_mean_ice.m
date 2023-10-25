% Plot monthly mean conc/h CICE output - 
% extracted in mean_ice_h_conc.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt=110;

yr1=1993;
yr2=2016;
mo=9;
dd=7; 

% Domain: AO and SPNA
i1=35;
j1=246;
i2=1600;
j2=1916;

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_icemat/';

fmat    = sprintf('%sarc08_%3.3i_cice_mean_mnth%2.2i.mat',pthmat,expt,mo);

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
HH=HH(j1:j2,i1:i2);
[mm,nn]=size(HH);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2
Acell=Acell(j1:j2,i1:i2);

fprintf('Loading %s\n',fmat);
load(fmat);

% Specify Arctic Ocean region only:
IJ=[183         737
         311         678
         365         702
         436         766
         482         781
         526         803
         526         821
         707         849
         720         810
         746         720
         856         704
         907         728
        1039         711
        1080         704
        1107         643
        1241         504
        1402         357
        1556         355
        1558        1665
         157        1665];

[II,JJ]=meshgrid([1:nn],[1:mm]);
P=inpolygon(II,JJ,IJ(:,1),IJ(:,2));
IN=find(P==1 & HH<0);


nyr=length(ICE);
% Plot overall mean hi, ui
for ik=1:nyr
  hi=ICE(ik).hice;
  ai=ICE(ik).aice;
  ui=ICE(ik).uice;
  vi=ICE(ik).vice;
  if ik==1
    h0=hi*0;
    u0=ui*0;
    v0=vi*0;
    a0=ai*0;
  end
  
  h0=h0+hi;
  a0=a0+ai;
  u0=u0+ui;
  v0=v0+vi;
  
% Overall volume:
  Vol(ik)=nansum(hi(IN).*ai(IN).*Acell(IN));
% Oveall area:
  Area(ik)=nansum(ai(IN).*Acell(IN));
% Oveall ext:
  I=find(ai<0.15);
  ai(I)=nan;
  Ext(ik)=nansum(ai(IN).*Acell(IN));
end;

hi=h0/ik;
ai=a0/ik;
ui=u0/ik;
vi=v0/ik;

Lmsk=HH*0;
Lmsk(HH<0)=1;

cmp=flipud(colormap_cold(200));
figure(1); clf;
axes('Position',[0.1 0.2 0.8 0.7]);
pcolor(hi); shading flat;
colormap(cmp);
hold
contour(ai,[0.15 0.15],'r-','linewidth',1.6);
contour(hi,[1:10],'k-','linewidth',1);
caxis([0 4]);
axis('equal');
set(gca,'xtick',[],...
	'ytick',[],...
	'xlim',[1 nn],...
	'ylim',[250 mm],...
	'color',[0.8 0.8 0.8]);

cb=colorbar;
set(cb,'Position',[0.91 0.25 0.02 0.6],...
       'Fontsize',14,...
       'TickLength',0.035);

yr1=ICE(1).Year;
yr2=ICE(end).Year;
stl=sprintf('0.08 HYCOM-CICE %3.3i, mean Hice, %2.2i, %i-%i',...
	    expt,mo,yr1,yr2);
title(stl);

btx='plot_mean-ice.m';
bottom_text(btx,'pwd',1,'position',[0.05 0.15 0.4 0.04]);


% Plot sea ice volume and extent
Ext=Ext*1e-6;  % m2->km2
Vol=Vol*1e-9;  % m3-> km3

YR=[yr1:yr2];
e19=4.32e6; % NSIDC, Sept 2019 extent 
e16=4.17e6; % NSIDC, Sept 2016 extent 
e8310=6.41e6; % mean 1983-2010 sept
v19=4200; % km3, PIOMAS estimate
v8090=15000; % approximate mean vol 1980-1990, PIOMAS


figure(2); clf;
axes('Position',[0.08 0.6 0.85 0.32]);
hb=bar(YR,Ext);
set(hb,'Facecolor',[0.4 0.4 0.4]);
hold on;
plot([yr1-1 yr2+1],[e19 e19],'r--','linewidth',2,'Color',[1 0.5 0]);
plot([yr1-1 yr2+1],[e8310 e8310],'r--','linewidth',2,'Color',[0 1 0.5]);
set(gca,'tickdir','out',...
	'xlim',[yr1-0.5 yr2+0.5],...
	'xtick',[yr1:2:yr2],...
	'ylim',[0 8e6],...
	'Fontsize',12);
stl=sprintf('0.08 HYCOM-CICE %3.3i, mean extent km^2, %2.2i, %i-%i',...
	    expt,mo,yr1,yr2);
title(stl);

axes('Position',[0.08 0.1 0.85 0.32]);
hb=bar(YR,Vol);
set(hb,'Facecolor',[0. 0.4 0.6]);
hold on;
plot([yr1-1 yr2+1],[v19 v19],'r--','linewidth',2,'Color',[1 0.5 0]);
plot([yr1-1 yr2+1],[v8090 v8090],'r--','linewidth',2,'Color',[0 1 0.5]);
set(gca,'tickdir','out',...
	'xlim',[yr1-0.5 yr2+0.5],...
	'xtick',[yr1:2:yr2],...
	'ylim',[0 30e3],...
	'Fontsize',12);
stl=sprintf('0.08 HYCOM-CICE %3.3i, mean Volume km^3, %2.2i, %i-%i',...
	    expt,mo,yr1,yr2);
title(stl);

bottom_text(btx,'pwd',1,'position',[0.05 0.02 0.4 0.04]);


