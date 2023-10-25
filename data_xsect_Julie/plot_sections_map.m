% Vertical sections output archive files:
% for most recent sensit. experiments with CICE5
% old simulations 
%
% Fields are long-term average over specified time
% 
% Vertical transects of T and S
% Across the Arctic, Bering-to-North Pole-to-Nordic Seas along 

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;

%xname = 'BerNatl'; % Bering - N.Atl
%xname = 'BeaufNatl'; % Beaufort Shelf - South Norway
%xname = 'LaptSea'; % Laptev Sea
xname = 'BeaufIcel';  % Beaufort Sea to Iceland via N Pole
%xname = 'BafnFram'; % section around Greenland from Baffin to Fram str 


pthout = '/nexsan/people/ddmitry/hycom/data_extracted/';

% ARCc0.04 010 similar to ARCc0.08 110 - tracers, river climatology, no Greenland
% ARCc0.04 012 similar to ARCc0.08-112 - tracers, river monthly data, Greenland on
% 022 - GOFS3.5 CICEv5

%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
% 
% 
% Specify segments of the x-section
XSCT = sub_xsections(xname,0.04);
IJs  = XSCT.IJs;
IIs  = IJs(:,1);
JJs  = IJs(:,2);


% Find the overall length of the whole section:
Isgm = XSCT.Isgm;
Jsgm = XSCT.Jsgm;
nsgm = size(Isgm,1)-1;
for ii=1:nsgm
  i0 = Isgm(ii);
  j0 = Jsgm(ii);
  i1 = Isgm(ii+1);
  j1 = Isgm(ii+1);
  lon0 = LON(j0,i0);
  lat0 = LAT(j0,i0);
  lon1 = LON(j1,i1);
  lat1 = LAT(j1,i1);

  Dsgm(ii) = distance_spheric_coord(lat0,lon0,lat1,lon1)*1e-3;  % spheric distance, km 
end
Ltot = sum(Dsgm);


nS=length(IIs);
clear Xl Yl
for ii=1:nS
  i0=IJs(ii,1);
  j0=IJs(ii,2);
  Xl(ii,1)=LON(j0,i0);
  Yl(ii,1)=LAT(j0,i0);
  Hb(ii,1)=HH(j0,i0);
end;
INDs=sub2ind(size(HH),JJs,IIs);

% Distance normalize and scale to match total length of the section
XL = sub_distance_section(Xl,Yl,Ltot);

% Allign smaller segments of POP with HYCOM
dL = 1000;
icc = 0;
for ll=0:dL:max(XL)
  dmm = abs(XL-ll);
  ii=find(dmm==min(dmm),1);
  i0=IIs(ii);
  j0=JJs(ii);
  icc = icc+1;
  Chck(icc,1) = i0;
  Chck(icc,2) = j0;
  Chck(icc,3) = LON(j0,i0);
  Chck(icc,4) = LAT(j0,i0);
  Chck(icc,5) = ll;
end
icc=icc+1;
i0 = IIs(end);
j0 = JJs(end);
Chck(icc,1) = i0;
Chck(icc,2) = j0;
Chck(icc,3) = LON(j0,i0);
Chck(icc,4) = LAT(j0,i0);
Chck(icc,5) = Ltot;

f_chck = 0;
if f_chck == 1
  save('pop_chck_pnts.mat','Chck');
end

%
% POP section
pthout = '/nexsan/people/ddmitry/data_POP_Julie/';
icc = 0;
YR=2017;
mo=1;
fmat_out = sprintf('%sCrossArctic_MonthlyAvg_%4.4i-%2.2i.mat',pthout,YR,mo);
PSCT = load(fmat_out);

XPl  = PSCT.lon;
YPl  = PSCT.lat;

I=find(XPl>180);
XPl(I)=XPl(I)-360;
%XPL = sub_distance_section(XPl,YPl,Ltot);
XPL = sub_scale_POPdist(XPl, YPl, Chck);

% Find HYCOM indices for POP section:
fprintf('Finding HYCOM indices for POP section ...\n');
np = length(XPl);
IPs = [];
JPs = [];
f_extr = 0;
if f_extr == 1
  i1=1;
  for ii=i1:np
    if mod(ii,50)==0,
      fprintf('  Processed %5.2f%% ...\n',ii/np*100);
    end
    xp = XPl(ii);
    yp = YPl(ii);

    dd = distance_spheric_coord(LAT,LON,yp,xp)*1e-3;
    [jmin,imin] = find(dd == min(min(dd)));
    JPs(ii,1) = jmin;
    IPs(ii,1) = imin;
  end
  save('pop_line_indx.mat','IPs','JPs');
else
  load('pop_line_indx.mat');
end

% There are messed up grid points around the NP in POP
% Need to cut out these segments:
XPlo = XPl;
YPlo = YPl;
XPLo = XPL;
IPso = IPs;
JPso = JPs;

IPs = IPs(:);
JPs = JPs(:);
XPL = XPL(:);
ii1 = 389;
ii2 = 407;
XPl = [XPl(1:ii1); XPl(ii2:end)];
YPl = [YPl(1:ii1); YPl(ii2:end)];
XPL = [XPL(1:ii1); XPL(ii2:end)];
IPs = [IPs(1:ii1); IPs(ii2:end)];
JPs = [JPs(1:ii1); JPs(ii2:end)];




% Plot:
fprintf('Plotting ...\n');
Lmsk = HH*0;
Lmsk(HH<0)=1;
cmsk = [0,0,0; 1,1,1];

btx = 'plot_sections_map.m';
figure('Position',[1000         483         752         833]); 
clf;
%set(gcf,'Position',[1000         483         752         833]);
hold on
%pcolor(Lmsk); shading flat;
%colormap(cmsk);
%caxis([0,1])
contour(HH,[0 0],'k','Linewidth',2);
contour(HH,[-5000:500:-100],'Color',[0.7 0.7 0.7]);
contour(HH,[-5000:1000:-100],'Color',[0.4 0.4 0.4]);
%contour(HH,[-5000:500:-100]);
colormap(jet(200));
plot(IIs,JJs,'b.-'); % HYCOM
plot(IPs,JPs,'g.-'); % POP
axis('equal');
title('HYCOM(blue) and POP(green) Sections');
%
% Plot distance marks along the section

dL = 250;
for ll=0:dL:max(XL)
  dmm = abs(XL-ll);
  ii=find(dmm==min(dmm),1);
  i0=IIs(ii);
  j0=JJs(ii);

  if mod(ll,1000)==0
    plot(i0,j0,'r.','Markersize',20);
    lts=sprintf('%i',ll);
    text(i0+5,j0,lts,'Fontsize',12);
  else
    plot(i0,j0,'r.','Markersize',12);
  end
end

for ll=0:dL:max(XPL)
  dmm = abs(XPL-ll);
  ii=find(dmm==min(dmm),1);
  i0=IPs(ii);
  j0=JPs(ii);

  if mod(ll,1000)==0
    plot(i0,j0,'m.','Markersize',20);
    lts=sprintf('%i',ll);
    text(i0+5,j0,lts,'Fontsize',12,'Color',[0.8,0,0.5]);
  else
    plot(i0,j0,'m.','Markersize',12);
  end
end

set(gca,'xlim',[500 3000],'ylim',[700 3800])
%set(gca,'xlim',[1400 1900],'ylim',[2400 2900])
%set(gca,'xlim',[1750 1900],'ylim',[2400 2570])
%set(gca,'xlim',[1800 1880],'ylim',[2460 2530])

bottom_text(btx,'pwd',1);




 
