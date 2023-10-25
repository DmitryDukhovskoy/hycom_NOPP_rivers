% Plot monhtly statistics of CICE output 
% extracted in monthly_seaice.m
%
% Area of sea ice, volume of sea ice, 
% Bering/Fram sea ice area/vol fluxes
%
% Arctic Ocean 0.04 HYCOM-CICEv5 GOFS3.5
% Greenland runoff, no passive tracers
% analyze monthly mean sea ice fields
% extract from daily instantenous output
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 22;
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat = pthout;
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

btx = 'stat_month_seaice.m';


ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

c1=0;
c2=4;
pos=[0.12 0.1 0.82 0.82];
xl1= 50;
xl2= nn;
yl1= 500;
yl2= 4000;
%CMP = colormap_PBYR(200,c1,c2);
CMP = create_colormapBGY(200,c1,c2);
cmp = CMP.colormap;


hibin = [0:1:7]; % ice H bins
hsbin = [0:0.05:0.6]; %h snow

% Plot maps:
yr1=2016;
yr2=2017;
icc=0;
for YR=yr1:yr2
  pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%i_cice/',expt,YR);

  for im=1:12
    fmatout = sprintf('%scice_monthly_%i%2.2i.mat',pthmat,YR,im);
    if ~exist(fmatout), continue; end;
    fprintf('Loading %s\n',fmatout);
    load(fmatout);

    icc=icc+1;

    ci = ICEMO.Ice_conc;       % conc
    si = ICEMO.Snow_thck;      % snow  thck
    hi = ICEMO.Ice_thck;       % ice thck
    hi(hi==0)=nan;
    ci(HH>=0)=nan;
    si(si==0)=nan;

% Sea ice area: total sea ice
    A = nansum(nansum(ci.*Acell*1e-6)); % km2
    V = nansum(nansum(hi.*Acell*1e-9)); % km3
    I = find(~isnan(hi));
    NN = histcounts(hi(I),hibin); 
    J = find(~isnan(si));
    NS = histcounts(si(J),hsbin);

    ISTAT.Area(icc) = A;   
    ISTAT.Vol(icc)  = V;
    ISTAT.hist_hi(icc,:) = NN;
    ISTAT.hist_hibins = hibin;
    ISTAT.hist_hs(icc,:) = NS;
    ISTAT.hist_hsbin  = hsbin;
 
  end
end

dx=0.4;
dy=0.35;
POS = [0.08 0.6 dx dy; ...
       0.55 0.6 dx dy; ...
       0.08 0.15 dx dy; ...
       0.55 0.15 dx dy];

% Sea ice 
YR = [2016; 2017];
Nyr=2;
figure(1); clf;
for ik=1:Nyr
  i1=(ik-1)*12+1;
  i2=i1+11;
  if ik==2, i2=19; end;
  AA = ISTAT.Area(i1:i2);

  pos=POS(ik,:);
  axes('Position',pos);
  hba = bar(AA,0.98);
  set(hba,'FaceColor',[0 0.8 0.5]);
  ttl = sprintf('Ice Area, km^2, %i',YR(ik));

  set(gca,'xlim',[0 13],...
          'xtick',[1:12],...
          'ylim',[0 12e6],...
          'ytick',[0:2e6:12e6],...
          'tickdir','out',...
          'Fontsize',12);
  ylabel('km^2');

  title(ttl);
 
end

bottom_text(btx,'pwd',1);

%
% Sea ice volume
%
figure(2); clf;
for ik=1:Nyr
  i1=(ik-1)*12+1;
  i2=i1+11;
  if ik==2, i2=19; end;
  AA = ISTAT.Vol(i1:i2);

  pos=POS(ik,:);
  axes('Position',pos);
  hba = bar(AA,0.98);
  set(hba,'FaceColor',[0.3 0. 0.9]);
  ttl = sprintf('Ice Vol, km^3, %i',YR(ik));

  set(gca,'xlim',[0 13],...
          'xtick',[1:12],...
          'ylim',[0 2.5e4],...
          'ytick',[0:0.5e4:2.5e4],...
          'tickdir','out',...
          'Fontsize',12);
  ylabel('km^2');

  title(ttl);
 
end

bottom_text(btx,'pwd',1);




  hmi = ISTAT.hist_hi(i1:i2);
  hbe = ISTAT.hist_hibins;
  hbx = hbe(1:end-1)+(hbe(2)-hbe(1))/2;
  hbr = bar(hbx,hmi);




