% T17 - for ARCc0.04
% ------------------------------------------
% Compare topo difference in the relax zone
% of ARCc 
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

PTH.topo='/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
T = 17;     % topo #
TV='17DD';  % new name 

fina = sprintf('%sdepth_ARCc0.04_17.a',PTH.topo);
finb = sprintf('%sdepth_ARCc0.04_17.b',PTH.topo);
fida = fopen(fina,'r');
flgrd = sprintf('%sregional.grid.ARCc0.04',PTH.topo);

GRD = read_grid_bath(flgrd,fina);
HH = GRD.Topo;
I=find(HH>1e20);
HH=-HH;
HH(I)=100;

fina = sprintf('%sdepth_ARCc0.04_17DD.a',PTH.topo); % T17DD
GRD = read_grid_bath(flgrd,fina);
HH17 = GRD.Topo;
I=find(HH17>1e20);
HH17=-HH17;  % Correceted T17 (T17D)
HH17(I)=100;




[mg,ng]=size(HH17);

f_pltH = 1; % topo maps
if f_pltH
  figure(1); clf;
  axes('Position',[0.05 0.08 0.4 0.84]);
  pcolor(HH); shading flat;
  hold on
  contour(HH,[0 0],'k');
  caxis([-1000 0]);
  set(gca,'Color',[0.4 0.4 0.4]);
  title('ARCc0.04 T17');

  axes('Position',[0.55 0.08 0.4 0.84]);
  pcolor(HH17); shading flat;
  hold on
  contour(HH17,[0 0],'k');
  caxis([-1000 0]);
  set(gca,'Color',[0.4 0.4 0.4]);
  title('T17D');


  %txtbtm='/hycom_arc08/NOPP_rivers/compare_T07_T09.m';
  txtbtm='compare_ARCc_T17vsT17D.m';
  bottom_text(txtbtm,'pwd',1);
end

f_pltDH = 0; % plot difference H
if f_pltDH
  dH = HH17-HH;
  figure(1); clf;
  axes('Position',[0.09 0.1 0.86 0.82]);
  pcolor(dH); shading flat;
  caxis([-5000 0]);
  set(gca,'Color',[0.4 0.4 0.4]);
  title('Topo difference');
  colorbar


  %txtbtm='/hycom_arc08/NOPP_rivers/compare_T07_T09.m';
  txtbtm='compare_ARCc_T09vsT17.m';
  bottom_text(txtbtm,'pwd',1);
  
end  

[mm,nn] = size(HH);

% Plot section:
sct='Pacif_ob2';
switch(lower(sct))
 case('we_carct');
% Central Arctic:
  jj1=1430;
  jj2=jj1;
  ii1=200;
  ii2=1450;
 case('atl_ob')
% Atl. OB:
  jj1=5;
  jj2=jj1;
  ii1=190;
  ii2=1050;
 case('pacif_ob');
% Pacific OB:
  jj1=2516;
  jj2=jj1;
  ii1=1;
  ii2=nn;
 case('pacif_ob2');
% Pacific OB:
  jj1=mm;
  jj2=jj1;
  ii1=1;
  ii2=nn;
end

[mm,nn]=size(HH);

h=HH(mm-100:end,400);
h17=HH17(mm-100:end,400);
figure(2); clf;
axes('Position',[0.08 0.56 0.9 0.48]);
plot(h,'b','linewidth',2);
hold on;
plot(h17,'r','linewidth',1.6);
stt=sprintf('%s, HYCOM topo 09 (b) vs topo 17 (r)',sct);
title(stt);

axes('Position',[0.08 0.08 0.9 0.4]);
plot(h-h17);
title('H difference');

figure(3); clf;
HH(isnan(HH))=100;
HH17(isnan(HH17))=100;
contour(HH,[0 0],'b');
hold on;
contour(HH17,[0 0],'r');

