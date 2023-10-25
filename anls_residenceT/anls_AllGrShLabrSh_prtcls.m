% Analysis of particles - histograms, mean transit time etc.
%
%  Gr Shelf & Labrador Shelf Particles 
%
% Combine all experiments with particles:
% N experiments with particles
% released at different depth levels
% 
% Compare different levels
%
% Particles are not added during the simulation
% All N prticles seeded at once at initial state
% For each particle: 
% no T or S is tracked, only time and location
% 
addpath /usr/people/ddmitry/codes/MyMatlab;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

fget=0;  % =0 - load previously loaded data

NSM=[1,2,3,4,5,6]; % each epxeriment released 100 particles at 1 depth
%NSM=[1,2,3,4];
nnm=length(NSM);

SLR=[10,15,23,31]; % depth levels of the particels
ZLR=[50,90,150,450];  % nominal depths of particles in deep ocean
nlr=length(SLR);

%slr = 15; % vertical layer
%s_par=1; % parallel session 
%f_plt=1; % plot prtcles
s_fig = 0;
%cc = 71; % last SAVED frame - if need to restart from frame # 
      % the next frame will be cc+1
      % = 0 - start from beginning

rg = 9806; 

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/%3.3i/fig_GrGprtcl/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
pthmat2 = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_WLabr_prt/';

txtb = 'anls_AllGrShLabrSh_prtcls.m';

% ------------------------
% TOPO
% ------------------------
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

xlim1 = 220;
xlim2 = 1250;
ylim1 = 30;
ylim2 = 1100;

% layer interface depths:
% in N Atl:
zz=[                     0
                        -1
         -2.80000007967061
         -6.04000002390118
         -10.7199998326917
         -15.6499996414823
         -21.4599995777458
         -28.3299994502728
         -36.3299994502728
         -44.3299994502728
         -52.3299994502728
         -60.3299994502728
         -68.3299994502728
         -76.3299994502728
         -84.3299994502728
         -92.3299994502728
         -100.329999450273
         -108.329999450273
         -116.329999450273
         -124.329999450273
         -132.329999450273
         -140.329999450273
         -148.329999450273
         -156.329999450273
         -166.329999450273
         -182.730000087638
         -218.650001234894
         -261.030001362367
         -311.050001872259
         -370.070002382151
         -439.709999577746
         -521.889997793124
         -642.060033995449
         -833.055521452108
         -1631.30343089531
         -2130.04884186818
          -2911.1784563899
         -3125.48785879659
         -3385.88090896996
         -3385.88090896996
         -3385.88090896996
         -3385.88090896996];

f_dps = 0;
if f_dps==1
% Read layer thicknesses:
  pthbin = '/nexsan/hycom/ARCc0.08_112/data/1993/';
  fina = sprintf('%s112_archm.1993_002_12.a',pthbin);
  finb = sprintf('%s112_archm.1993_002_12.b',pthbin);

  [ZM,ZZ] = sub_zz_zm(fina, finb,HH);
% on the shelf depth of the layer:
 ZM(31,539,529)

  fspg='SPG_noNorth_indx.mat';  % get IGR indices of the region
  load(fspg);
  IGR(end+1,:)=IGR(1,:);

  [XX,YY] = meshgrid((1:nn),(1:mm));

  ING = inpolygon(XX,YY,IGR(:,1),IGR(:,2));
  Inn = find(ING==1 & HH<-500);

  SLR2 = [10; 15; 23; 29];
 
  for ik=1:length(SLR);
    izz=SLR2(ik); 
    z1 = squeeze(ZM(izz,Inn));
    zm_lrs(ik)=mean(mean(z1));
    zz1 = squeeze(ZZ(izz,Inn));
    zz2 = squeeze(ZZ(izz+1,Inn));
    zz_lrs(ik,1)=mean(mean(zz1));
    zz_lrs(ik,2)=mean(mean(zz2));
  end

end

%keyboard

YRPLT=[1993:2019];
nplt = length(YRPLT);

fnmout  = sprintf('%sGrShprt_comb.mat',pthmat);
fnmout2 = sprintf('%sWLabrShprt_comb.mat',pthmat2);


% Combine all experiments with Lagrangian particles together 
% GrShelf and Labr shelf separately
% for missing dates - use the previous number and locations of particles
PPG = sub_GrSh(fnmout);  % Gr shelf
PPL = sub_GrSh(fnmout2); % Labr shelf


[a1,a2]=size(PPG(1).Xp);
cmp=colormap_blue(a2);

% Colorcode for floats at different depths:
CLRZ=[0 0.5 0.8;...
     0.9 0.4 0;...
     0. 0.8 0.2;...
     0.5 0.1 1];  

Hmsk=HH*0;
Hmsk(HH<0)=1;

% Region of interest - double check
% with dS_FWC_timeseries_SubpolarGyre.m
fspg='SPG_noNorth_indx.mat';  % get IGR indices of the region
load(fspg);
IGR(end+1,:)=IGR(1,:);


%keyboard  
% Analyze # of particles in the SPNA by years
% SPNA region is defined in FWC_regions_GreenlExp.m
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% N Atl regions
%BX = sub_define_boxes(HH,LON,LAT,1);
[XX,YY] = meshgrid((1:nn),(1:mm));

Iyr = [1:27];
% Calculate probability of part. being inside SPNA at least for N years after release
PRBG = sub_prob_inSPNA(PPG,SLR,IGR,Iyr);
TMy=PRBG.TMy;

PRBL = sub_prob_inSPNA(PPL,SLR,IGR,Iyr);



% Same plot differently presented
dpx=0.38;
dpy=0.32;
POS = [0.09 0.6 dpx dpy; ...
       0.59 0.6 dpx dpy; ...
       0.09 0.15 dpx dpy; ...
       0.59 0.15 dpx dpy];

Ti=[0:0.02:max(Iyr)]';
T0=[0:max(Iyr)]';


% ==================
%
%  Age distribution
%
% =====================

% Histograms for GrSh and LabrSh:
%
bb = PRBG.PRb;
figure(10); clf;
set(gcf,'Position',[854 503 1687 828]);
for ilr=1:nlr
%  clr=CLRZ(ilr,:);
  pos=POS(ilr,:);
  axes('Position',pos);
  hold on;

  clr1=[0.5 0.7 1];
  clr2=[1 0.7 0.5];
  clr11=[0 0.4 0.8];
  clr22=[0.8 0.4 0];

  br1 = PRBG.PRyr(ilr,:);   % Prob T>= N years - something ? should be 1-CDF, CDF=P(T<=t)
%  br1 = PRBG.TRyr(ilr,:);
%  br1 = PRBG.PRb(ilr,:);   % true probability: P(T=t years)
% Normalize - not for cumulative density functions, 
%  br1=br1./sum(br1);
  bmax=max(br1);
 
  dx=0.4; 
  for ixx=1:length(br1);
    i1=Iyr(ixx)-dx;
    i2=Iyr(ixx);
    patch([i1 i1 i2 i2],[0 br1(ixx) br1(ixx) 0],clr1);
  end

  br1 = PRBL.PRyr(ilr,:);
% Normalize
%  br1=br1./sum(br1);
  bmax=max([bmax,max(br1)]);

  for ixx=1:length(br1);
    i1=Iyr(ixx);
    i2=Iyr(ixx)+dx;
    patch([i1 i1 i2 i2],[0 br1(ixx) br1(ixx) 0],clr2);
  end
% 
% Plot P(Age<=Nyrs):
% Get probability of Age <= N years:
  PRyr = PRBG.PRyr;
  Page = 1-PRyr(ilr,:);
  Page = Page(:);
  Page = [0; Page];
  PGi = interp1(T0,Page,Ti,'pchip');
  plot(Ti,PGi*bmax,'-','Color',clr11,'linewidth',2); 

  PRyr = PRBL.PRyr;
  Page = 1-PRyr(ilr,:);
  Page = Page(:);
  Page = [0; Page];
  PLi = interp1(T0,Page,Ti,'pchip');
  plot(Ti,PLi*bmax,'-','Color',clr22,'linewidth',2); 

  dxtk=0.3;
  xlim2=28-dxtk;
  ylim2=bmax*1.1;
  Ytick=[0:bmax*0.2:bmax];
  plot([xlim2 xlim2],[0 ylim2],'-',...
       'Linewidth',1,...
       'Color',[0 0 0]);
  for ik=2:length(Ytick)
    plot([xlim2 xlim2+dxtk],[Ytick(ik) Ytick(ik)],'-','Color',[0 0 0]);
    text(xlim2+1.02*dxtk,Ytick(ik),sprintf('%3.1f',(ik-1)*2/10));
  end

  stl = sprintf('Bar: P(Age>=N), Line: P(Age<=N), Lr %i, Z=%i',ilr,ZLR(ilr));
  title(stl);

  set(gca,'Tickdir','out',...
   'Fontsize',14,...
   'xlim',[0 xlim2+dxtk],...
   'xtick',Iyr,...
   'ylim',[0 ylim2])

% Save median and percentiles:
% of partilce ages:
% ~38% - flow from Baffin Bay and ~62% - flow from S Gr and Davis Strait
  i1=max(find(PGi<=0.25)); % 25%
  PRCTG(ilr,1)=Ti(i1);
  i1=max(find(PGi<=0.5));
  PRCTG(ilr,2)=Ti(i1);
  i1=max(find(PGi<=0.75));
  PRCTG(ilr,3)=Ti(i1);
  i1=max(find(PGi<=0.1)); % 10%
  PRCTG(ilr,4)=Ti(i1);
  i1=max(find(PGi<=0.9)); % 90%
  PRCTG(ilr,5)=Ti(i1);
  

  i1=max(find(PLi<=0.25)); % 25%
  PRCTL(ilr,1)=Ti(i1);
  i1=max(find(PLi<=0.5));
  PRCTL(ilr,2)=Ti(i1);
  i1=max(find(PLi<=0.75));
  PRCTL(ilr,3)=Ti(i1);
  i1=max(find(PLi<=0.1)); % 10%
  PRCTL(ilr,4)=Ti(i1);
  i1=max(find(PLi<=0.9)); % 90%
  PRCTL(ilr,5)=Ti(i1);
% 
% Weighted average of particles ages:
  wGr=0.7;
  for jj=1:3
    PRCT(ilr,jj)=wGr*PRCTG(ilr,jj)+(1-wGr)*PRCTL(ilr,1);
  end

%
% Estimate mean age and 95% CI:
% Central Limit thm, for large n
% mean(x) ~N(0,1)
  

end

axes('Position',[0.01 0.42 0.07 0.2]);
text(0.1,0.6,'WGrSh','Color',clr11,'Fontsize',14);
text(0.1,0.5,'NLabrSh','Color',clr22,'Fontsize',14);
set(gca,'xlim',[0.1 0.6],'ylim',[0.2 0.7]);
set(gca,'Visible','off');

bottom_text(txtb,'pwd',1);

keyboard



% Plot age for Gr and Labr shelves
% separately
figure(13); clf;
axes('Position',[0.2 0.45 0.4 0.4]);
hold on;
for ilr=1:nlr
% Plot Median and percentiles
  yp0=1.12*max(PRCT(:,3));
  pg1=PRCTG(ilr,1);
  pg2=PRCTG(ilr,3);
  pg4=PRCTG(ilr,4);
  pg5=PRCTG(ilr,5);
  mgd=PRCTG(ilr,2);

  pl1=PRCTL(ilr,1);
  pl2=PRCTL(ilr,3);
  pl4=PRCTL(ilr,4);
  pl5=PRCTL(ilr,5);
  mld=PRCTL(ilr,2);
  
% GrSh particles
  dl=0.1;
  i1=ilr-2*dl;
  i2=ilr-dl;
  i00=0.5*(i1+i2);
  patch([i1 i1 i2 i2],[pg1 pg2 pg2 pg1],clr1);
  plot([i1 i2],[mgd mgd],'Linewidth',2.6,'Color',clr11);
  plot([i00 i00],[pg1 pg4],'Linewidth',2.,'Color',clr1);
  plot([i00 i00],[pg2 pg5],'Linewidth',2.,'Color',clr1);

% LabrSh particles
  dl=0.15;
  i1=ilr+2*dl;
  i2=ilr+dl;
  i00=0.5*(i1+i2);
  patch([i1 i1 i2 i2],[pl1 pl2 pl2 pl1],clr2);
  plot([i1 i2],[mld mld],'Linewidth',2.6,'Color',clr22);
  plot([i00 i00],[pl1 pl4],'Linewidth',2.,'Color',clr2);
  plot([i00 i00],[pl2 pl5],'Linewidth',2.,'Color',clr2);

  plot([ilr+0.5 ilr+0.5],[0 32],':','Color',[0.8 0.8 0.8]);
end
set(gca,'tickdir','out',...
        'xlim',[0.5 4.6],...
        'ylim',[-1 28],...
        'xtick',[1:4],...
        'ytick',[0:2:30],...
        'ygrid','on',...
        'Fontsize',14);
title('GrSh (blue) and NLabrSh (red) particles, 25/75th prc Age Prtcl');
xlabel('Layers');
ylabel('Age, Years');
bottom_text(txtb,'pwd',1,'position',[0.02 0.2 0.4 0.04]);

%
% Plot mean and 95CI
figure(14); clf;
axes('Position',[0.2 0.45 0.4 0.4]);
hold on;
for ilr=1:nlr
  xmn = PRBG.Age_mn(ilr);
  ci1 = PRBG.CI95(ilr,1);
  ci2 = PRBG.CI95(ilr,2);

% GrSh particles
  dl=0.1;
  i1=ilr-dl;
  plot([i1 i1],[ci1 ci2],'-','Linewidth',2.2,'Color',clr1);
  plot(i1,xmn,'.','Markersize',20,'Color',clr11);


% LabrSh particles
  xmn = PRBL.Age_mn(ilr);
  ci1 = PRBL.CI95(ilr,1);
  ci2 = PRBL.CI95(ilr,2);

  i1=ilr+dl;
  plot([i1 i1],[ci1 ci2],'-','Linewidth',2.2,'Color',clr2);
  plot(i1,xmn,'.','Markersize',20,'Color',clr22);

  plot([ilr+0.5 ilr+0.5],[0 32],':','Color',[0.8 0.8 0.8]);
end
set(gca,'tickdir','out',...
        'xlim',[0.5 4.6],...
        'ylim',[-1 18],...
        'xtick',[1:4],...
        'ytick',[0:2:30],...
        'ygrid','on',...
        'Fontsize',14);
title('Mean Age GrSh(blue) and NLabrSh (red), 95CI');
xlabel('Layers');
ylabel('Age, Years');

axes('Position',[0.62 0.28 0.35 0.55]);
text(0.05,1,'Mean, 95CI','Fontsize',12);
for ilr=1:nlr
  i1=0.05;
  j1=1-ilr/10;

  xmn = PRBG.Age_mn(ilr);
  ci1 = PRBG.CI95(ilr,1);
  ci2 = PRBG.CI95(ilr,2);
  text(i1,j1,sprintf('GrSh Lr %i: %4.1f [%4.1f, %4.1f]',ilr,xmn,ci1,ci2),'Color',clr11,'Fontsize',12);

  j1=0.5-ilr/10;
  xmn = PRBL.Age_mn(ilr);
  ci1 = PRBL.CI95(ilr,1);
  ci2 = PRBL.CI95(ilr,2);
  text(i1,j1,sprintf('LarbSh Lr %i: %4.1f [%4.1f, %4.1f]',ilr,xmn,ci1,ci2),'Color',clr22,'Fontsize',12);

end
set(gca,'visible','off');

bottom_text(txtb,'pwd',1,'position',[0.02 0.2 0.4 0.04]);



% Summarize in boxplot:
% Overall stat for all particles pooled together
% Plot boxplots 
figure(11); clf;
axes('Position',[0.1 0.5 0.4 0.4]);
hold on;
for ilr=1:nlr
% Plot Median and percentiles
  yp0=1.12*max(PRCT(:,3));
  p1=PRCT(ilr,1);
  p2=PRCT(ilr,3);
  md=PRCT(ilr,2);

  clr0=[1 0 0];
  plot([ilr ilr],[p1 p2],'Linewidth',2,'Color',[0.6 0.6 0.6]);
  plot(ilr,md,'.','Markersize',30,'Color',[0 0 0]);

end
set(gca,'tickdir','out',...
        'xlim',[0.5 4.5],...
        'ylim',[0 17],...
        'xtick',[1:4],...
        'ytick',[0:2:30],...
        'ygrid','on',...
        'Fontsize',14);
title('All Particles Median, 25/75th prc Age Prtcl');
xlabel('Layers');
ylabel('Years');
bottom_text(txtb,'pwd',1,'position',[0.02 0.2 0.4 0.04]);




%
% Calculate average Tmean_age as weighted average (0.38*Tz(Labr)+0.62*Tz(GreenlSh))
% Tz is depth-average T mean_age
Zlr = [0; -50; -90; -150; -450; -1000];
dZ = abs(diff(Zlr));
Tmn = PRBL.Age_mn;
Tlow = PRBL.CI95(:,1);
Tup = PRBL.CI95(:,2);
TzLabr = Tmn(1)*dZ(1);
TzLabr1 = Tlow(1)*dZ(1);
TzLabr2 = Tup(1)*dZ(1);
for ik=2:length(Tmn);
  TzLabr  = TzLabr+0.5*(Tmn(ik-1)+Tmn(ik))*dZ(ik);
  TzLabr1 = TzLabr1+0.5*(Tlow(ik-1)+Tlow(ik))*dZ(ik);
  TzLabr2 = TzLabr2+0.5*(Tup(ik-1)+Tup(ik))*dZ(ik);
end
TzLabr = TzLabr+Tmn(end)*dZ(end);
TzLabr = TzLabr/abs(Zlr(end));
TzLabr1 = TzLabr1+Tlow(end)*dZ(end);
TzLabr1 = TzLabr1/abs(Zlr(end));
TzLabr2 = TzLabr2+Tup(end)*dZ(end);
TzLabr2 = TzLabr2/abs(Zlr(end));




Tmn = PRBG.Age_mn;
Tlow = PRBG.CI95(:,1);
Tup = PRBG.CI95(:,2);
TzGrsh = Tmn(1)*dZ(1);
TzGrsh1 = Tlow(1)*dZ(1);
TzGrsh2 = Tup(1)*dZ(1);
for ik=2:length(Tmn);
  TzGrsh  = TzGrsh+0.5*(Tmn(ik-1)+Tmn(ik))*dZ(ik);
  TzGrsh1 = TzGrsh1+0.5*(Tlow(ik-1)+Tlow(ik))*dZ(ik);
  TzGrsh2 = TzGrsh2+0.5*(Tup(ik-1)+Tup(ik))*dZ(ik);
end
TzGrsh = TzGrsh+Tmn(end)*dZ(end);
TzGrsh = TzGrsh/abs(Zlr(end));
TzGrsh1 = TzGrsh1+Tlow(end)*dZ(end);
TzGrsh1 = TzGrsh1/abs(Zlr(end));
TzGrsh2 = TzGrsh2+Tup(ik)*dZ(end);
TzGrsh2 = TzGrsh2/abs(Zlr(end));

ww1=0.389;
ww2=1-ww1;
Tmn_age = ww1*TzLabr+ww2*TzGrsh;
Tmn_low = ww1*TzLabr1+ww2*TzGrsh1;
Tmn_up  = ww1*TzLabr2+ww2*TzGrsh2;


fprintf('Depth-integrated over %i m Tmean_age = %6.4g yrs\n',abs(Zlr(end)),Tmn_age);
fprintf('             low/upper 95CI: %6.4g, %6.4g\n',Tmn_low,Tmn_up);













