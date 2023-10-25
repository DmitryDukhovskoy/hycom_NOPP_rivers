% Analysis of particles - histograms, mean transit time etc.
%
%  Gr Shelf Particles 
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
addpath /Net/Movies0/ddmitry/MyMatlab
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
pthfig = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_GrGprtcl/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
txtb = 'anls_AllGrSh_prtcls.m';

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
 
end


YRPLT=[1993:2019];
nplt = length(YRPLT);

fnmout=sprintf('%sGrShprt_comb.mat',pthmat);

if fget>0
  % Combine all trajectories:
  ik=0;
  for ip=1:nlr
    for it=1:nnm
      ik=ik+1;
      PP(ik).Xp=[];
      PP(ik).Yp=[];
      PP(ik).TM(1)=datenum(1993,1,1); % missing initial date

    end
  end


  for ii=1:nplt
    yr = YRPLT(ii);
  %  fmat = sprintf('%sGG_particles-lr%2.2i_%i.mat',pthmat,slr,yr);
  % Combine all simulations
    icc=0;  % all nlevels*nexpts
    for ilr=1:nlr
      slr=SLR(ilr);
      for iip=1:nnm
	       nsim=NSM(iip);
       	fmat = sprintf('%sGrSh_part-lr%2.2i_%i-%2.2i.mat',pthmat,slr,yr,nsim);
	
								if ~exist(fmat,'file')
										fprintf('Not found %s\n',fmat);
										continue;
								end
								
								fprintf('Loading saved %s\n',fmat);
								load(fmat);


								TR = PRTCL.TRACK;
								nr = length(TR);

								icc=icc+1;
								
								Xp=PP(icc).Xp;
								Yp=PP(icc).Yp;
								TM=PP(icc).TM;
								for it=1:nr    
										X = TR(it).I;
										Y = TR(it).J;
										Xp=[Xp,X];
										Yp=[Yp,Y];
										TM=[TM;TR(it).TM];
								end
								PP(icc).Xp=Xp;
								PP(icc).Yp=Yp;
								PP(icc).TM=TM;
								PP(icc).layer=slr;
								PP(icc).nsim=nsim;
      end  % simulations in 1 layer
    end  % layers
  end

  fprintf('Saving %s\n',fnmout);
  save(fnmout,'PP');
  
else
  fprintf('Loading %s\n',fnmout);
  load(fnmout);
end


% Combine all experiments with Lagrangian particles together 
% for missing dates - use the previous number and locations of particles
nll=length(PP);
ntm=0;
for iip=1:nll
  tm=PP(iip).TM;
  nt=length(tm);
  ntm=max([ntm,nt]);
  if ntm==nt
    Li=iip;
  end
end

TM0=PP(Li).TM;
% Exclude repeated dates:
dTm=diff(TM0);
I=find(dTm>0);
I=[I;ntm];
TM0=TM0(I);
ntm=length(TM0);

% Fill in missing dates in the experiments 
% so that all experiments have same # of days (ntm)
clear nnm
nll=length(PP);
for iip=1:nll
  Xp=PP(iip).Xp;
  Yp=PP(iip).Yp;
  Tm=PP(iip).TM;
  zl=PP(iip).layer;
  [a1,a2]=size(Xp);

  Xn=[];
  Yn=[];
  for it=1:ntm
    t0=TM0(it);
    jt=find(Tm==t0,1);

    if ~isempty(jt),
      Xn(:,it)=Xp(:,jt);
      X0=Xp(:,jt);
      Yn(:,it)=Yp(:,jt);
      Y0=Yp(:,jt);
    else
      Xn(:,it)=X0;
      Yn(:,it)=Y0;
    end
  end

  PP(iip).Xp=Xn;
  PP(iip).Yp=Yn;
  PP(iip).ZL=ones(a1,1)*zl;
end


% Combine all experiments      
XP=[];
YP=[];
ZL=[]; 
TM=TM0;

for iip=1:nll
  Xp=PP(iip).Xp;
  Yp=PP(iip).Yp;
  zl=PP(iip).ZL;
  XP=[XP;Xp];
  YP=[YP;Yp];
  ZL=[ZL;zl];
end;




[a1,a2]=size(XP);
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


%INP = inpolygon(XX,YY,IGR(:,1),IGR(:,2));
%INGr= find(INP==1 & HH<0);
%INdeep = find(INP==1 & HH<-800);
%II = find(HH<0);



TMd=TM-TM(1);
TMy=TMd./365.24;

clear INdp


% Calculate probability of part. being inside SPNA at least for N years after release
NDout=90;
NDays=[];
for ilr=1:nlr
  lr=SLR(ilr);
  IZ=find(ZL==lr);
  npp=length(IZ);
  dTMy=diff(TMy);
  dTMy(end+1)=dTMy(end);

  XPz=XP(IZ,:);
  YPz=YP(IZ,:);
  [a1,a2]=size(XPz);
  for ipp=1:npp
    X0=XPz(ipp,:);
    Y0=YPz(ipp,:);
    II=inpolygon(X0,Y0,IGR(:,1),IGR(:,2));
    dtm=dTMy.*II';

    i0=min(find(II==0));
    if isempty(i0), i0=a2; end;
    while i0<a2
      dmm=II;
      dmm(1:i0)=0;
      i2=min(find(dmm>0));
      ddays=sum(dTMy(i0:i2-1))*365;
      if isempty(i2) | ddays>NDout; break; end;
      dmm=II;
      dmm(1:i2)=1;
      i0=min(find(dmm==0));
      if isempty(i0), i0=a2; end;
    end
    NDays(ilr,ipp) = sum(dtm(1:i0));

  end
end
 
% Calculate probability of partlce staying within the SPNA for at least N years
% Probability of transit time (exclude points that never left the domain)
Iyr = [1:27];
nnr = length(Iyr);
NPyr=[];
TRyr=[];
Dmax = sum(dTMy); % max possible transit time for particle, otherwise it staied in the domain
for jj=1:nnr
  yr0=Iyr(jj);
  for ilr=1:nlr
    aa=NDays(ilr,:);
    IP = find(aa>=yr0); 
    NPyr(ilr,jj)=length(IP);
    IT = find(aa>=yr0 & aa<=Dmax);
    TRyr(ilr,jj)=length(IT);
  end
end    

% Probability:
PRyr = NPyr./npp;
TRyr = TRyr./npp;

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

figure(10); clf;
set(gcf,'Position',[854 503 1687 828]);
for ilr=1:nlr
%  clr=CLRZ(ilr,:);
  pos=POS(ilr,:);
  axes('Position',pos);
  hold on;
  br = PRyr(ilr,:);
  hb = bar(Iyr,br);
  set(hb,'Facecolor',[0.6 0.6 0.6]);
% 
% Plot P(Age<=Nyrs):
% Get probability of Age <= N years:
  Page = 1-PRyr(ilr,:);
  Page = Page(:);
  Page = [0; Page];
  Pi = interp1(T0,Page,Ti,'pchip');
  plot(Ti,Pi,'-','Color',[0 0 0],'linewidth',2); 

% Save median and percentiles:
  i1=max(find(Pi<=0.1)); % 10%
  PRCT(ilr,1)=Ti(i1);
  i1=max(find(Pi<=0.5));
  PRCT(ilr,2)=Ti(i1);
  i1=max(find(Pi<=0.9));
  PRCT(ilr,3)=Ti(i1);


  stl = sprintf('Bar: Prob(Age>=N), Line: P(Age<=N), Lr %i, Z=%i',ilr,ZLR(ilr));
  title(stl);

  set(gca,'Tickdir','out',...
   'Fontsize',14,...
   'xlim',[0 28],...
   'xtick',Iyr,...
   'ylim',[0 1],...
   'ytick',[0:0.1:1],...
   'xgrid','on',...
   'ygrid','on');

end

bottom_text(txtb,'pwd',1);

% Summarize in boxplot:
% Plot Median and percentiles:
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
        'ylim',[0 30],...
        'xtick',[1:4],...
        'ytick',[0:5:30],...
        'ygrid','on',...
        'Fontsize',14);
title('Median, 10/90th prc Age Prtcl');
xlabel('Layers');
ylabel('Years');
bottom_text(txtb,'pwd',1,'position',[0.08 0.34 0.4 0.04]);


% =====================
% 
%  Transit time:
%  get rid of particles that never 
%  leave the domain
% =====================
figure(12); clf;
set(gcf,'Position',[854 503 1687 828]);
for ilr=1:nlr
%  clr=CLRZ(ilr,:);
  pos=POS(ilr,:);
  axes('Position',pos);
  hold on;
  br = TRyr(ilr,:);
  hb = bar(Iyr,br);
  set(hb,'Facecolor',[0.6 0.6 0.6]);
% 
% Plot P(Age<=Nyrs):
% Get probability of Age <= N years:
  Page = 1-TRyr(ilr,:);
  Page = Page(:);
  Page = [0; Page];
  Pi = interp1(T0,Page,Ti,'pchip');
  plot(Ti,Pi,'-','Color',[0 0 0],'linewidth',2); 

% Save median and percentiles:
  i1=max(find(Pi<=0.1)); % 10%
  PRTR(ilr,1)=Ti(i1);
  i1=max(find(Pi<=0.5));
  PRTR(ilr,2)=Ti(i1);
  i1=max(find(Pi<=0.9));
  PRTR(ilr,3)=Ti(i1);


  stl = sprintf('Bar: Prob(TransT>=N), Line: P(TransT<=N), Lr %i, Z=%i',ilr,ZLR(ilr));
  title(stl);

  set(gca,'Tickdir','out',...
   'Fontsize',14,...
   'xlim',[0 28],...
   'xtick',Iyr,...
   'ylim',[0 1],...
   'ytick',[0:0.1:1],...
   'xgrid','on',...
   'ygrid','on');

end

bottom_text(txtb,'pwd',1);

% Summarize in boxplot:
% Plot Median and percentiles:
figure(13); clf;
axes('Position',[0.1 0.5 0.4 0.4]);
hold on;
for ilr=1:nlr
% Plot Median and percentiles
  yp0=1.12*max(PRCT(:,3));
  p1=PRTR(ilr,1);
  p2=PRTR(ilr,3);
  md=PRTR(ilr,2);

  clr0=[1 0 0];
  plot([ilr ilr],[p1 p2],'Linewidth',2,'Color',[0.6 0.6 0.6]);
  plot(ilr,md,'.','Markersize',30,'Color',[0 0 0]);

end
set(gca,'tickdir','out',...
        'xlim',[0.5 4.5],...
        'ylim',[0 30],...
        'xtick',[1:4],...
        'ytick',[0:5:30],...
        'ygrid','on',...
        'Fontsize',14);
title('Median, 10/90th prc Transit T');
xlabel('Layers');
ylabel('Years');
bottom_text(txtb,'pwd',1,'position',[0.08 0.34 0.4 0.04]);



keyboard
% Below is obsolete analaysis - histograms etc.


% 
% Histogram of time inside the SPG by each particles
% Plot total time and plot time until the particle leaves SPG
% Time is considered to leave SPG if it stays > NDout
% Estimate mean age = average age of parcels in domain at time T
%
% Transit time = average age of water parcles leaving the domain
%
% Residnce time = Transit time or mean age
% For steady state should be same
%
NDout=90;
clear TMprt TMpD
for ilr=1:nlr
  fprintf('Calculating mean part age, ilr=%i\n',ilr);
% Mean Age
  lr=SLR(ilr);
  IZ=find(ZL==lr);
  npp=length(IZ);
  dTMy=diff(TMy);
  dTMy(end+1)=dTMy(end);
  
  XPz=XP(IZ,:);
  YPz=YP(IZ,:);  
  [a1,a2]=size(XPz);
  for ipp=1:npp
    X0=XPz(ipp,:);
    Y0=YPz(ipp,:);
    II=inpolygon(X0,Y0,IGR(:,1),IGR(:,2));
    dtm=dTMy.*II';
    TMprt(ilr,ipp)=nansum(dtm); % all days in SPG - not very useful 
    i0=min(find(II==0));
    if isempty(i0), i0=a2; end;
    while i0<a2
      dmm=II;
      dmm(1:i0)=0;
      i2=min(find(dmm>0));
      ddays=sum(dTMy(i0:i2-1))*365;
      if isempty(i2) | ddays>NDout; break; end;
      dmm=II;
      dmm(1:i2)=1;
      i0=min(find(dmm==0));
      if isempty(i0), i0=a2; end;
    end
    TMpD(ilr,ipp)=sum(dtm(1:i0));   % Mean Age
 
  end
end

clear a1 Nhst Dhst
ddx=2;
xb=[0:ddx:28];
for ilr=1:nlr
  lr=SLR(ilr);
  IZ=find(ZL==lr);
  npp=length(IZ);
%  NN=hist(TMprt(ilr,:),xb);
  NN=histcounts(TMprt(ilr,:),xb);
  NN=NN/npp; % fraction of all particles
  Nhst(ilr,:)=NN;
% Years before particles leave - no return
%  LD=hist(TMpD(ilr,:),xb);
  LD=histcounts(TMpD(ilr,:),xb);
  LD=LD/npp;
  Dhst(ilr,:)=LD;
  
  ST(ilr).Layer=SLR(ilr);
  ST(ilr).ZL=ZLR(ilr);
  ST(ilr).Mean=mean(TMpD(ilr,:));
  ST(ilr).Median=median(TMpD(ilr,:));
  ST(ilr).p10=prctile(TMpD(ilr,:),10);
  ST(ilr).p90=prctile(TMpD(ilr,:),90);
  ST(ilr).AllMean=mean(TMprt(ilr,:));
  ST(ilr).AllMedian=median(TMprt(ilr,:));
  ST(ilr).Allp10=prctile(TMprt(ilr,:),10);
  ST(ilr).Allp90=prctile(TMprt(ilr,:),90);
end
llb=length(xb)-1;
dlb=ddx/nlr*0.9;
dfst=(ddx-dlb*nlr)/2;

figure(3); clf;
axes('Position',[0.08 0.5 0.9 0.4]);
hold on;
for ilr=1:nlr
  clr=CLRZ(ilr,:);
  for ik=1:llb
    yy=Nhst(ilr,ik);
    x0=xb(ik)+dfst+dlb*(ilr-1);
    xv=[x0,x0,x0+dlb,x0+dlb];
    yv=[0,yy,yy,0];
    patch(xv,yv,clr,'edgecolor','none');
  end
end
  
%hbr=bar(xb+ddx/2,Nhst',1.,'grouped');
%colormap(CLRZ);

title('Total Time Lagr. prt in SPG (return allowed)');
set(gca,'Tickdir','out',...
	'Fontsize',14,...
	'xlim',[0 xb(end)],...
	'xtick',xb,...
	'ylim',[0 1.1*max(max(Nhst))],...
	'ytick',[0:0.05:1.2],...
	'xgrid','on',...
	'ygrid','on');

xlabel('Years');
ylabel('Probability');

axes('Position',[0.6 0.25 0.3 0.2]);
hold

for ilr=1:nlr
  clr=CLRZ(ilr,:);
  zlr=ZLR(ilr);
  dy=1.2-ilr/nlr;
  plot([0.2 0.23],[dy dy],'-','Linewidth',20,'Color',clr);
  text(0.25,dy,sprintf('%i m',zlr),'Fontsize',14);
end
set(gca,'xlim',[0 0.4],...
	'ylim',[0 1.5],...
	'visible','off')

% Stat:
axes('Position',[0.1 0.1 0.3 0.25]);
hold
for ilr=1:nlr
  clr=CLRZ(ilr,:);
  p1=ST(ilr).Allp10;
  p2=ST(ilr).Allp90;
  mn=ST(ilr).AllMean;
  md=ST(ilr).AllMedian;
  plot(ilr,md,'.','Markersize',40,'Color',clr);
  plot([ilr ilr],[p1 p2],'Linewidth',2,'Color',clr);
end
set(gca,'tickdir','out',...
        'xlim',[0.5 nlr+0.5],...
        'ylim',[0 30],...
        'xtick',[1:nlr],...
        'ytick',[0:5:30],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
title('10, 50, 90th prctiles');

bottom_text(txtb,'pwd',1,'Position',[0.4 0.06 0.4 0.04]);


figure(4); clf;
axes('Position',[0.08 0.5 0.9 0.4]);
hold on;
for ilr=1:nlr
  clr=CLRZ(ilr,:);
  for ik=1:llb
    yy=Dhst(ilr,ik);
    x0=xb(ik)+dfst+dlb*(ilr-1);
    xv=[x0,x0,x0+dlb,x0+dlb];
    yv=[0,yy,yy,0];
    patch(xv,yv,clr,'edgecolor','none');
  end
end
  
%hbr=bar(xb+ddx/2,Nhst',1.,'grouped');
%colormap(CLRZ);

title('Age Part. (Time prt in SPG - no return >90 days)');
set(gca,'Tickdir','out',...
	'Fontsize',14,...
	'xlim',[0 xb(end)],...
	'xtick',xb,...
	'ylim',[0 1.1*max(max(Nhst))],...
	'ytick',[0:0.05:1.2],...
	'xgrid','on',...
	'ygrid','on');

xlabel('Years');
ylabel('Probability');

axes('Position',[0.6 0.25 0.3 0.2]);
hold

for ilr=1:nlr
  clr=CLRZ(ilr,:);
  zlr=ZLR(ilr);
  dy=1.2-ilr/nlr;
  plot([0.2 0.23],[dy dy],'-','Linewidth',20,'Color',clr);
  text(0.25,dy,sprintf('%i m',zlr),'Fontsize',14);
end
set(gca,'xlim',[0 0.4],...
	'ylim',[0 1.5],...
	'visible','off')

% Stat:
axes('Position',[0.1 0.1 0.3 0.25]);
hold
for ilr=1:nlr
  clr=CLRZ(ilr,:);
  p1=ST(ilr).p10;
  p2=ST(ilr).p90;
  mn=ST(ilr).Mean;
  md=ST(ilr).Median;
  plot(ilr,md,'.','Markersize',40,'Color',clr);
  plot([ilr ilr],[p1 p2],'Linewidth',2,'Color',clr);
end
set(gca,'tickdir','out',...
        'xlim',[0.5 nlr+0.5],...
        'ylim',[0 30],...
        'xtick',[1:nlr],...
        'ytick',[0:5:30],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
title('10, 50, 90th prctiles');

bottom_text(txtb,'pwd',1,'Position',[0.4 0.06 0.4 0.04]);

%
% Same plot differently presented
dpx=0.38;
dpy=0.32;
POS = [0.09 0.6 dpx dpy; ...
       0.59 0.6 dpx dpy; ...
       0.09 0.15 dpx dpy; ...
       0.59 0.15 dpx dpy];

figure(5); clf;
for ilr=1:nlr
%  clr=CLRZ(ilr,:);
  pos=POS(ilr,:);
  axes('Position',pos);
  hold on; 
 
%  HS = histogram(TMpD(ilr,:),xb,'Normalization','probability');
%  set(HS,'FaceColor',[0.6 0.6 0.6]);

  LD=histcounts(TMpD(ilr,:),xb);
  npp=length(find(TMpD(ilr,:)>0));
  LD=LD/npp;
  
  for ik=1:llb
    yy=LD(ik);
    x0=xb(ik);
    x1=xb(ik+1);
    xv=[x0,x0,x1,x1];
    yv=[0,yy,yy,0];
    clr = [0.6 0.6 0.6];
    patch(xv,yv,clr);
  end

  stl=sprintf('Pprobab Distrib. Part Age in the SPNA (no retrun), L=%i',SLR(ilr));
  title(stl);

		set(gca,'Tickdir','out',...
			'Fontsize',14,...
			'xlim',[-0.1 30],...
			'xtick',xb,...
   'ylim',[0 1.2*max(LD)],...
			'ytick',[0:0.05:1.2],...
			'xgrid','on',...
			'ygrid','on');

% Plot Median and percentiles
%  yp0=1.12*max(LD);
%  p1=ST(ilr).p10;
%  p2=ST(ilr).p90;
%  md=ST(ilr).Median;
  
%  clr0=[1 0 0];
%  plot(md,yp0,'.','Markersize',30,'Color',clr0);
%  plot([p1 p2],[yp0 yp0],'Linewidth',2,'Color',clr0);

  xlabel('Years');
  ylabel('Probability');
end

bottom_text(txtb,'pwd',1);

% Plot Median and percentiles:
figure(6); clf;
axes('Position',[0.1 0.5 0.4 0.4]);
hold on;
for ilr=1:nlr
% Plot Median and percentiles
  yp0=1.12*max(LD);
  p1=ST(ilr).p10;
  p2=ST(ilr).p90;
  md=ST(ilr).Median;

  clr0=[1 0 0];
  plot([ilr ilr],[p1 p2],'Linewidth',2,'Color',[0.6 0.6 0.6]);
  plot(ilr,md,'.','Markersize',30,'Color',[0 0 0]);

end
set(gca,'tickdir','out',...
        'xlim',[0.5 4.5],...
        'ylim',[0 30],...
        'xtick',[1:4],...
        'ytick',[0:5:30],...
        'ygrid','on',...
        'Fontsize',14);
title('Median, 10/90th prc Transit Time');
xlabel('Layers');
ylabel('Years');
bottom_text(txtb,'pwd',1,'position',[0.08 0.34 0.4 0.04]);









