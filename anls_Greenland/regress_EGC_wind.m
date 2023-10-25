% HYCOM 0.08
%
% For wind roses (map boxes, etc.) 
% see ../anls_atm/calc_Ekman_trpt.m
%
% Regress along-isobath wind component
% on EGCC: hb (bottom depth of EGCC front)
% width, mean current (depth-integrated)
%
%Something similar to Sutherland and Pickard, JPO
% 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

irs=5; % select SE Gr Shelf point wind rose that matches with SE section

Sc0=33.8;  % Salinity that defines the front
%Sc0=33.5;  % Salinity that defines the front
%Sc0=33.6;  % Salinity that defines the front
%Sc0=33.7;  % Salinity that defines the front

%regn='ARCc0.08';
%expt=112;
regn='ARCc0.04';
expt=012;


rhoa = 1.2;
Cd   = 0.0013;
rhow = 1027;

%s_mat = 1;

pthdat = '/Net/data/ccmp/v02.0/';
pthdt  = '/Net/tholia/ddmitry/ovwst/data_mat/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pth8='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthdat1 = '/Net/kronos/ddmitry/ncep_cfsr/';
%pthdat2 = '/Net/kronos/ddmitry/ncep_cfsv2/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';
pthmat = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/data_GrSect/',regn,expt);

%pthout  = ('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/');
%fmat    = sprintf('%sarc008_Greenl_upwl_CFSR_month.mat',pthout);
fmat = sprintf('%swrose_Greenland_ccmp020rss_1997_2016.mat',pthdt);
btx='regress_EGC_wind.m';

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
fprintf('Getting bathymetry %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

fprintf('Loading %s\n',fmat);
load(fmat); % WCLIM WRS

nbx=length(WRS);

% Find general directions 
% of the coast for each box
% by specifing endponts of
% the vector 
% c/clckwise direction
% around Greenland such that
% positive wind projection = upwelling
% Indices are for 0.08 grid! 
% Only needed to find North direction
% at Wind roses sites
% ok to use with hycom 0.04 output
%
Dr(1).XY=[924   787
   949   961];
Dr(2).XY=[882   645
   922   796];
Dr(3).XY=[753   559
   878   645];
Dr(4).XY=[656   500
   789   586];
Dr(5).XY=[617   415
   651   497];
Dr(6).XY=[581   424
	 619   415];
Dr(7).XY=[535   537
	  565   440];
Dr(8).XY=[535   660
	  538   518];
Dr(9).XY=[582   780
	  525   638];
Dr(10).XY=[609   882
	   570   756];
Dr(11).XY=[581   957
	   614   914];

% Calculate orientation of the
% local isobath directions
% wrt to 0 degree = Eastward axis
for ik=1:length(Dr)
  i1=Dr(ik).XY(1,1);
  i2=Dr(ik).XY(2,1);
  j1=Dr(ik).XY(1,2);
  j2=Dr(ik).XY(2,2);
  ln1=LON(j1,i1);
  ln2=LON(j2,i2);
  lt1=LAT(j1,i1);
  lt2=LAT(j2,i2);
%  L=distance_spheric_coord(lt2,ln2,lt1,ln1);
  Lx=dsphc(lt1,ln2,lt1,ln1);
  Lx=sign(ln2-ln1)*Lx;
  Ly=dsphc(lt2,ln1,lt1,ln1);
  Ly=sign(lt2-lt1)*Ly;
  
 % plot([i1 i2],[j1 j2],'-');
%
% Vector in 0 direction - East:
%  Nx=1;
%  Ny=0;
  thtd=atan2d(Ly,Lx);
  Dr(ik).Theta_dgr=thtd;
end;

% Regression for SE Gr. Shelf
ik=irs;
U=WRS(ik).U;
V=WRS(ik).V;
S = sqrt(U.^2+V.^2);

% Wind-drag - Large and Pond 1981
Cd=U*0+1.2e-3; % |U|<11 m/s
I=find(S>=11 & S<25);
Cd(I)=(0.49+0.065*S(I))*1e-3;
I=find(S>=25);
Cd(I)=2.115*1e-3;

tx = rhoa*Cd.*S.*U;
ty = rhoa*Cd.*S.*V;

% Upwelling favorable winds are negative
% Project winds on the isobath:
thtD=Dr(ik).Theta_dgr;
Pcst=[cosd(thtD), sind(thtD)]; % unit vector in dir of isobath
Tproj=[tx,ty]*Pcst'; % Stress Projection on the isobath
Wproj=[U,V]*Pcst'; % Wind vector Projection on the isobath
  
TM=WRS(1).TM;
DV=datevec(TM);

% -----------------------------
% HYCOM currents:
% Get EGC U,T,S section
% -----------------------------
switch (regn)
 case('ARCc0.04');
%  foutp=sprintf('%sarc04_%3.3i_anls_EGC.mat',pthmat,expt);
  foutp=sprintf('%sarc04_%3.3i_anls_EGC%3.3i.mat',pthmat,expt,Sc0*10);
 case('ARCc0.08');
  foutp=sprintf('%sarc08_%3.3i_anls_EGC.mat',pthmat,expt);
  foutp=sprintf('%sarc08_%3.3i_anls_EGC%3.3i.mat',pthmat,expt,Sc0*10);
end

fprintf('Loading %s\n',foutp);
load(foutp);

Hb=PLM.Hbottom;
Ish=max(find(Hb>-220));
TMsh=PLM.TM;
Hp=PLM.Front_foot_depth; % depth of EGCC front foot
Wbtm=PLM.DistBtm; % width of the bottom part
Wsrf=PLM.DistSrf; % width surface part
Sbanm=nanmean(PLM.Sbtm_anom(:,1:Ish),2); % S bottom anomaly
Umn=PLM.Umean_plume;
day1=TMsh(1);
day2=TMsh(end);
Wsrf=Wsrf(:);
Wbtm=Wbtm(:);
Umn=Umn(:);
Hp=Hp(:);



% Create n-day averaged wind stress for regression
nav=3; % # of previous days + the current day
cc=0;
clear Tau
for ic=1:length(TMsh)
  iday=TMsh(ic);
  iTm1=find(TM==iday-nav); % 2 previous days
  iTm2=find(TM==iday+1)-1; % include the current day

%  w10=Wproj(iTm1:iTm2); % wind vector
  w10=Tproj(iTm1:iTm2); % stress
  aw10=abs(w10);
  mw10=max(aw10);
  ii=find(aw10==mw10);
  mxw10=sign(w10(ii))*mw10;
  avw10=mean(w10);
  
% dS = bottom-surface
  dmm=PLM.dltSvrt;
  dS=nanmean(dmm,2);
  
  cc=cc+1;
%  Tau(cc,1)=mxw10;
  Tau(cc,1)=avw10;
end

% Separate by seasons:

DV=datevec(TMsh);
Is=find(DV(:,2)>=6 & DV(:,2)<=8); % summer
Iw=find(DV(:,2)<=3 | DV(:,2)>=10); % fall


POS=[0.1 0.6 0.3 0.3;...
     0.6 0.6 0.3 0.3;...
     0.1 0.2 0.3 0.3;...
     0.6 0.2 0.3 0.3];

xl='Tau, Pa';  
f_all=0;
if f_all==1
  figure(1); clf;
  ps=POS(1,:);
  axes('Position',ps);
  plot(Tau,Umn,'.');
  xlabel(xl);
  ylabel('Umean');

  ps=POS(2,:);
  axes('Position',ps);
  plot(Tau,Hp,'.');
  xlabel(xl);
  ylabel('Hp, m');

  ps=POS(3,:);
  axes('Position',ps);
  plot(Tau,log(Wsrf),'.');
  xlabel(xl);
  ylabel('log(Wsrf), m');

  ps=POS(4,:);
  axes('Position',ps);
  plot(Tau,log(Wbtm),'.');
  xlabel(xl);
  ylabel('log(Wbtm), m');
end

% Winter
stl=sprintf('%s-%3.3i, Winter, %i-%i',regn,expt,DV(1,1),DV(end,1));
figure(2); clf;
ps=POS(1,:);
axes('Position',ps);
plot(Tau(Iw),Umn(Iw),'.');
xlabel(xl);
ylabel('Umean');
title(stl);

ps=POS(2,:);
axes('Position',ps);
%plot(Tau(Iw),Hp(Iw),'.');
%xlabel(xl);
%ylabel('Hp, m');
plot(Tau(Iw),dS(Iw),'.');
ylabel('dlt S');

ps=POS(3,:);
axes('Position',ps);
plot(Tau(Iw),log(Wsrf(Iw)),'.');
xlabel(xl);
ylabel('log(Wsrf), m');

ps=POS(4,:);
axes('Position',ps);
%plot(Tau(Iw),log(Wbtm(Iw)),'.');
%xlabel(xl);
%ylabel('log(Wbtm), m');
plot(Tau(Iw),Sbanm(Iw),'.');
xlabel(xl);
ylabel('Sbtm anom');


sinf{1}=sprintf('Wind avrg = %i days',nav);
sinf{2}=sprintf('Front def S=%4.2f',Sc0);

axes('Position',[0.5 0.03 0.45 0.05]);
text(0,0.5,sinf);
set(gca,'visible','off');
bottom_text(btx,'pwd',1);


% Summer
stl=sprintf('%s-%3.3i, Summer, %i-%i',regn,expt,DV(1,1),DV(end,1));
figure(3); clf;
ps=POS(1,:);
axes('Position',ps);
plot(Tau(Is),Umn(Is),'.');
xlabel(xl);
ylabel('Umean');
title(stl);

ps=POS(2,:);
axes('Position',ps);
%plot(Tau(Is),Hp(Is),'.');
xlabel(xl);
%ylabel('Hp, m');
plot(Tau(Is),dS(Is),'.');
ylabel('dlt S');

ps=POS(3,:);
axes('Position',ps);
plot(Tau(Is),log(Wsrf(Is)),'.');
xlabel(xl);
ylabel('log(Wsrf), m');

ps=POS(4,:);
axes('Position',ps);
%plot(Tau(Is),log(Wbtm(Is)),'.');
%xlabel(xl);
%ylabel('log(Wbtm), m');
plot(Tau(Is),Sbanm(Is),'.');
xlabel(xl);
ylabel('Sbtm anom');

axes('Position',[0.5 0.03 0.45 0.05]);
text(0,0.5,sinf);
set(gca,'visible','off');

bottom_text(btx,'pwd',1);



% Winter tau-U mean
xx=Tau(Iw);
xx1=xx(1:end-1);
xx2=xx(2:end);
xx1=[xx1;nan];
xx2=[xx2;nan];
Xb=[ones(length(xx),1),xx1,xx2];
X=[ones(length(xx),1),xx];
Y=Umn(Iw);
[b,bint,r,rint,stats]=regress(Y,X);
uwf=X*b;
Bw_tu=b;

% Summer tau-U mean
xx=Tau(Is);
X=[ones(length(xx),1),xx];
Y=Umn(Is);
[b,bint,r,rint,stats]=regress(Y,X);
Bs_tu=b;

% Winter tau-log(Wsrf) - width surface front
xx=Tau(Iw);
X=[ones(length(xx),1),xx];
Y=log(Wsrf(Iw));
[b,bint,r,rint,stats]=regress(Y,X);
Bw_ts=b;

% Summer tau-log(Wsrf) - width surface front
xx=Tau(Is);
X=[ones(length(xx),1),xx];
Y=log(Wsrf(Is));
[b,bint,r,rint,stats]=regress(Y,X);
Bs_ts=b;

% Winter tau-Sbtm anom 
xx=Tau(Iw);
xx1=xx(1:end-1);
xx2=xx(2:end);
xx1=[xx1;nan];
xx2=[xx2;nan];
Xb=[ones(length(xx),1),xx1,xx2];
X=[ones(length(xx),1),xx];
Y=Sbanm(Iw);
[b,bint,r,rint,stats]=regress(Y,X);
[b,bint,r,rint,stats]=regress(Y,Xb);
Bw_ts=b;

% Summer tau- Sbtm anom
xx=Tau(Is);
X=[ones(length(xx),1),xx];
Y=Sbanm(Is);
[b,bint,r,rint,stats]=regress(Y,X);
Bs_ts=b;


