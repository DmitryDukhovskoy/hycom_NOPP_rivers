% Plot fluxes 
% extracted data are in extr_TSVdaily_straits08.m
% The code uses improved code for collocating
% T,S, dh with U,V points along the segments
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=112;
TV=11;
YR1=2005;
YR2=2016;
dday=7;

pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='plot_TSVdaily_straits.m';

Tref1 = -1.8; % Ref T to calc. H flux
Tref2 = 0;    % Ref T to calc. H flux
Sref1 = 34.8; 
Sref2 = 34.9; 

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
%[DX,DY]=sub_dx_dy(LON,LAT);

% Fram Section is close to Moorings ~79N
%SCT = sub_define_sections(HH,LON,LAT);
%nsct = length(SCT);


CLR=[0 0.4 0.8; ...
     0.8 0.4 0; ...
     1 0.2 0; ...
     0 1 0; ...
     0.8 0 0.6; ...
     0.8 1 0];
%
% Average all years:
isct=5;
VFlx=[];
FWf1=[];
FWf2=[];
FWf3=[];
Hf1=[];
Hf2=[];
UU=[];
SS=[];
TT=[];
TM=[];
cc=0;
ctot=0;
for YR=YR1:YR2
  yr=YR;      
  dE=datenum(yr,12,31);
  dJ1=datenum(yr,1,1);
  ndays=dE-dJ1+1;
  tdays=dJ1+[1:dday:365]-1;  
  TM=[TM; tdays];

  fmatout=sprintf('%shycom008_%3.3i_StraitFluxesDay_%4.4i.mat',...
                    pthmat,expt,YR);
  fprintf('Loading %s\n',fmatout);
  load(fmatout);

  cc=cc+1;
  if ~exist('ZZi','var')
    ZZi=SCT(isct).ZZintrp;
    Slon=SCT(isct).long;
    Slat=SCT(isct).latd;
    dL=SCT(isct).segm_dL;
    IIs=SCT(isct).I;
    JJs=SCT(isct).J;
  end
  nm=SCT(isct).Name;
  AA=SCT(isct).Unrm;
  dmm=squeeze(mean(AA,1));
  UU(cc,:,:)=dmm;
 
  AA=SCT(isct).T;
  dmm=squeeze(mean(AA,1));
  TT(cc,:,:)=dmm;
 
  AA=SCT(isct).S;
  dmm=squeeze(mean(AA,1));
  SS(cc,:,:)=dmm;

%  VFlx=[VFlx;SCT(isct).VolFlx_m3s];
%  FWf1=[FWf1;SCT(isct).FWflx1_m3s];  % FW ref S1
%  FWf2=[FWf2;SCT(isct).FWflx2_m3s];  % FW ref S2
%  FWf3=[FWf3;SCT(isct).FWflx1_s1_m3s];  % FW ref S1 only S<=S1
%  Hf1=[Hf1;SCT(isct).Hflx1_W];  % Heat Flux ref T1
%  Hf2=[Hf2;SCT(isct).Hflx2_W];  % Heat Flux ref T2
  nrc=length(SCT(isct).VolFlx_m3s); 
  VFlx(cc,1:nrc)=SCT(isct).VolFlx_m3s;
  FWf1(cc,1:nrc)=SCT(isct).FWflx1_m3s;     % FW ref S1
  FWf2(cc,1:nrc)=SCT(isct).FWflx2_m3s;     % FW ref S2
  FWf3(cc,1:nrc)=SCT(isct).FWflx1_s1_m3s;  % FW ref S1 only S<=S1
  Hf1(cc,1:nrc)=SCT(isct).Hflx1_W;  % Heat Flux ref T1
  Hf2(cc,1:nrc)=SCT(isct).Hflx2_W;  % Heat Flux ref T2


end

% Overall mean:
UUm=squeeze(mean(UU,1));
TTm=squeeze(mean(TT,1));
SSm=squeeze(mean(SS,1));

% Get HYCOM bottom
ni=length(IIs);
for ii=1:ni
  i1=IIs(ii);
  j1=JJs(ii);

  Hb(ii)=HH(j1,i1);
end;






%keyboard

VFlx=VFlx*1e-6; % m3/s-> Sv
FWf1=FWf1*1e-3; % mSv
FWf2=FWf2*1e-3; % mSv
FWf3=FWf3*1e-3; % mSv
Hf1=Hf1*1e-12;  % TW
Hf2=Hf2*1e-12;  % TW


% Time series: 
% monthly climatology
[a1,a2]=size(TM);
for ii=1:a1
  for jj=1:a2
    dv=datevec(TM(ii,jj));
    MM(ii,jj)=dv(2);
  end
end



for im=1:12
  I=find(MM==im);

  mVflx(im)=nanmean(VFlx(I));
  Vflx_p1(im)=prctile(VFlx(I),25);
  Vflx_p2(im)=prctile(VFlx(I),75);

  mFW1(im)=nanmean(FWf1(I));
  FW1_p1(im)=prctile(FWf1(I),25);
  FW1_p2(im)=prctile(FWf1(I),75);

  mFW2(im)=nanmean(FWf2(I));
  FW2_p1(im)=prctile(FWf2(I),25);
  FW2_p2(im)=prctile(FWf2(I),75);

  mFW3(im)=nanmean(FWf3(I));
  FW3_p1(im)=prctile(FWf3(I),25);
  FW3_p2(im)=prctile(FWf3(I),75);

  mHF1(im)=nanmean(Hf1(I));
  HF1_p1(im)=prctile(Hf1(I),25);
  HF1_p2(im)=prctile(Hf1(I),75);

  mHF2(im)=nanmean(Hf2(I));
  HF2_p1(im)=prctile(Hf2(I),25);
  HF2_p2(im)=prctile(Hf2(I),75);

end

POS=[0.08 0.7 0.85 0.25; ...
     0.08 0.4 0.85 0.25;...
     0.08 0.1 0.85 0.25];


imo=[1:12];
figure(1); clf;

% Vol flux:
ps=POS(1,:);
axes('Position',ps);
hold on;
clr=[0 0.4 0.7];
plot(imo,mVflx,'-','Color',clr,'Linewidth',2);
plot(imo,mVflx,'.','Color',clr,'Markersize',14);
pc1=Vflx_p1;
pc2=Vflx_p2;
clrp=[0.6 0.6 0.6];
plot(imo,pc1,'--','Color',clrp,'Linewidth',1);
plot(imo,pc2,'--','Color',clrp,'Linewidth',1);
%for im=1:12
%  plot([im im],[pc1(im) pc2(im)],'Linewidth',1.6,'Color',clrp);
%  plot(im,pc1(im),'.','Color',clrp,'Markersize',14);
%  plot(im,pc2(im),'.','Color',clrp,'Markersize',14);
%end

set(gca,'tickdir','out',...
        'xlim',[0.8 12.2],...
        'ylim',[0 6],...
        'xtick',[1:12],...
        'ytick',[0:5],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',12);
stl=sprintf('HYCOM 0.08-112, VolF, Sv, BSO, %i-%i',YR1,YR2);
title(stl);
Vmn=nanmean(mVflx);
Smn=std(nanmean(VFlx,2)); % std of annual mean
tll=sprintf('VFlux=%4.1f +/-%4.1f Sv',Vmn,Smn);
ttx=text(1,1,tll,'Fontsize',12);

% Heat Flux 1:
ps=POS(2,:);
axes('Position',ps);
hold on;
clr=[0.9 0.4 0.];
plot(imo,mHF1,'-','Color',clr,'Linewidth',2);
plot(imo,mHF1,'.','Color',clr,'Markersize',14);
pc1=HF1_p1;
pc2=HF1_p2;
clrp=[0.4 0.4 0.4];
plot(imo,pc1,'--','Color',clrp,'Linewidth',1);
plot(imo,pc2,'--','Color',clrp,'Linewidth',1);
set(gca,'tickdir','out',...
        'xlim',[0.8 12.2],...
        'ylim',[0 250],...
        'xtick',[1:12],...
        'ytick',[0:50:300],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',12);
stl=sprintf('HYCOM 0.08-112, HeatF (T=%3.1fC), TW, BSO, %i-%i',Tref1,YR1,YR2);
title(stl);
Hmn=nanmean(mHF1);
Sth=std(nanmean(Hf1,1));
tll=sprintf('HFlux=%4.1f +/-%4.1f TW',Hmn,Sth);
ttx=text(1,30,tll,'Fontsize',12);

% Heat Flux 2:
ps=POS(3,:);
axes('Position',ps);
hold on;
clr=[0.9 0.4 0.];
plot(imo,mHF2,'-','Color',clr,'Linewidth',2);
plot(imo,mHF2,'.','Color',clr,'Markersize',14);
pc1=HF2_p1;
pc2=HF2_p2;
clrp=[0.4 0.4 0.4];
plot(imo,pc1,'--','Color',clrp,'Linewidth',1);
plot(imo,pc2,'--','Color',clrp,'Linewidth',1);
set(gca,'tickdir','out',...
        'xlim',[0.8 12.2],...
        'ylim',[0 250],...
        'xtick',[1:12],...
        'ytick',[0:50:300],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',12);
stl=sprintf('HYCOM 0.08-112, HeatF (T=%3.1fC), TW, BSO, %i-%i',Tref2,YR1,YR2);
title(stl);
Hmn=nanmean(mHF2);
Sth=std(nanmean(Hf2,1));
tll=sprintf('VFlux=%4.1f +/-%4.1f TW',Hmn,Sth);
ttx=text(1,30,tll,'Fontsize',12);


bottom_text(btx,'pwd',1);



% ==============================
%  FW Fluxes
% ==============================
figure(2); clf;

% Vol flux:
ps=POS(1,:);
axes('Position',ps);
hold on;
clr=[0 0.4 0.9];
plot(imo,mFW1,'-','Color',clr,'Linewidth',2);
plot(imo,mFW1,'.','Color',clr,'Markersize',14);
pc1=FW1_p1;
pc2=FW1_p2;
clrp=[0.4 0.4 0.4];
plot(imo,pc1,'--','Color',clrp,'Linewidth',1);
plot(imo,pc2,'--','Color',clrp,'Linewidth',1);

set(gca,'tickdir','out',...
        'xlim',[0.8 12.2],...
        'ylim',[-20 30],...
        'xtick',[1:12],...
        'ytick',[-30:10:50],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',12);
stl=sprintf('HYCOM 0.08-112, FWFlx (%5.2f), mSv, BSO, %i-%i',Sref1,YR1,YR2);
title(stl);
Fmn=nanmean(mFW1);
Sfw=std(nanmean(FWf1,2)); % std of annual mean
tll=sprintf('FWF=%5.1f +/-%5.1f mSv',Fmn,Sfw);
ttx=text(1,-15,tll,'Fontsize',12);


ps=POS(2,:);
axes('Position',ps);
hold on;
clr=[0 0.4 0.9];
plot(imo,mFW2,'-','Color',clr,'Linewidth',2);
plot(imo,mFW2,'.','Color',clr,'Markersize',14);
pc1=FW2_p1;
pc2=FW2_p2;
clrp=[0.4 0.4 0.4];
plot(imo,pc1,'--','Color',clrp,'Linewidth',1);
plot(imo,pc2,'--','Color',clrp,'Linewidth',1);

set(gca,'tickdir','out',...
        'xlim',[0.8 12.2],...
        'ylim',[-10 45],...
        'xtick',[1:12],...
        'ytick',[-20:10:50],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',12);
stl=sprintf('HYCOM 0.08-112, FWFlx (%5.2f), mSv, BSO, %i-%i',Sref2,YR1,YR2);
title(stl);
Fmn=nanmean(mFW2);
Sfw=std(nanmean(FWf2,2)); % std of annual mean
tll=sprintf('FWF=%5.1f +/-%5.1f mSv',Fmn,Sfw);
ttx=text(1,30,tll,'Fontsize',12);


ps=POS(3,:);
axes('Position',ps);
hold on;
clr=[0 0.4 0.9];
plot(imo,mFW3,'-','Color',clr,'Linewidth',2);
plot(imo,mFW3,'.','Color',clr,'Markersize',14);
pc1=FW3_p1;
pc2=FW3_p2;
clrp=[0.4 0.4 0.4];
plot(imo,pc1,'--','Color',clrp,'Linewidth',1);
plot(imo,pc2,'--','Color',clrp,'Linewidth',1);

set(gca,'tickdir','out',...
        'xlim',[0.8 12.2],...
        'ylim',[-10 45],...
        'xtick',[1:12],...
        'ytick',[-20:10:50],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',12);
stl=sprintf('HYCOM 0.08-112, FWFlx (S<%5.2f), mSv, BSO, %i-%i',Sref1,YR1,YR2);
title(stl);
Fmn=nanmean(mFW3);
Sfw=std(nanmean(FWf3,2)); % std of annual mean
tll=sprintf('FWF=%5.1f +/-%5.1f mSv',Fmn,Sfw);
ttx=text(1,25,tll,'Fontsize',12);
  
bottom_text(btx,'pwd',1);


%=================================

cl1=flipud(colormap_blue(200));
cl2=colormap_red(200);
for ik=1:5
  cl1(end-(ik-1),:)=[1 1 1];
  cl2(ik,:)=[1 1 1];
end
cmp1=[cl1;cl2];
cmp1=smooth_colormap(cmp1,5);


% To avoid "zig-zigzagging" project onto Y-axis
dst=cumsum(dL);
XX=[0;dst];
dI=diff(IIs);
Ib=find(dI>0);
Ing=find(dI==0);


UUm=sub_fill_bottom_nans(UUm);
SSm=sub_fill_bottom_nans(SSm);
TTm=sub_fill_bottom_nans(TTm);

Hb0=Hb(Ib);
Um=UUm(:,Ib);
Tm=TTm(:,Ib);
Sm=SSm(:,Ib);
%X=XX(Ib);
X=Slat(Ib);

Xbtm=[X(1),X,X(end)];
Zbtm=[-1000,Hb0,-1000];
Z=ZZi;


xl1=70.1;
xl2=76.6;
yl1=-500;
yl2=0;


% Plot 2D sections:
figure(3); clf;
axes('Position',[0.1 0.4 0.8 0.5]);

pcolor(X,Z,Um); shading interp;
colormap(cmp1);
caxis([-0.15 0.15]);
hold on;
fill(Xbtm,Zbtm,[0 0 0]);
hb=colorbar;
set(hb,'Fontsize',12,...
       'Position',[0.92 0.4 0.01 0.5]);

set(gca,'tickdir','out',...
        'xtick',[70:1:78],...
        'xdir','reverse',...
        'ytick',[-500:50:0],...
        'xlim',[xl1 xl2],...
        'ylim',[yl1 yl2],...
        'Fontsize',12);
stl=sprintf('HYCOM 0.08-112, BSO, U m/s %i-%i',YR1,YR2);
title(stl,'Fontsize',12);

bottom_text(btx,'pwd',1,'Position',[0.08 0.3 0.4 0.05]);



nint= 200;
CMP = create_colormap5(nint,c1,c2);
cmp0= CMP.colormap;
for k=1:round(0.06*nint);
  cmp0(k,:)=[1 0.4 1];
end
nav = 15;
cmp = smooth_colormap(cmp0,nav);
%cmp=cmp0;
cnt = CMP.intervals;


s0=35.6;
flg=1;
AA=salin2exp(s0,Sm,flg);
c1=min(min(AA));
c2=max(max(AA));
cntrS=[34.8:0.05:35];  
cntr=exp((35.6-cntrS).^(-1));

figure(4); clf;
axes('Position',[0.1 0.4 0.8 0.5]);

pcolor(X,Z,AA); shading interp;
colormap(cmp);
caxis([c1 c2]);
hold on;
fill(Xbtm,Zbtm,[0 0 0]);
hb=colorbar;
set(hb,'Fontsize',12,...
       'Position',[0.92 0.4 0.01 0.5]);
  tkk=get(hb,'Ticks');
  for ik=1:length(tkk)
    stk(ik)=salin2exp(s0,tkk(ik),-1);
    stk(ik)=(round(stk(ik)*100))/100;
    lbstk{ik}=sprintf('%4.2f',stk(ik));
  end;

  set(hb,'TickLabels',lbstk);


set(gca,'tickdir','out',...
        'xtick',[70:1:78],...
        'xdir','reverse',...
        'ytick',[-500:50:0],...
        'xlim',[xl1 xl2],...
        'ylim',[yl1 yl2],...
        'Fontsize',12);

stl=sprintf('HYCOM 0.08-112, BSO, S %i-%i',YR1,YR2);
title(stl,'Fontsize',12);

bottom_text(btx,'pwd',1,'Position',[0.08 0.3 0.4 0.05]);




figure(5); clf;
axes('Position',[0.1 0.4 0.8 0.5]);

pcolor(X,Z,Tm); shading interp;
colormap(cmp);
caxis([0 8]);
hold on;
fill(Xbtm,Zbtm,[0 0 0]);
hb=colorbar;
set(hb,'Fontsize',12,...
       'Position',[0.92 0.4 0.01 0.5]);

set(gca,'tickdir','out',...
        'xtick',[70:1:78],...
        'xdir','reverse',...
        'ytick',[-500:50:0],...
        'xlim',[xl1 xl2],...
        'ylim',[yl1 yl2],...
        'Fontsize',12);
stl=sprintf('HYCOM 0.08-112, T, U m/s %i-%i',YR1,YR2);
title(stl,'Fontsize',12);

bottom_text(btx,'pwd',1,'Position',[0.08 0.3 0.4 0.05]);




