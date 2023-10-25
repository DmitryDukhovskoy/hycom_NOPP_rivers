% Note: double check the code 
% It seems to be errors in flux calculation
% at the "steps" segments
%
% Correct code - see greenl_fluxes_POPboxes.m
%
% The problem is that when say V>0 and U<0 
% the flow crosses the segment 
%
%         *-------
%         |
%       <---  V<0  negative flux (toW) should be 
%         |        positive in sign convention
%    -----*        as in Fram Str. positive flux means
%                  Vol/Fw transp towards N Pole
%    South
%
%
% Quick fix of the error of wrong
% transport sign at the step segments
% need to change signs of the error segments
% but probably not all of them?  see Jerr index of all 
% wrong oriented transports
%
%
% Calculate volume, heat, fw fluxes
% seasonality etc.
% use 3z interpolated fields
% see interp2z_UTS_Fram.m
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear

expt = 112;
TV   = 11;
sctnm= 'Fram';

YR1=2005;
YR2=2016;

Sref=34.9; 
rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
btx='anls_UTS_Fram.m';

btx='anls_trnsp_Fram.m';

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


% combine fields:
cc=0;
VV=struct;
for yr=YR1:YR2
  fmat_out = sprintf('%sarc08-%3.3i_UTSonZ_daily_%s_%4.4i.mat',...
                 pthmat,expt,sctnm,yr);
  fprintf('Loading %s\n',fmat_out);
  load(fmat_out);
% SCTZ 
% Compute monthly means
% For plotting, find segments that go east (j=const)
% to avoid jumps in staircases of the section
  if ~exist('Jpr','var');
    JJ=SCTZ.J;
    dJ=diff(JJ);
    Jpr=find(dJ==0);
    Jerr=find(abs(dJ)>0);
  end

  U=SCTZ.NVel;
  TM=SCTZ.TM;
  ndays=size(U,1);
  DV=datevec(TM);
% Monthly means:
  for im=1:12
    cc=cc+1;
    I=find(DV(:,2)==im);
    aa=U(I,:,:);
    aa=squeeze(nanmean(aa,1));
    FRAM(cc).V=aa;
    
    aa=SCTZ.Saln(I,:,:);
    aa=squeeze(nanmean(aa,1));
    FRAM(cc).S=aa; 
    
    aa=SCTZ.Temp(I,:,:);
    aa=squeeze(nanmean(aa,1));
    FRAM(cc).T=aa; 
    
    TMe(cc)=datenum(yr,im,15);
  end
end

% Overall U,T,S Mean:
nrc=cc;
aa=FRAM(1).V*0;
ss=aa;
tt=aa;
for ii=1:nrc
  aa=aa+FRAM(ii).V;
  ss=ss+FRAM(ii).S;
  tt=tt+FRAM(ii).T;
end
aa=aa./nrc;
ss=ss./nrc;
tt=tt./nrc;
VM=aa(:,Jpr);
SM=ss(:,Jpr);
TM=tt(:,Jpr);

% Std:
aa=VM*0;
ss=aa;
tt=aa;
for ii=1:nrc
  dmm=(FRAM(ii).V(:,Jpr)-VM).^2;
  aa=aa+dmm;
  dmm=(FRAM(ii).S(:,Jpr)-SM).^2;
  ss=ss+dmm;
  dmm=(FRAM(ii).T(:,Jpr)-TM).^2;
  tt=tt+dmm;
end
aa=aa./nrc;
ss=ss./nrc;
tt=tt./nrc;
sgmV=sqrt(aa);
sgmS=sqrt(ss);
sgmT=sqrt(tt);

ZZ=SCTZ(1).ZZ;
ZM=SCTZ(1).ZM;
dZ=abs(diff(ZZ));

X=SCTZ.LON(Jpr); % for plotting only X segments
Z=ZZ(2:end);

Xfull=SCTZ.LON; 

cl1=flipud(colormap_blue(200));
cl2=colormap_red(200);
for ik=1:5
  cl1(end-(ik-1),:)=[1 1 1];
  cl2(ik,:)=[1 1 1];
end
cmp1=[cl1;cl2];
cmp1=smooth_colormap(cmp1,5);

% Bottom
JJ=SCTZ.J;
II=SCTZ.I;
dJ=diff(JJ);
cc=0;
Hb=[];
for jj=1:length(dJ)
  i0=II(jj);
  j0=JJ(jj);
  if dJ(jj)==0
    cc=cc+1;
    Hb(cc,1)=HH(j0,i0);
  end
end


VM=sub_fill_bottom_nans(VM);
SM=sub_fill_bottom_nans(SM);
TM=sub_fill_bottom_nans(TM);

Xbtm=[X(1),X,X(end)];
Zbtm=[-5000,Hb',-5000];


SBE = Fram_moorAWI;
nf=length(SBE);


%keyboard

% ================
% Plot mean U and std dv
% ================
figure(1); clf;
pcolor(X,Z,VM); shading flat;
colormap(cmp1);
caxis([-0.2 0.2]);
hold on;
contour(X,Z,sgmV,[0:0.01:0.5],'Color',[0 0 0]);
contour(X,Z,sgmV,[0.01 0.01],'Color',[0 0 0],'Linewidth',1.6);
%plot(X,Hb,'Color',[0 0 0],'Linewidth',4);
fill(Xbtm,Zbtm,[0 0 0]);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  XY=[x,y];
  IJ = sub_XY2indx(XY,LON,LAT);
  zb=HH(IJ(2),IJ(1));
%  plot(IJ(1),IJ(2),'r.','Markersize',16);
  plot(x,zb-100,'r^','Markersize',9,'MarkerFaceColor',[1 0 0]);
  nm=SBE(ik).Name;
  text(x,zb-200,nm,'Fontsize',12,'Color',[1 0 0]);
end


hb=colorbar;
set(hb,'Fontsize',14,...
       'Position',[0.92 0.1 0.02 0.8]);

set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[-3500:500:0],...
	'xlim',[-19.2 11.9],...
	'ylim',[-3150 0],...
	'Fontsize',14);

title('ARCc0.08-112 mean V, StDev (0.01), 2005-2016');
bottom_text(btx,'pwd',1);


% ================
% Plot mean S
% ================
figure(2); clf;
cs1=32;
cs2=35.;
%CMP = colormap_sclr2(200,cs1,cs2);
%cmpS1=CMP.colormap;
cmpS=colormap(parula(400));

pcolor(X,Z,SM); shading flat;
colormap(cmpS);
caxis([cs1 cs2]);
hold on;
%contour(X,Z,sgmS,[0:0.1:1],'Color',[0 0 0]);
%contour(X,Z,sgmS,[0.01:0.02:0.099],'--','Color',[0 0 0]);
%contour(X,Z,sgmS,[0.01 0.01],'--','Color',[0 0 0],'Linewidth',1.6);
%contour(X,Z,SM,[34.9:0.01:34.99],'k-','Color',[0.4 0.4 0.4]);
%contour(X,Z,SM,[34.95 34.95],'k-','Color',[0.4 0.4 0.4],'Linewidth',1.6);
contour(X,Z,SM,[34.8:0.02:36.0],'k-','Color',[0.3 0.3 0.3]);
contour(X,Z,SM,[34.94 34.94],'k-','Color',[0.3 0.3 0.3],'Linewidth',1.6);
contour(X,Z,SM,[34.5 34.5],'k-','Color',[0.6 0.6 0.6],'Linewidth',1);
contour(X,Z,SM,[32.5:0.5:34.],'k-','Color',[0.6 0.6 0.6],'Linewidth',1);

%plot(X,Hb,'Color',[0 0 0],'Linewidth',4);
fill(Xbtm,Zbtm,[0 0 0]);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  XY=[x,y];
  IJ = sub_XY2indx(XY,LON,LAT);
  zb=HH(IJ(2),IJ(1));
%  plot(IJ(1),IJ(2),'r.','Markersize',16);
  plot(x,zb-100,'r^','Markersize',9,'MarkerFaceColor',[1 0 0]);
  nm=SBE(ik).Name;
  text(x,zb-200,nm,'Fontsize',12,'Color',[1 0 0]);
end


hb1=colorbar;
set(hb1,'Fontsize',14,...
       'Position',[0.92 0.1 0.02 0.8]);

set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[-3500:500:0],...
	'xlim',[-19.2 11.9],...
	'ylim',[-3150 0],...
	'Fontsize',14);

title('ARCc0.08-112 mean S, StDev, 2005-2016');
bottom_text(btx,'pwd',1);


% ==========================  
%   Plot Mean T 
% ==========================  
figure(5); clf;
ct1=-2;
ct2=5;
nint=400;
CMP = create_colormap5(nint,ct1,ct2);
cmpT = CMP.colormap;

pcolor(X,Z,TM); shading flat;
colormap(cmpT);
caxis([ct1 ct2]);
hold on;
%contour(X,Z,sgmT,[0:0.5:2],'Color',[0 0 0]);
%contour(X,Z,sgmT,[0.0:0.05:0.49],'--','Color',[0 0 0]);
%contour(X,Z,sgmT,[0.05 0.05],'--','Color',[0 0 0],'Linewidth',1.6);
contour(X,Z,TM,[-2:1:8],'k-','Color',[0.6 0.6 0.6]);
contour(X,Z,TM,[0 0],'k-','Color',[0.6 0.6 0.6],'Linewidth',1.6);

%plot(X,Hb,'Color',[0 0 0],'Linewidth',4);
fill(Xbtm,Zbtm,[0 0 0]);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  XY=[x,y];
  IJ = sub_XY2indx(XY,LON,LAT);
  zb=HH(IJ(2),IJ(1));
%  plot(IJ(1),IJ(2),'r.','Markersize',16);
  plot(x,zb-100,'r^','Markersize',9,'MarkerFaceColor',[1 0 0]);
  nm=SBE(ik).Name;
  text(x,zb-200,nm,'Fontsize',12,'Color',[1 0 0]);
end


hb1=colorbar;
set(hb1,'Fontsize',14,...
       'Position',[0.92 0.1 0.02 0.8]);

set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[-3500:500:0],...
	'xlim',[-19.2 11.9],...
	'ylim',[-3150 0],...
	'Fontsize',14);

title('ARCc0.08-112 mean T, StDev, 2005-2016');
bottom_text(btx,'pwd',1);


% ==========================  
%
%  Analysis of FW/Vol fluxes
%
% ==========================  

% Monthly means:
%  for im=1:12
%    I=find(DV(:,2)==im);
%    aa=U(I,:,:);
%    aa=squeeze(nanmean(aa,1));
%    cc=cc+1;
%    VV(cc).V=aa(:,Jpr);
%    TMe(cc)=datenum(yr,im,15);
%  end

% Volume flux:
ZZ=SCTZ(1).ZZ;
ZM=SCTZ(1).ZM;
dZ=abs(diff(ZZ));
Is=SCTZ(1).I;
Js=SCTZ(1).J;
clear Acell;
% Area of the grid cells:
[ntm,nlyrs,nsgm]=size(U);
for isgm=1:nsgm
  i1=Is(isgm);
  i2=Is(isgm+1);
  j1=Js(isgm);
  j2=Js(isgm+1);
  if j1==j2
    dL=DX(j1,i1);
  else
    dL=DY(j1,i1);
  end
  Acell(:,isgm)=dL*dZ;
end


% Calculate climatology
% and total transport:
% and FW flux, Sref
cc=0;
cs=0;
cw=0;
clear FWmo
FWmo=struct;
Vmo=struct;
for im=1:12
  FWmo(im).ccm=0;
  Vmo(im).ccm=0;
end

nrc=length(FRAM);
clear avVb Tr cTr FW FWS FWW cumFW
for im=1:nrc
  cc=cc+1;
% Quick fixing an error of wrong
% transport sign at the step segments
% need to change signs of the error segments
% but probably not all of them?
  V=FRAM(im).V;
  V(:,Jerr)=-V(:,Jerr);
  aa=nanmean(V);
%  saa=sign(aa);
%  naa=length(saa);
%  cer=0;
%  clear Jerr2
%  for ik=2:naa-1
%    if saa(ik)~=saa(ik-1) & saa(ik)~=saa(ik+1)
%      cer=cer+1;
%      Jerr2(cer)=ik;
%    end
%  end
  
  
  S=FRAM(im).S;
  Fvol=V.*Acell;
  a=nansum(Fvol);
%
%  a(Jerr)=-a(Jerr);
  avVb(cc,:)=a(Jpr);  % depth-integr flow
%  Tr(cc,1)=nansum(a(Jpr));
%  cTr(cc,:)=cumsum(a(Jpr));
  Vfull=a;  % depth-integr flow
  Tr(cc,1)=nansum(a);
  cTr(cc,:)=cumsum(a);
  fwf=Fvol.*(Sref-S)/Sref;
  fwf(S>Sref)=nan;
  fwf(S<=1e-10)=nan;
%  fwf(Jerr)=-fwf(Jerr);
  fwf=nansum(fwf);
%  cumFW(cc,:)=cumsum(fwf(Jpr));
  FW(cc,:)=fwf(Jpr); % depth-integrated FWFlux
  cumFW(cc,:)=cumsum(fwf);
%  FW(cc,:)=fwf; % depth-integrated FWFlux
  
  dv=datevec(TMe(im));
  jmo=dv(2);
  ccm=FWmo(jmo).ccm;
  ccm=ccm+1;
%  FWmo(jmo).FW(ccm,:)=fwf(Jpr); % fwf flux by months
  FWmo(jmo).FW(ccm,:)=fwf; % fwf flux by months
  FWmo(jmo).ccm=ccm;
  Vmo(jmo).Vol(ccm,:)=Vfull; % Vol flux by months
% Season:
  if dv(2)>=6 & dv(2)<=10
    cs=cs+1;
    FWS(cs,:)=fwf;
  else
    cw=cw+1;
    FWW(cw,:)=fwf;
  end
end

cl2=colormap_orange(200);
cl1=flipud(colormap_blue(200));
for ik=1:8
  cl2(ik,:)=[1 1 1];
  cl1(end-(ik-1),:)=[1 1 1];
end
cmpa=[cl1;cl2];
cmpa=smooth_colormap(cmpa,15);

cl1=flipud(colormap_cold(200));
cl2=colormap_green(200);
for ik=1:10
  cl2(ik,:)=[1 1 1];
  cl1(end-(ik-1),:)=[1 1 1];
end
cmpF=[cl1;cl2];
cmpF=smooth_colormap(cmpF,15);

% =================================
% Hovm diagr of depth-integr flow
% ===============================
avVb(avVb==0)=nan;
DVe=datevec(TMe);
Yrm=[DVe(1,1):1/12:DV(end,1)+0.999];

figure(3); clf;
axes('Position',[0.2 0.1 0.4 0.6]);
hold on;
pcolor(X,Yrm,avVb*1e-6); shading flat;
colormap(cmpa);
caxis([-0.6 0.6]);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  plot([x x],[Yrm(1) Yrm(end)],'k--','Color',[0.5 0.5 0.5]);
end

set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[2000:2018],...
	'xlim',[-19.2 11.9],...
	'ylim',[min(Yrm) max(Yrm)+0.1],...
	'Fontsize',14);
		 
hb=colorbar;
set(hb,'Position',[0.63 0.1 0.02 0.6],...
       'Fontsize',14);


% Plot bottom profile
axes('Position',[0.2 0.71 0.4 0.1]);
hold on;
fill(Xbtm,Zbtm,[0 0 0]);
set(gca,'tickdir','out',...
	'xtick',[],...
	'ytick',[-3000:1000:0],...
	'xlim',[-19.2 11.9],...
	'ylim',[-3200 0],...
	'Fontsize',14);
for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  XY=[x,y];
  IJ = sub_XY2indx(XY,LON,LAT);
  zb=HH(IJ(2),IJ(1));
%  plot(IJ(1),IJ(2),'r.','Markersize',16);
  plot(x,zb-100,'r^','Markersize',6,'MarkerFaceColor',[1 0 0]);
  nm=SBE(ik).Name;
  text(x,zb-200,nm,'Fontsize',8,'Color',[1 0 0]);
end

title('Depth-intgr U, Sv, HYCOM 0.08-112');
bottom_text(btx,'pwd',1);
% 


% ------------------------
%
%   Hovm diagr of FW flux
%
% -------------------------
figure(23); clf;
axes('Position',[0.2 0.1 0.4 0.6]);
hold on;
pcolor(X,Yrm,FW*1e-3); shading flat;
colormap(cmpF);
caxis([-4 4]);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  plot([x x],[Yrm(1) Yrm(end)],'k--','Color',[0.5 0.5 0.5]);
end

set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[2000:2018],...
	'xlim',[-19.2 11.9],...
	'ylim',[min(Yrm) max(Yrm)+0.1],...
	'Fontsize',14);
		 
hb=colorbar;
set(hb,'Position',[0.63 0.1 0.02 0.6],...
       'Fontsize',14);

% Plot bottom profile
axes('Position',[0.2 0.71 0.4 0.1]);
hold on;
fill(Xbtm,Zbtm,[0 0 0]);
set(gca,'tickdir','out',...
	'xtick',[],...
	'ytick',[-3000:1000:0],...
	'xlim',[-19.2 11.9],...
	'ylim',[-3200 0],...
	'Fontsize',14);
for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  XY=[x,y];
  IJ = sub_XY2indx(XY,LON,LAT);
  zb=HH(IJ(2),IJ(1));
%  plot(IJ(1),IJ(2),'r.','Markersize',16);
  plot(x,zb-100,'r^','Markersize',6,'MarkerFaceColor',[1 0 0]);
  nm=SBE(ik).Name;
  text(x,zb-200,nm,'Fontsize',8,'Color',[1 0 0]);
end

title('Depth-intgr FWFlux mSv, HYCOM 0.08-112');
bottom_text(btx,'pwd',1);


% ======================
% Mean Tr 
% ======================
cTr=cTr*1e-6; %m3/s -> Sv

McTr=mean(cTr,1);
p10=prctile(cTr,10,1);
p90=prctile(cTr,90,1);

% Fraction of transport not measured:
for ik=1:nf
  Ymr(ik)=SBE(ik).lat_lon(1);
  Xmr(ik)=SBE(ik).lat_lon(2);
end
x=min(Xmr);

%ixW=max(find(X<=x)); % westernmost location
%TrW=cTr(:,ixW-1); % transp in the western shelf section not covered by obs
ixW=max(find(Xfull<=x)); % westernmost location
TrW=cTr(:,ixW-1); % transp in the western shelf section not covered by obs

x=min(Xmr(Xmr>-8)); % excluding F17 (shelf mooring)
ixW2=max(find(Xfull<=x));
TrW2=cTr(:,ixW2-1); % transp in the western shelf section not covered by obs

x=max(Xmr);
ixE=min(find(Xfull>=x)); % easternmost location in x-sgements arrays only
TrE=cTr(:,end)-cTr(:,ixE);
% 
% Trt within observed moorings
TrObs=cTr(:,ixE)-cTr(:,ixW-1);

% Indices in full x-y segments across Fram
x=min(Xmr);
ixyW=max(find(Xfull<=x));
x=min(Xmr(Xmr>-8)); % excluding F17 (shelf mooring)
ixyW2=max(find(Xfull<=x));
x=max(Xmr);
ixyE=min(find(Xfull>=x)); % easternmost location in x-sgements arrays only


figure(4); clf;
axes('Position',[0.1 0.55 0.82 0.35]);
hold on;
plot(X,McTr(Jpr),'Linewidth',2.2);
plot(X,p10(Jpr),'Color',[0 .8 1],'Linewidth',1.6);
plot(X,p90(Jpr),'Color',[0 .8 1],'Linewidth',1.6);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
%  plot([x x],[Yrm(1) Yrm(end)],'k--','Color',[0.7 0.7 0.7]);
  plot(x,min(min(cTr)),'r^','Markersize',6,'MarkerFaceColor',[1 0 0]);
end


set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[-20:2:2],...
	'xlim',[-19.2 11.9],...
	'xgrid','on',...
	'ygrid','on',...
	'ylim',[min(min(cTr)) max(max(cTr))],...
	'Fontsize',14);
title('Cumulative Vol Transport, Sv, Fram');


%
axes('Position',[0.1 0.1 0.1 0.3]);
hx=boxplot(TrW);
set(hx,'LineWidth',1);
set(gca,'tickdir','out',...
	'xtick',[],...
	'ytick',[-1.8:0.2:1],...
	'ylim',[-1.4 0.2],...
	'Fontsize',14);
title('Trnsp Shore-F17');

axes('Position',[0.32 0.1 0.1 0.3]);
hx=boxplot(TrW2);
set(hx,'LineWidth',1);
set(gca,'tickdir','out',...
	'xtick',[],...
	'ytick',[-1.8:0.2:1],...
	'ylim',[-1.4 0.2],...
	'Fontsize',14);
title('Trnsp Shore-F14');

axes('Position',[0.55 0.1 0.1 0.3]);
hx=boxplot(TrObs);
set(hx,'LineWidth',1);
set(gca,'tickdir','out',...
	'xtick',[],...
	'ytick',[-10:2:3],...
	'ylim',[-9 3],...
	'Fontsize',14);
title('Trnsp F17-F1');

axes('Position',[0.8 0.1 0.1 0.3]);
hx=boxplot(TrE);
set(hx,'LineWidth',1);
set(gca,'tickdir','out',...
	'xtick',[],...
	'ytick',[-0.5:0.05:0.5],...
	'ylim',[-0.2 0.3],...
	'Fontsize',14);
title('Trnsp F1-shore');

bottom_text(btx,'pwd',1);

% =========================================
%  Fresh water flux
% =========================================
%
cumFW=cumFW*1e-3; % mSv
FWW=FWW*1e-3;
FWS=FWS*1e-3;

Mfwf=mean(cumFW,1);
p10=prctile(cumFW,10,1);
p90=prctile(cumFW,90,1);
Wfwf=cumsum(mean(FWW,1));
Sfwf=cumsum(mean(FWS,1));

% Fraction of transport not measured:
fwW=FWW(:,ixW-1); % winter transp in the western shelf section not covered by obs
fwS=FWS(:,ixW-1); % summer transp in the western shelf section not covered by obs

fwW2=FWW(:,ixW2-1); % F17-F14
fwS2=FWS(:,ixW2-1); % 

feW=FWW(:,end)-FWW(:,ixE);  % unobserved on Sptitsb shelf winter
feS=FWS(:,end)-FWS(:,ixE);  % Spitsb shelf summer
% 
% FW Trt within observed moorings
FWObs=FWW(:,ixE)-FWW(:,ixW-1);
FSObs=FWS(:,ixE)-FWS(:,ixW-1);


% cumFWlux
figure(6); clf;
axes('Position',[0.1 0.55 0.82 0.35]);
hold on;
plot(X,Mfwf(Jpr),'Color',[0. 0. 0.],'Linewidth',2.2);
plot(X,Wfwf(Jpr),'Color',[0. 0.5 1],'Linewidth',2.);
plot(X,Sfwf(Jpr),'Color',[1.0 0.6 0.],'Linewidth',2.);
plot(X,p10(Jpr),'Color',[0.6 .6 .6],'Linewidth',1.6);
plot(X,p90(Jpr),'Color',[0.6 .6 .6],'Linewidth',1.6);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
%  plot([x x],[Yrm(1) Yrm(end)],'k--','Color',[0.7 0.7 0.7]);
  plot(x,min(min(cumFW)),'r^','Markersize',6,'MarkerFaceColor',[1 0 0]);
end

set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[-200:40:100],...
	'xlim',[-19.2 11.9],...
	'xgrid','on',...
	'ygrid','on',...
	'ylim',[min(min(cumFW)) max(max(cumFW))],...
	'Fontsize',14);
title('Mean, Summer, winter Cumulative FW Flux, mSv, Fram');

  
axes('Position',[0.1 0.1 0.1 0.3]);
ll=length(fwW);
dmm=fwS;
dmm(end:ll)=nan;
A=[fwW,dmm];
hx=boxplot(A);
set(hx,'LineWidth',1);
set(gca,'tickdir','out',...
	'xtick',[1 2],...
	'ytick',[-80:10:20],...
	'ylim',[-75 10],...
	'xticklabel',{'W';'S'},...
	'Fontsize',14);
title('cumFWlux Shore-F17');

axes('Position',[0.32 0.1 0.1 0.3]);
dmm=fwS2;
dmm(end:ll)=nan;
A=[fwW2,dmm];
hx=boxplot(A);
set(hx,'LineWidth',1);
set(gca,'tickdir','out',...
	'xtick',[1 2],...
	'ytick',[-100:10:0],...
	'ylim',[-95 0],...
	'xticklabel',{'W';'S'},...
	'Fontsize',14);
title('cumFWlux Shore-F14');

axes('Position',[0.55 0.1 0.1 0.3]);
dmm=FSObs;
dmm(end:ll)=nan;
A=[FWObs,dmm];
hx=boxplot(A);
set(hx,'LineWidth',1);
set(gca,'tickdir','out',...
	'xtick',[1 2],...
	'ytick',[-120:10:-20],...
	'ylim',[-118 -55],...
	'xticklabel',{'W';'S'},...
	'Fontsize',14);
title('cumFWlux F17-F1');

axes('Position',[0.8 0.1 0.1 0.3]);
dmm=feS;
dmm(end:ll)=nan;
A=[feW,dmm];
hx=boxplot(A);
set(hx,'LineWidth',1);
set(gca,'tickdir','out',...
	'xtick',[1 2],...
	'ytick',[-0.9:0.2:0.9],...
	'ylim',[-0.55 0.65],...
	'xticklabel',{'W';'S'},...
	'Fontsize',14);
title('cumFWlux F1-shore');

bottom_text(btx,'pwd',1);



% Plot monthly FW transport
% and contributions for different segments
for imo=1:12
  aa=FWmo(imo).FW*1e-3; % mSv
  a1=nansum(aa(:,1:ixyW-1),2);  % West segm on Greenland shelf
  mF1(imo)=mean(a1);
  a2=nansum(aa(:,ixyW:ixE-1),2); % covered with moorings
  mF2(imo)=mean(a2);
  a3=nansum(aa(:,ixyE:end),2);
  mF3(imo)=mean(a3);
  
  ftot=nansum(aa,2);
  mF(imo)=mean(ftot);
  p1(imo)=prctile(ftot,25);
  p2(imo)=prctile(ftot,75);
  
end

FF=[mF1', mF2', mF3'];
pp=mF2./mF; % fraction of observed FW flux

cmp1=[0. 0. 0.4; ...
	  0. 0.4 0.8; ...
	  0  0.8 1];

% Volume transport:
for imo=1:12
  aa=Vmo(imo).Vol*1e-6; 
  a1=nansum(aa(:,1:ixyW-1),2);  % West segm on Greenland shelf
  mV1(imo)=mean(a1);
  a2=nansum(aa(:,ixyW:ixyE-1),2); % covered with moorings
  mV2(imo)=mean(a2);
  a3=nansum(aa(:,ixyE:end),2);
  mV3(imo)=mean(a3);
  
  ftot=nansum(aa,2);
  mV(imo)=mean(ftot);
end  
FVV=[mV1',mV2',mV3'];
ppv=abs(mV2)./(abs(mV1)+abs(mV2)+abs(mV3)); % fraction of observed volume

cmp2=[0 0.3 0; ...
      0 0.7 0.3;...
      0 1 0.6];


figure(10); clf;
axes('Position',[0.1 0.55 0.8 0.35]);
hold on;
HB=bar(FF,'stacked');
colormap(cmp1);
for ib=1:3
  HB(ib).EdgeColor='none';
  HB(ib).FaceColor=cmp1(ib,:);
  HB(ib).BarWidth=0.95;
end
for imo=1:12
  stx=sprintf('%3.2f',pp(imo));
  text(imo-0.3,0.7*mF(imo),stx,'Fontsize',14);
end

set(gca,'tickdir','out',...
	'Fontsize',14,...
	'xlim',[0.5 12.5],...
	'ylim',[-150 0],...
	'xtick',[1:12],...
	'ytick',[-150:25:0]);


hcb=colorbar;
set(hcb,'Position',[0.91 0.55 0.014 0.35],...
	'Ticks',[1/6 1/2 5/6],...
	'TickLabels',{'S1','S2','S3'},...
	'Fontsize',14);
stl=sprintf('0.08-112 HYCOM, FW Flux %i-%i',YR1,YR2);
title(stl);

bottom_text(btx,'pwd',1,'position',[0.08 0.4 0.4 0.05]);


% Volume flux - Note there is 
% an error in the flux calculation
% need to recompute
% transport should be ~-2 Sv
figure(11); clf;
axes('Position',[0.1 0.55 0.8 0.35]);
hold on;
HB=bar(FVV,'stacked');
colormap(cmp2);

for ib=1:3
  HB(ib).EdgeColor='none';
  HB(ib).FaceColor=cmp2(ib,:);
  HB(ib).BarWidth=0.95;
end
for imo=1:12
  stx=sprintf('%3.2f',ppv(imo));
  text(imo-0.3,0.7*mV(imo),stx,'Fontsize',14);
end

set(gca,'tickdir','out',...
	'Fontsize',14,...
	'xlim',[0.5 12.5],...
	'ylim',[-4 0],...
	'xtick',[1:12],...
	'ytick',[-4:0.5:0]);


hcb=colorbar;
set(hcb,'Position',[0.91 0.55 0.014 0.35],...
	'Ticks',[1/6 1/2 5/6],...
	'TickLabels',{'S1','S2','S3'},...
	'Fontsize',14);
stl=sprintf('0.08-112 HYCOM, VolTrt %i-%i',YR1,YR2);
title(stl);

bottom_text(btx,'pwd',1,'position',[0.08 0.4 0.4 0.05]);


%
% Find monthly regression
% of obs at F17 and Gr shelf
cff=1; % scale coefficient for better fit
for imo=1:12
  aa=FWmo(imo).FW*1e-3; % mSv
  a1=nansum(aa(:,1:ixyW-1),2)/cff;  % West segm on Greenland shelf
  a17=aa(:,ixyW);
%  a2=nansum(aa(:,ixyW:ixE-1),2)/cff; % covered with moorings
  
  XX=[ones(length(a17),1),a17];
  [bb,bint,r,rint,stats]=regress(a1,XX);
  
  CF(imo,1)=bb(1);
  CF(imo,2)=bb(2);
  CF(imo,3)=stats(1); % R2
  CF(imo,4)=stats(3); % p-value
end
  





