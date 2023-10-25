% Calculate FW flux
% fluxes through straits
% FW tracers
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

rg  = 9806;
hgg = 1e20; %
Sref= 34.8;

expt = '110';
TV   = 11;  % topo version
%segm = 'BeringS';
segm = 'FramS';
%segm = 'BarOp';
%segm = 'DavisS';

YR1  = 1993;
YR2  = 2016;

fprintf('Analyzing %s\n',segm);

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);

btx = 'calc_FWFlux_straits.m';


fmat = sprintf('%s%s_UV_straits_%4.4i_%4.4i.mat',...
	       pthmat,expt,YR1,YR2);
fprintf('Loading %s\n',fmat);
load(fmat);  % UV
nuv=length(UV);

for ik=1:nuv
  nm=UV(ik).Name;
  if strncmp(nm,segm,4),
    break;
  end
end

% Choose positive flux, + to AO:
switch(segm)
 case('BeringS')
  fpos = -1;
 case('FramS')
  fpos = 1;
 case('BarOp');
  fpos = 1;
 case('DavisS');
  fpos = 1;
end


% Load T/S for straits
fmatS = sprintf('%s%s_TS_straits_%4.4i_%4.4i.mat',...
	       pthmat,expt,YR1,YR2);
fprintf('Loading %s\n',fmatS);
load(fmatS);
TMs=SGMTS(1).TM;

% Calculate normal flux (here, + is upward, - downward)
% need to swap +/- for Bering Str. to
% make positive flux northward (downward on the grid)
TM=UV(ik).TM;
uv=UV(ik).UV_normal;
ZZ=UV(ik).ZZ;
dx=UV(ik).Dist;
nrc=length(TM);
DV=datevec(TM);
xyr = (TM-TM(1))/365.25+DV(1,1);
[a1,ll,nn] = size(uv);

ss=SGMTS(ik).Saln;
[a2,lS,nS] = size(ss);

if lS~=ll | nS~=nn
  error('Check UV & S Segments - dimensions do not agree ...');
end

% Check that time seris match:
if TM(1)~=TMs(1) | TM(end)~=TMs(end)
  fprintf('S and UV time series do not match ...\n');
  error('  Check TIME - time series ...');
end


nlv=ll;
xx = cumsum(dx);
[DX,dmm]=meshgrid(dx,(1:nlv));
[XX,dmm]=meshgrid(xx,(1:nlv));

% Find months:
clear IMM
for im=1:12
%  IMM(im).im=find(DV(:,2)==im);
  IMM(im).cnc=0;
  IMM(im).V=0;
  IMM(im).S=0;
end

% Deriving mean V (m/s) and S
% by months
FWF = [];
for it=1:nrc
  zz=squeeze(ZZ(it,:,:));
  zm=zz(2:end,:);
  zm(zm>0)=0;
  dz=abs(diff(zz,1,1));
  agrd=abs(dz.*DX);
  v=squeeze(uv(it,:,:));
  s=squeeze(ss(it,:,:));
  
  v(v==0)=nan;
  s(s==0)=nan;
  v(dz<1e-3)=nan;
  s(dz<1e-3)=nan;
  
  I=find(~isnan(s));
  Smn = nansum(s(I).*dz(I).*DX(I))/sum(dz(I).*DX(I));
  I=find(~isnan(v));
  Vmn = nansum(v(I).*dz(I).*DX(I))/sum(dz(I).*DX(I));

  im= DV(it,2);
  IMM(im).cnc = IMM(im).cnc+1;
  IMM(im).V   = IMM(im).V+Vmn; % Transport m3/s
  IMM(im).S   = IMM(im).S+Smn;
  
  dS = (Sref-s);
  Trp= fpos*v.*dz.*DX; % fpos - pos/neg normal, makes flux >0 into basin
  FWF(it,1) = nansum(nansum(dS/Sref.*Trp)); % m3/s
  IMM(im).Trp = nansum(nansum(Trp)); % total transport
  
  f_plt=0;
  if f_plt==1
    figure(1); clf;
    axes('Position',[0.08 0.53 0.83 0.4]);
    pcolor(xx,zm,v); shading flat;
    caxis([-0.8 0.8]);
    colorbar

    axes('Position',[0.08 0.08 0.83 0.4]);
    pcolor(xx,zm,s); shading flat;
    caxis([30.5 32.5])
    colorbar
  end
  
end

% Monthly climatology
clear Trp
for im=1:12
  IMM(im).V = IMM(im).V/IMM(im).cnc;
  IMM(im).S = IMM(im).S/IMM(im).cnc;
  II = find(DV(:,2)==im);
  FWFm(im,1) = sum(FWF(II))/length(II);
  Trp(im,1) = IMM(im).Trp;
end  

% Find monthly means:
cc=0;
for iyr=YR1:YR2
  for imo=1:12
    cc=cc+1;
    I=find(DV(:,1)==iyr & DV(:,2)==imo);
    FWflx_mo(cc,1) = mean(FWF(I)); % m3/s
  end
end
% Long-term mean FW transport (km3/yr):
Fmn_yr = mean(FWflx_mo)*3600*24*365*1e-9; % km3/yr
Wn = 1/6; % annual mean
[Bf,Af] = butter(9,Wn,'low');
fFWflx_mo = filtfilt(Bf,Af,FWflx_mo); % 

tmx = [1993:1/12:2016.99];

mnn = min(FWflx_mo);
mxx = max(FWflx_mo);
yl1 = 1.02*min([0, mnn]);
yl2 = 1.02*max([0, mxx]);

figure(1); clf;
axes('Position',[0.08 0.45 0.88 0.43]);
plot(tmx,FWflx_mo);
hold on;
plot(tmx,fFWflx_mo,'r');
set(gca,'tickdir','out',...
	'ylim',[yl1 yl2],...
	'xlim',[1993 2017],...
	'xtick',[1993:2017],...
	'xminortick','on',...
	'xgrid','on',...
	'ygrid','on');
stl = sprintf('ARCc0.08-110, FWFlux (S=%3.1f), m3/sec %s, overall mean: %8.1f km3/yr',...
	      Sref,segm,Fmn_yr);
title(stl);

bottom_text(btx,'pwd',1,'position',[0.02 0.34 0.9 0.05]);


% Monthly climatology FWF flux
% mSv
figure(3); clf;
axes('Position',[0.08 0.55 0.83 0.38]);
plot(FWFm*1e-6); % mSv, + to the Arctic
stl = sprintf('%s, FWFlux, mSv, + to AO, Sref=%3.1f',segm,Sref);
title(stl);
set(gca,'xlim',[1 12]);

dmo = [31,28,31,30,31,30,31,31,30,31,30,31]';
Fkm = FWFm*3600*24.*dmo*1e-9;
Fyr = sum(Fkm);
axes('Position',[0.08 0.08 0.83 0.38]);
plot(Fkm); % km3/mo, + to the Arctic
stl = sprintf('FWFlux, km3/mo, + to AO, km3/yr=%5.1f',Fyr);
title(stl);
set(gca,'xlim',[1 12]);

bottom_text(btx,'pwd',1);

% Plot monthly fluxes:






