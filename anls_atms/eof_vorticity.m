% Calculate EOF of vorticity fields 
% prepared for area-mean vorticity wind fields in 
% calc_mean_vort
%
% from NCEP CFSR and CFSv2
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig = 0;

%iof = 1;
iof = 2;

RR   = 50; % ring size (diamter), km  
%regn = 'natl'; % region
%regn = 'arctic'; % region
regn = 'arctB'; % region

pthdat1 = '/Net/kronos/ddmitry/ncep_cfsr/';
pthdat2 = '/Net/kronos/ddmitry/ncep_cfsv2/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/data_mat/';
pth72   = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/topo_grid/';
pthout  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/fig_atm/';
%fmat    = sprintf('%smonthly_areamean_vort_%s.mat',pthout,regn);
fmat    = sprintf('%smonthly_areamean_vort_%s%3.3ikm.mat',pthout,regn,RR);
btx = 'eof_vorticity.m';

fprintf('Loading vorticity %s\n',fmat);
load(fmat);

% Get indices for interpolating onto HYCOM ARCc0.72 grid:
fgrd = sprintf('%scfsr_gridindx_arc072_nghb.mat',pthmat);
load(fgrd);
mc = INDX.dim_rows;
nc = INDX.dim_colms;

% Topo for ARCc0.72
fsv=[pth72,'new_bath072.mat'];
load(fsv);
LN=elon;
LT=alat;
HH=hnew;
HH(isnan(HH))=10;
clear hnew elon alat;
[mm,nn]=size(HH);
[DX,DY]=sub_dx_dy(LN,LT);
[II,JJ] = meshgrid([1:nn],[1:mm]);

VR   = VRTC.VRT;
Iocn = VRTC.Index_Ocean;
TM   = VRTC.TM;
dT   = TM(2)-TM(1); % output freq., days
rr0  = VRTC.Ring_Size;
if rr0~=RR,
  error('%s ring radius does not match title',fld);
end

% detrend time series:
% not QSCAT due to nans
[aa,bb]=size(VR);
%for kk=1:bb
%  cmm=VR(:,kk);
%  mn=mean(cmm);
%  VR(:,kk)=VR(:,kk)-mn;
%end
fprintf('Detrending\n');
VR=detrend(VR,1);

J0=find(isnan(VR));
VR(J0)=0;

f_plt_vort=0;
if f_plt_vort>0
  inn = 100;
  dnmb=TM(inn);
  DV=datevec(dnmb);
  stl=sprintf('%s vorticity x 1e5, %i/%i/%i',regn,DV(3),DV(2),DV(1));
  Vort=HH*nan;
  Iocn=VRTC.Index_Ocean;
  Vort(Iocn)=VR(inn,:)*1e5;

  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on
  pcolor(Vort); shading flat;
  title(stl);
  bottom_text(btx,'pwd',1);

end

A=VR;
fprintf('Calculating SVD ...\n');
[Up,Sp,Vp] = svd(A);

% Squared diagonal values of Sp (e/values of SVD) 
% are the e/values of covariance matrix F'*F
% e/vectors are in Vp
% also Lmb(i) = pc'*pc;
% where pc - is i-th PC of the e/vector, i.e. amplitude time series 
Lmb = diag(Sp).^2;

iL1=find(Lmb==max(Lmb));
if iL1~=1
  fprintf('Flipping e/values & e/vectors\n');
  Lmb=flipud(Lmb);
  Vp=fliplr(Vp);
end

% e/vectors:
V = Vp(:,iof);

Eof=HH*nan;
Iocn=VRTC.Index_Ocean;
Eof(Iocn)=V;

% Make EOF in Irminger Sea >0 
% for ease of visual comparison
dmm = Eof(40,80);
if dmm<0
  Vp=-Vp;
  V=Vp(:,iof);
  Eof(Iocn)=V;
end

% Pr. Components (PC):
pc=A*Vp(:,iof);

pr=Lmb(iof)/sum(Lmb)*100;

cl2 = colormap_red(100);
cl1 = colormap_blue(100);
for ik=1:2;
  cl2(ik,:) = [1 1 1];
  cl1(ik,:) = [1 1 1];
end
cl1 = flipud(cl1);
cmp = [cl1;cl2];
cmp = smooth_colormap(cmp,10);

c1=-0.04;
c2=0.04;

switch(regn),
 case('natl');
  xl1=45;
  xl2=122;
  yl1=12;
  yl2=112;
  cr1=-0.05;
  cr2=0.05;
  dcr=0.005;
 case('arctic');
  xl1=35;
  xl2=168;
  yl1=70;
  yl2=210;
  cr1=-0.05;
  cr2=0.05;
  dcr=0.005;
 case('arctB');
  xl1=35;
  xl2=168;
  yl1=70;
  yl2=210;
  cr1=-0.05;
  cr2=0.05;
  dcr=0.005;
end

  
% Spatial filtering
Hmsk = HH*0;
Hmsk(Iocn) = 1;
pgrd=9;
Eof = sub_fltr(Eof,pgrd,Hmsk);
%dmm = sub_fltr(dmm,pgrd,Hmsk);

Hmsk = HH*0;
Hmsk(HH>0)=1;
cmm = [1,1,1; 0,0,0];

dmm = LN;

figure(1); clf;
axes('Position',[0.05 0.1 0.8 0.8]);
pcolor(Hmsk); shading flat;
colormap(cmm);
caxis([0 1]);
hold on;
freezeColors;

%contour(HH,[0 0],'k','linewidth',2);
%hold on;
pcolor(Eof); shading flat;
colormap(cmp);
caxis([c1 c2]);

if ~isempty(cr1);
  contour(Eof,[cr1:dcr:cr2],'k','linewidth',1);
end
LN(LN>140) = nan;
LN(LN<-140)= nan;
LN(LT>89)  = nan;
contour(LN,[-140:20:140],'color',[0.8 0.8 0.8],'linewidth',1);
LN = dmm;
LN(LN<0) = LN(LN<0)+360;
LN(LT>89)  = nan;
LN(LN>360) = nan;
LN(LN<140) = nan;
contour(LN,[160:20:200],'color',[0.8 0.8 0.8],'linewidth',1);
contour(LT,[40:10:89],'color',[0.8 0.8 0.8],'linewidth',1);
axis('equal');
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[],...
	'ytick',[]);

ch = colorbar;
set(ch,'position',[0.86 0.1 0.018 0.8],...
       'TickLength',0.022);

ttl=sprintf('EOF-%i, %i%%, R=%ikm, maxEOF=%6.4f, CFSR/CFSv2',...
	    iof,round(pr),RR,max(max(Eof)));
title(ttl);
bottom_text(btx,'pwd',1,'position',[0.02 0.075 0.8 0.08]);

if s_fig
  fgnm=sprintf('%s%s_eof%i_vort%i',pthfig,regn,iof,RR);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r250',fgnm);
end;


DV=datevec(TM);
Is=find(DV(:,2)>3 & DV(:,2)<10); % summer
Iw=find(DV(:,2)>=10 | DV(:,2)<=3);

pcW=mean(pc(Iw));
pcW90=prctile(pc(Iw),90);
pcS=mean(pc(Is));
pcS90=prctile(pc(Is),10);

nrc = length(TM); 
yrs = [0:nrc-1]/12+DV(1,1);
% Annual mean PC:
for yr=1993:2016
  I=find(DV(:,1)==yr);
  pcY(I)=mean(pc(I));
end


% PC time series:
Wn = 1/6;
[Bf,Af] = butter(9,Wn,'low');
pcf = filtfilt(Bf,Af,pc); % m3/s

figure(2); clf;
axes('Position',[0.08 0.52 0.85  0.4]);
plot(yrs,pc,'linewidth',2);
hold
plot(yrs,pcY,'r-');
plot([yrs(1) yrs(end)],[0 0],'k-','Linewidth',1);
%plot([yrs(1) yrs(end)],[pcW pcW],'r-','Linewidth',1);
%plot([yrs(1) yrs(end)],[pcW90 pcW90],'r--','Linewidth',1);
%plot([yrs(1) yrs(end)],[pcS pcS],'r-','Linewidth',1);
%plot([yrs(1) yrs(end)],[pcS90 pcS90],'r--','Linewidth',1);
stt{1}=sprintf('Summer PC=%d',pcS);
stt{2}=sprintf('Winter PC=%d',pcW);
%text(yrs(10),max(pc),stt,'Color',[1 0 0],'Fontsize',12);
yll=1.1*max(abs(pc));

set(gca,'xlim',[yrs(1) yrs(end)],...
	'ylim',[-yll yll],...
	'tickdir','out',...
	'xtick',[1993:2017],...
	'xgrid','on',...
	'xminortick','on',...
	'ygrid','on',...
	'Fontsize',12);

%datetick('x','mm/yy','keeplimits','keepticks');
set(gca,'Fontsize',12);
ttl=sprintf('%s PC-%i, %i%%, R=%ikm',regn,iof,round(pr),RR);
title(ttl);

% Plot first N eigenvalues
% Do some statistics following
% Jonsson, 1991, Wind stress curl over the Nordic Sea
% Eig. values have to be > 2 dltL 
% (called "sampling errors of eigenvectors") in order to 
% form separated eigenvectors
% otherwise they are degenerate
% Note that the number of points is actually the number of independent
% points, estimated after 2-3 months - not correlated
M=round(size(VR,1)/3); % # of points in 1 time series
dltL=Lmb*((2/M)^0.5);
Nb=10;
xx=(1:Nb);

axes('Position',[0.1 0.05 0.12 0.35]);
hb=bar(xx,Lmb(1:Nb),0.8);
set(hb,'EdgeColor','none','FaceColor',[0.8 0.8 0.8]);
hold on;
for ix=1:Nb
  dl=dltL(ix);
  y1=Lmb(ix)-dl;
  y2=Lmb(ix)+dl;
  plot([ix ix],[y1 y2],'k','Linewidth',1.2);
  if ix==1, ylm1=y2; end;
end


set(gca,'tickdir','out',...
	'xlim',[0.1 Nb+0.7],...
	'xtick',[0:Nb],...
	'ylim',[0 1.05*ylm1]);
ttl=sprintf('First %i Eigenvalues with Sampling Errors',Nb);
htt=title(ttl);
set(htt,'Position',[37.5 ylm1 1]);

bottom_text(btx,'pwd',1,'position',[0.35 0.2 0.6 0.1]);

if s_fig
  fgnm2=sprintf('%s%s_eofPC%i_vort%i',pthfig,fld,iof,RR);
  fprintf('Saving %s\n',fgnm2);
%    print('-dpng','-r150',fgnm);
  print('-depsc2',fgnm2);
end;
