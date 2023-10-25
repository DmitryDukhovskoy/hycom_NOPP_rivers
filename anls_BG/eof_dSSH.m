% Calc EOF of SSH anomaly fields
% from HYCOM simulation
% extracted in calc_dSSH.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

% Experiments:
% 110 - HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

iof = 1;
s_eof = 2; % save eof vectors, it takes a while to to SVD
s_fig=0;

regn = 'ARCc0.08';
expt = 110;  
%pthfig = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig  = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmat    = sprintf('%sLapl_SSH_BG.mat',pthmat);
feof    = sprintf('%seof_ssh.mat',pthmat);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);


[DX,DY]=sub_dx_dy(LON,LAT);
[II,JJ]=meshgrid([1:nn],[1:mm]);

hmsk=HH;
hmsk(HH<0)=nan;

% Define interior Arctic Ocean
hmin = -800;
ARC = sub_arctic_domain(HH,hmin);

fprintf('Loading %s\n',fmat);
load(fmat);

Iocn = LAPLE.Iocn;
TMd = LAPLE.TM;
%dv = 1;
VR   = double(LAPLE.SSH); % SSH demeaned

[VR,TM] = sub_day2month(VR,TMd);

% Subsample more:
% Get smaller domain - only BG region
% otherwise, EOF is too noisy and doesn't pick up SSH
% variabilty in BG very well

dii=4;
Iocn=Iocn(1:dii:end);
VR=VR(:,1:dii:end);


fprintf('Detrending\n');
%VRd=detrend(VR,1);
clear VRd
for k=1:length(Iocn);
  amm=VR(:,k);
  na=length(amm);
  X=[ones(na,1),[1:na]',([1:na].^2)'];
  [b,bint,rr]=regress(amm,X);
  VRd(:,k)=rr;
end

J0=find(isnan(VRd));
VRd(J0)=0;

%A=VR;
if s_eof==1
  fprintf('Calculating SVD, may take a while ...\n');
  [Up,Sp,Vp] = svd(VRd);
  fprintf('Saving %s\n',feof);
% Save first 10 e/vectors
% and Lmb

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

  Vp = Vp(:,1:10);
  save(feof,'Lmb','Vp');
else
  fprintf('Loading %s\n',feof);
  load(feof); %--> Lmb, Vp(:,1:10)
end  


% e/vectors:
V = Vp(:,iof);

% Fill in missing grid points:
dv=4;
nI=length(Iocn);
Iall=ARC.IN;
nALL=length(Iall);
cc=0;
%for ii=1:nI
%  cc=cc+1;
%  if cc>nALL; break; end;
%  Vf(cc)=V(ii);
%  for jj=1:3
%    cc=cc+1;
%    if cc>nALL; break; end;
%    Vf(cc)=V(ii);
%  end
%end
clear Vf
for ii=1:nI-1
  i1=Iocn(ii);
  i2=Iocn(ii+1);
  iall1 = find(Iall==i1);
  iall2 = find(Iall==i2);
  Vf(iall1:iall2-1)=V(ii);
end
ii=nI;
iall1=iall2;
iall2=length(Iall);
Vf(iall1:iall2)=V(ii);

Eof=HH*nan;
Eof(Iall)=Vf;

% Make EOF in Irminger Sea >0 
% for ease of visual comparison
%dmm = Eof(40,80);
%if dmm<0
%  Vp=-Vp;
%  V=Vp(:,iof);
%  Eof(Iocn)=V;
%end

% Pr. Components (PC):
pc=VRd*Vp(:,iof);

pr=Lmb(iof)/sum(Lmb)*100;

cl2 = colormap_red(100);
cl1 = colormap_blue(100);
for ik=1:2;
  cl2(ik,:) = [1 1 1];
  cl1(ik,:) = [1 1 1];
end
cl1(end,:) = [0 0 0];
cl2(end,:) = [0 0 0];
cl1 = flipud(cl1);
cmp = [cl1;cl2];
cmp = smooth_colormap(cmp,10);

c1=-0.012;
c2=0.012;
xl1=250;
xl2=1300;
yl1=700;
yl2=2000;
%cr1=[];
cr1=-0.016;
cr2=0.016;
dcr=0.002;
  
% Spatial filtering
Hmsk = HH*0;
Hmsk(Iall) = 1;
pgrd=15;
Eof = sub_fltr(Eof,pgrd,Hmsk);
%dmm = sub_fltr(dmm,pgrd,Hmsk);

Hmsk = HH*0;
Hmsk(HH>0)=1;
cmm = [1,1,1; 0,0,0];

LN = LON;
LT = LAT;
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

ttl=sprintf('EOF-%i, %i%%, maxEOF=%6.4f, SSH',...
	    iof,round(pr),max(max(Eof)));
title(ttl);
btx = 'eof_dSSH.m';
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


% PC time series:
figure(2); clf;
axes('Position',[0.08 0.52 0.85  0.4]);
plot(yrs,pc,'linewidth',2);
hold
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
	'ytick',[-10:2:10],...
	'tickdir','out',...
	'xtick',[1993:2017],...
	'xgrid','on',...
	'xminortick','on',...
	'ygrid','on',...
	'Fontsize',12);

%datetick('x','mm/yy','keeplimits','keepticks');
set(gca,'Fontsize',12);
ttl=sprintf('%s PC-%i, %i%%',regn,iof,round(pr));
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
  fgnm2=sprintf('%s%s_eofPC%i_dSSH_BG%i',pthfig,fld,iof,RR);
  fprintf('Saving %s\n',fgnm2);
%    print('-dpng','-r150',fgnm);
  print('-depsc2',fgnm2);
end;
