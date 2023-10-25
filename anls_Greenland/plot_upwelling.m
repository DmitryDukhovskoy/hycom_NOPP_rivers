% Upwelling index is calculated using Bakun method:
% M = 1/(f*rho_w)*(tau x k)
% see anls_atm/calc_upwelling.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

rhoa = 1.2;
Cd   = 0.0013;
rhow = 1025;

s_mat = 1;

pthdat1 = '/Net/kronos/ddmitry/ncep_cfsr/';
pthdat2 = '/Net/kronos/ddmitry/ncep_cfsv2/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';

pthout  = ('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/');
fmat    = sprintf('%sarc008_Greenl_upwl_CFSR_month.mat',pthout);
fprintf('Loading %s\n',fmat);
load(fmat); % <--- GRC - contour, UPWL - upw. index

% Topo for ARCc0.08
ftopo = sprintf('%s/depth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LN = nc_varget(ftopo,'Longitude');
LT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);

% Load Greenland coast
% with local normal unit vectors
%GRC = sub_greenl_coast08(HH);
%ngr = length(GRC.X);

TM = UPWL.TM;
nc = length(TM);

% Select pnts where plot 
% hist of upwelling index
% averaged by months and over N pnts along contour
IJ=[   607   831
   532   535
   577   433
   657   492
   738   559
   838   610
   922   838];
npp=size(IJ,1);
Nav=40;

Xc  = GRC.Iocn;
Yc  = GRC.Jocn;
clear IP
for ip=1:npp
  i0=IJ(ip,1);
  j0=IJ(ip,2);
  D=sqrt((Xc-i0).^2+(Yc-j0).^2);
  I=find(D==min(D),1);
  IP(ip,1)=I;
  IP(ip,2)=I-round(Nav/2);
  IP(ip,3)=I+round(Nav/2);
  plot(Xc(IP(ip,2)),Yc(IP(ip,2)),'c*');
  plot(Xc(IP(ip,3)),Yc(IP(ip,3)),'m*');
  text(i0,j0,sprintf('%2.2i',ip));
end

MNRM    = UPWL.Normal_m3s; % Ekman transport
[nt,np] = size(MNRM);
nyr = nt/12;
A=reshape(MNRM,[12,nyr,np]);

btx='plot_upwelling.m';

for ipp=1:npp
  i1=IP(ipp,2);
  i2=IP(ipp,3);
  aa=A(:,:,i1:i2);
  aa=nanmean(aa,3);
  
  clear BB
  for im=1:12
    b=aa(im,:);
%    bm=mean(b);
%    bmd=median(b);
%    b1=prctile(b,10);
%    b2=prctile(b,90);
    BB(:,im)=b';
  end
  figure(ipp); clf;
  axes('Position',[0.09 0.4 0.85 0.5]);
  plot([1 12],[0 0],'k--');
  hold on
  boxplot(BB,'symbol','k.');
  ha=findobj(gca,'type','line');
  set(ha,'linew',1.6);
  yl1=min(min(BB));
  yl2=max(max(BB));

%  set(gca,'tickdir','out',...
%	  'ylim',[yl1 yl2]);
  stl=sprintf('Upw Indx, pnt %2.2i',ipp);
  title(stl);

%  bottom_text(btx,'pwd',1,'position',[0.08 0.3 0.5 0.05]);
end
bottom_text(btx,'pwd',1,'position',[0.08 0.3 0.5 0.05]);

% Mean DJF
% Plot contour
dmm=A([1,2,12],:,:);
W=squeeze(mean(dmm,1));
W=mean(W,1);


stl = sprintf('Upwelling index, CFSR/CFSv2, DJF 1993-2016, m3/s per 1m');
fn = 2;
sub_plot_upwelling(HH,LN,LT,fn,W,GRC,stl)
bottom_text(btx,'pwd',1,'position',[0.08 0.08 0.5 0.05]);

% Mean JJA
% Plot contour
dmm=A(6:8,:,:);
S=squeeze(mean(dmm,1));
S=mean(S,1);


stl = sprintf('Upwelling index, CFSR/CFSv2, JJA 1993-2016, m3/s per 1m');
fn = 3;
sub_plot_upwelling(HH,LN,LT,fn,S,GRC,stl)
bottom_text(btx,'pwd',1,'position',[0.08 0.08 0.5 0.05]);


