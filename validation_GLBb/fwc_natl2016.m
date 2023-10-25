% Calculate FW content for monthly 2016
% to create annual mean
% Sref
% integrating to the depth of Sref
% Following Haine
% across straits 
% /nexsan/GLBa0.08/expt_90.9/data
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 1; % 
s_fig  = 0;


rg=9806;  % convert pressure to depth, m
Sref=35; % N.Atl. is too saline
TV = '07';

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % rean. 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % rean. 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
fmat = sprintf('%sFWC_natl_2016monthly.mat',pthmat);

YRPLT=[];
cnc=0;
for yr=2016:2016
  cnc=cnc+1;
  YRPLT(cnc,1)=yr;
end

f_getH=0; % prepare sub-region
if f_getH>0
  expt=190;
  sub_get_GLBgrid(pthglb,expt,ftopo);
end
fprintf('Loading topo %s\n',ftopo);
load(ftopo);
% Mask of region of interest:
LMSK=HH;
LMSK(LMSK>=-500)=0;
LMSK(LMSK<0)=1;
LMSK(1:100,1:100)=0;
LMSK(:,1:58)=0;
LMSK(461:end,1:146)=0;
LMSK(842:end,:)=0;
LMSK(800:end,500:end)=0;
LMSK(:,900:end)=0;
LMSK(522:630,233:563)=0;
IDEEP=find(LMSK>0);
INAN =find(LMSK==0);

fprintf('FWC calculation, N.Atl, save mat=%i\n',s_mat);


%np=length(YRPLT);
%iy1=1;
cc=0;
for im=1:12
  tic;
  yr=2016;
  fprintf('Calculating FWC for yr=%i\n',yr);
  
  if im<4
    expt = 911;
  else
    expt = 912;
  end
  pthbin=sprintf('/nexsan/GLBa0.08/expt_%3.1f/data/meanstd/',expt/10);
  
  fnm=sprintf('%s%3.3i_archMN.%4.4i_%2.2i',pthbin,expt,yr,im);
  fina=[fnm,'.a'];
  finb=[fnm,'.b'];
  
  ind1=IND.i1;
  ind2=IND.i2;
  jnd1=IND.j1;
  jnd2=IND.j2;
  
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  DP=F(:,jnd1:jnd2,ind1:ind2)./rg; % subset
  DP(DP<0.1)=nan; % 0-m layers, vanished

% Interf. depths and depths of the middle of the layers (m)
% NOTE: sign convention: depths are negative
% dP - Pa or m
  [ZZ,ZM] = sub_thck2dpth(DP); % 

  fld='salin';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  S=F(:,jnd1:jnd2,ind1:ind2);


% Follow ~Haine et al., 2015 definition of FWC
% The better way would be to interpolate
% between S values to find the exact depth
%  ZZi=[(0:-1:-20)';(-22:-2:-50)';(-60:-10:-500);(-550:-50:-1000)'];
  [ll,mm,nn]=size(S);
  Fwc=zeros(mm,nn);
%keyboard
%  tic;
  for k=1:ll-1
    dz=abs(squeeze(DP(k,:,:)));
    dz(isnan(dz))=0;
    ss1=squeeze(S(k,:,:));
    ss2=squeeze(S(k+1,:,:));
    ss1(INAN)=nan;
    ss2(INAN)=nan;
    Ib1=find(ss1>Sref);
    if k==1, Isrf=Ib1; end; % exclude region with surf S>Sref
    ss1(Ib1)=nan;
    ss2(Ib1)=nan;
    Ib2=find(ss2>Sref); % CHeck next layer see if S > Sref there
 % Interpolate to find depth of Sref
    fwc=dz.*(Sref-ss1)/Sref; % m
% Take care about layers where S becomes S>Sref
% add FWC below layer interface from depth of Sref
% up to the bottom lyaer k, S=0.5(S(k)+Sref)
    if ~isempty(Ib2)
      zm1 = squeeze(ZM(k,:,:));
      zm2 = squeeze(ZM(k+1,:,:));
      dZm = zm2-zm1;
      dS = ss2-ss1;
      zmSref = zm1+dZm./dS.*(Sref-ss1);
      zz = squeeze(ZZ(k+1,:,:)); % lower interface
      dltZs  = abs(zz-zmSref);
      Smn=0.5*(ss1+Sref); % for integration from 
      dltFwc = dltZs.*(Sref-Smn)/Sref;
      fwc(Ib2) = fwc(Ib2)+dltFwc(Ib2);
    end
    fwc(Ib1)=0;    
    Fwc=Fwc+fwc;
  end
  Fwc(Isrf)=0; 
%  Fwc(Ilnd)=nan; % exclude not needed regions
  fprintf('FWC calculation %6.2f min\n',toc/60);
  
  fprintf('%i: Subpolar Gyre, FWC=%8.2f m\n',yr,Fwc(200,300)); 
  fprintf('%i: GreenGyre, FWC=%8.2f m\n',yr,Fwc(600,650)); 
%keyboard
  cc=cc+1;
  FWC(cc).Title = 'Annual monthly FWC  from GLBb0.08';
  FWC(cc).Source=pthbin;
  FWC(cc).Year  = yr;
  FWC(cc).Fwc_m = Fwc;
  FWC(cc).Sref = Sref;
% saving data
  if s_mat>0
    fprintf('Saving %s\n\n',fmat);
    save(fmat,'FWC');
  end;
end; % months

% Calculate annual mean and add to 1993-2015
f_add=0;
if f_add>0
% load fmat;
  FWC2 = FWC;
  fmat1 = sprintf('%sFWC_natl_1993-2015.mat',pthmat);
  load(fmat1);
  smm= FWC2(im).Fwc_m*0;
  for im=1:12
    A = FWC2(im).Fwc_m;
    smm=smm+A;
  end
  A = smm/12;
  nrc=length(FWC);
  nrc=nrc+1;
  FWC(nrc).Title = FWC(nrc-1).Title;
  FWC(nrc).Source = 'fwc_natl2016.m';
  FWC(nrc).Year = 2016;
  FWC(nrc).Fwc_m  = A;
  FWC(nrc).Sref = 35;
  
  fmat2 = sprintf('%sFWC_natl_1993-2016.mat',pthmat);
  fprintf('Saving %s\n',fmat2);
  save(fmat2,'FWC');
  
end

  