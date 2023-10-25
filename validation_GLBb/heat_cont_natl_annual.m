% Calculate Heat content
% relative to Sref
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
Tref=-1.9; % Reference T, C
Zref=-100; % depth for Heat Cont. calc. 
Cp = 4186; % heat capacity, J*kg/C
TV = '07';

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % rean. 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % rean. 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
fmat = sprintf('%sHeatCont%4.4im_natl_1993-2015.mat',pthmat,abs(round(Zref)));

YRPLT=[];
cnc=0;
for yr=1993:2015
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
LMSK(730:end,500:end)=0;
LMSK(:,900:end)=0;
LMSK(522:630,233:563)=0;
IDEEP=find(LMSK>0);
INAN =find(LMSK==0);

fprintf('HEAT CONT calculation, N.Atl, save mat=%i\n',s_mat);


np=length(YRPLT);
iy1=1;
cc=iy1-1;
for iyr=iy1:np
  tic;
  yr=YRPLT(iyr);
  fprintf('Calculating HTC for yr=%i\n',yr);
  
  if yr<1995
    pthbin=pth190;
    expt=190;
  elseif yr>=1995 & yr <2013
    pthbin=pth191;
    expt=191;
  elseif yr==2013
    pthbin='/nexsan/GLBa0.08/expt_91.0/data/meanstd/';
    expt=910; % 
  elseif yr==2014 | yr==2015
    pthbin='/nexsan/GLBa0.08/expt_91.1/data/meanstd/';
    expt=911; % 
  end
  
  fnm=sprintf('%s%3.3i_archMN.%4.4i_01_%4.4i_12',pthbin,expt,yr,yr);
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

  fld='temp';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  T=F(:,jnd1:jnd2,ind1:ind2);
  fld='salin';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  S=F(:,jnd1:jnd2,ind1:ind2);
  Rho=sw_dens0(S,T);

  [ll,mm,nn]=size(S);
  HC=zeros(mm,nn);
%keyboard
%  tic;
  for k=1:ll-1
    dz=abs(squeeze(DP(k,:,:)));
    dz(isnan(dz))=0;
    tt1=squeeze(T(k,:,:));
    rho1=squeeze(Rho(k,:,:));
    zz1=squeeze(ZZ(k,:,:));
    IZ=find(zz1<=Zref);
    dz(IZ)=0;
    IT=find(tt1<Tref);
    tt1(IT)=Tref;
    hcnt=rho1.*Cp.*(tt1-Tref).*dz; % J/m2
    HC=HC+hcnt;
  end

%  Fwc(Ilnd)=nan; % exclude not needed regions
  fprintf('HEATC calculation %6.2f min\n',toc/60);
  
  fprintf('%i: Subpolar Gyre, HEATC=%8.2f J/m2\n',yr,HC(200,300)); 
  fprintf('%i: GreenGyre, HTC=%8.2f J/m2\n',yr,HC(600,650)); 
%keyboard
  cc=cc+1;
  HTC(cc).Title = 'Annual mean Heat Cont  from GLBb0.08';
  HTC(cc).Source=pthbin;
  HTC(cc).Year  = yr;
  HTC(cc).HeatCnt_J_m2 = HC;
  HTC(cc).Zref = Zref;
  HTC(cc).Tref = Tref;
% saving data
  if s_mat>0
    fprintf('Saving %s\n\n',fmat);
    save(fmat,'HTC');
  end;
end; % year
  