% To calculate eddy and mean fluxes <u>*<T>
% and <u'*T'> need mean U, T, S - reference means
% u_bar, T-bar, S_bar: u=u_bar+u'
% 
% along the Greenland contour
% Save monthly mean U, T 
% along the closed contour around Greenland
% Interpolated onto fixed Z-levels
% for ease of caclulations of U'
% Only normal components are saved
%
% For each segment (i,j)-(i+1,j), i.e. X-section
% only save V1=(v(i,j+1)+v(i,j))/2
%  and V2=...
% This will allow to compute flux in the middle of the segment
%
%
% Fluxes calculated:
%            v(i,j+1)        v(i+1,j+1)
%      |--------|-------|-------|---------|
%      |                |                 |
%      |                |                 |
%      |     V1         |       V2        |  Flx = Cp*rho*(T-Tref)V1*dH1*dX1/2+
%      -        *===============*         |        Cp*rho*(T-Tref)V2*dH2*dX2/2+
%u(i,j)|     T,S,dH     |                 |   where V1 and V2 are flux-averaged
%      |                |                 |   v(i,j+1)&v(i,j) and v(i+1,j+1)&v(i+1,j)
%      |                |                 |
%      |--------|-------|-------|---------|
%      |      v(i,j)         v(i+1,j)
%      |   
%
%  
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat = 1; 

YR1 = 2004;
rg  = 9806;  % convert pressure to depth, m
hgg = 1e20; 

plr=0; % highlight this interface
btx = 'ocn_hflx_greenl008.m';

fprintf('Save monthly U,T,S, Greenland Contour,  intrp==>Zlevels, %i\n',YR1);

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';
fmatu    = sprintf('%sarc08_expt%3.3i_greenl_contr_mnthUzlv_%i.mat',pthmat,expt,YR1);
fmatt    = sprintf('%sarc08_expt%3.3i_greenl_contr_mnthTzlv_%i.mat',pthmat,expt,YR1);
fmats    = sprintf('%sarc08_expt%3.3i_greenl_contr_mnthSzlv_%i.mat',pthmat,expt,YR1);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section
np = length(Hs);

% Specify fixed ZZ (interf) and ZM levels:
ZZf = [(0:-1:-10)';(-12:-2:-30)';(-35:-5:-100)';...
       (-110:-10:-1000)';(-1025:-25:-2500)';(-2550:-50:-5000)'];
kzz = length(ZZf);

dZf=diff(ZZf);
ZMf = [];
for ik=1:kzz-1
  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
end

UZGR = GC;
UZGR.Title = 'U norm, m/s, Greenl Contour, Z levels, monthly mean';
UZGR.ZZlevels = ZZf;
UZGR.ZM = ZMf;

TZGR = GC;
TZGR.Title = 'T, Greenl Contour, Z levels, monthly mean';
TZGR.ZZlevels = ZZf;
TZGR.ZM = ZMf;

SZGR = GC;
SZGR.Title = 'S, Greenl Contour, Z levels, monthly mean';
SZGR.ZZlevels = ZZf;
SZGR.ZM = ZMf;


dday = 5;
yr = YR1;
iyr = yr;
for imo=1:12
  d1=datenum(yr,imo,1);
  d2=d1+32;
  dv2=datevec(d2);
  dm=datenum(dv2(1),dv2(2),1)-d1;
  
  Usm = zeros(kzz-1,np);
  Tsm = zeros(kzz-1,np);
  Ssm = zeros(kzz-1,np);
  
  cc=0;
  for mdd=1:dday:dm
    dnmb = datenum(iyr,imo,mdd);
    DV   = datevec(dnmb);
    iday = dnmb-datenum(iyr,1,1)+1;
    pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
    
    if ~exist(fina,'file') | ~exist(finb,'file')
      fprintf('Not found %s or %s\n\n',fina,finb);
      continue;
    end
    
    cc=cc+1;

    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);
    
    tic;
    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
    F(F>hgg)=nan;
    T=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'salin');
    F(F>hgg)=nan;
    S=F;
%
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F(F>hgg)=0;
    U=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>hgg)=0;
    V=F;

%    fld='thknss';
%    [F,n,m,l] = read_hycom(fina,finb,fld);
%    [F,n,m,l] = read_hycom(fina,finb,fld,'r_layer',34);
%    F(F>1e18)=0;
%    F=F/rg;
%    F(F<1e-2)=0;
%    dH = F;
    
    
    [ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
%    ZZ(isnan(ZZ))=100;
%    ZM(isnan(ZM))=100;
    dH=abs(diff(ZZh,1));
    [Uz,Tz,Sz] = sub_interp2z(GC,U,V,T,S,ZZf,ZZh,HH,dH,DX,DY);    
    
    Usm = Usm+Uz;
    Tsm = Tsm+Tz;
    Ssm = Ssm+Sz;

    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);
    
  end  % month day
  
  Usm = Usm/cc;
  Tsm = Tsm/cc;
  Ssm = Ssm/cc;
  
  umx = max(max(Usm));
  umn = min(min(Usm));
  tmx = max(max(Tsm));
  tmn = min(min(Tsm));
  smx = max(max(Ssm));
  smn = min(min(Ssm));
  fprintf('End of Month: umx = %6.3f umn=%6.3f\n',umx,umn);
  fprintf('End of Month: tmx = %6.3f tmn=%6.3f\n',tmx,tmn);
  fprintf('End of Month: smx = %6.3f smn=%6.3f\n',smx,smn);

  UZGR.U(imo,:,:) = Usm;
  TZGR.T(imo,:,:) = Tsm;
  SZGR.S(imo,:,:) = Ssm;
  
  if s_mat==1
    fprintf('Saving %s\n',fmatu);
    save(fmatu,'UZGR');

    fprintf('Saving %s\n',fmatt);
    save(fmatt,'TZGR');

    fprintf('Saving %s\n',fmats);
    save(fmats,'SZGR');
  end  

end    % months

if s_mat==1
  fprintf('Saving %s\n',fmatu);
  save(fmatu,'UZGR');

  fprintf('Saving %s\n',fmatt);
  save(fmatt,'TZGR');

  fprintf('Saving %s\n',fmats);
  save(fmats,'SZGR');
end  
