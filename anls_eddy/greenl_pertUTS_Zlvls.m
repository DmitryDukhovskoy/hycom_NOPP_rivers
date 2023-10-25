% Calculate perturb of u', t', s' 
% for eddy flux calculations
% <u*T>=<umn>*<Tmn>+<u'*T'>
% 
% Saved fields monthly mean perturb terms:
% <u'*u'>, <u'*T'>, <u'*S'> and total <u*T>
%
% along Greenland contour
% see flux calculation description
% in greenl_meanUTS_Zlvls.m
% means are taken from 
% greenl_meanUTS_Zlvls.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat = 1; % =1 - start from Jan, =2 - start from last saved 

YR1 = 2005;
rg  = 9806;  % convert pressure to depth, m
hgg = 1e20; 

plr=0; % highlight this interface
btx = 'ocn_hflx_greenl008.m';

fprintf('Calculate daily perturb U,T,S, Greenland Contour, on Zlevels, %i\n',YR1);

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers

pthtopo= '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';
fmatu  = sprintf('%sarc08_expt%3.3i_greenl_contr_mnthUzlv_%i.mat',pthmat,expt,YR1);
fmatt  = sprintf('%sarc08_expt%3.3i_greenl_contr_mnthTzlv_%i.mat',pthmat,expt,YR1);
fmats  = sprintf('%sarc08_expt%3.3i_greenl_contr_mnthSzlv_%i.mat',pthmat,expt,YR1);
fmatup = sprintf('%sarc08_expt%3.3i_greenl_contr_prtrbUzlv_%i.mat',pthmat,expt,YR1);
fmattp = sprintf('%sarc08_expt%3.3i_greenl_contr_prtrbTzlv_%i.mat',pthmat,expt,YR1);
fmatsp = sprintf('%sarc08_expt%3.3i_greenl_contr_prtrbSzlv_%i.mat',pthmat,expt,YR1);
fmatut = sprintf('%sarc08_expt%3.3i_greenl_contr_meanUT_%i.mat',pthmat,expt,YR1);

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

% Load mean fields
%fprintf('Loading mean fields ...\n');
%load(fmatu);
%load(fmats);
%load(fmatt);

UpZGR = GC;
UpZGR.Title = 'U norm, m/s, Greenl Contour, Z levels, perturbation';
UpZGR.ZZlevels = ZZf;
UpZGR.ZM = ZMf;

TpZGR = GC;
TpZGR.Title = 'T, Greenl Contour, Z levels, perturbation';
TpZGR.ZZlevels = ZZf;
TpZGR.ZM = ZMf;

SpZGR = GC;
SpZGR.Title = 'S, Greenl Contour, Z levels, perturbation';
SpZGR.ZZlevels = ZZf;
SpZGR.ZM = ZMf;

UTZGR = GC;
UTZGR.Title = 'total <u*T>, Greenl Contour, Z levels, monthly';
UTZGR.ZZlevels = ZZf;
UTZGR.ZM = ZMf;

dday = 5;
yr = YR1;
iyr = yr;
imo1 = 1;

if s_mat==2
  fprintf('Starting from last saved month in %s\n',fmatup);
  fprintf('Starting from last saved month in %s\n',fmattp);
  fprintf('Starting from last saved month in %s\n',fmatsp);
  fprintf('Starting from last saved month in %s\n',fmatut);

  load(fmatup);
  load(fmattp);
  load(fmatsp);
  load(fmatut);

  TM=TpZGR.TM;
  tm1=TM(end);
  TM=SpZGR.TM;
  tm2=TM(end);
  TM=UpZGR.TM;
  tm3=TM(end);
  TM=UTZGR.TM;
  tm4=TM(end);
  
  
  dmn=min([tm1,tm2,tm3,tm4]);
  dv=datevec(dmn);
  
  imo1=dv(2)+1;
  fprintf('Last saved month: %i/%i\n',dv(1),dv(2));
  fprintf('Start from month: %i/%i\n',dv(1),imo1);
  
end

for imo=imo1:12
  d1=datenum(yr,imo,1);
  d2=d1+32;
  dv2=datevec(d2);
  dm=datenum(dv2(1),dv2(2),1)-d1;
% 
% Currently Updated monthly means - will need
% to take this out of the loop later
  fprintf('Loading mean fields ...\n');
  load(fmatu);
  load(fmats);  
  load(fmatt);

  Umn=squeeze(UZGR.U(imo,:,:));
  Tmn=squeeze(TZGR.T(imo,:,:));
  Smn=squeeze(SZGR.S(imo,:,:));
  
  cc=0;  % counter monthly
  Upr=[];
  Spr=[];
  Tpr=[];
  UUm=[];
  TTm=[];
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
    
    upr = Uz-Umn;
    tpr = Tz-Tmn;
    spr = Sz-Smn;

    cc=cc+1;
    Upr(cc,:,:)=upr;
    Tpr(cc,:,:)=tpr;
    Spr(cc,:,:)=spr;
    UUm(cc,:,:)=Uz; % total U
    TTm(cc,:,:)=Tz; % total T
    
    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);
    
  end  % month day

  utmx = max(max(upr.*tpr));
  utmn = min(min(upr.*tpr));
  usmx = max(max(upr.*spr));
  usmn = min(min(upr.*spr));
  fprintf('End of Month <u''*t''>: max = %6.3f min=%6.3f\n',utmx,utmn);
  fprintf('End of Month <u''*s''>: max = %6.3f min=%6.3f\n',utmx,utmn);

  dmm=Upr.*Upr;
  UpZGR.UpUp(imo,:,:) = squeeze(nanmean(dmm,1)); % save monthly <u'*u'>
  dmm=Upr.*Tpr;
  TpZGR.UpTp(imo,:,:) = squeeze(nanmean(dmm,1)); % save monthly <u'*T'>
  dmm=Upr.*Spr; % 
  SpZGR.UpSP(imo,:,:) = squeeze(nanmean(dmm,1)); % save monthly <u'*S'>
  dmm=UUm.*TTm;  % total U and T = <u*T>
  UTZGR.UmTm(imo,:,:) = squeeze(nanmean(dmm,1)); % <u*T> monthly
  UpZGR.TM(imo,1)=dnmb;
  TpZGR.TM(imo,1)=dnmb;
  SpZGR.TM(imo,1)=dnmb;
  
  if s_mat>0
    fprintf('Saving %s\n',fmatup);
    save(fmatup,'UpZGR');

    fprintf('Saving %s\n',fmattp);
    save(fmattp,'TpZGR');

    fprintf('Saving %s\n',fmatsp);
    save(fmatsp,'SpZGR');
    
    fprintf('Saving %s\n',fmatut);
    save(fmatut,'UTZGR');
  end  
  
end    % months

if s_mat==-99999
  fprintf('Saving %s\n',fmatup);
  save(fmatup,'UpZGR');

  fprintf('Saving %s\n',fmattp);
  save(fmattp,'TpZGR');

  fprintf('Saving %s\n',fmatsp);
  save(fmatsp,'SpZGR');
  
  fprintf('Saving %s\n',fmatut);
  save(fmatut,'UTZGR');
end  
