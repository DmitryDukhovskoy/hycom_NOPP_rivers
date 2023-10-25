% Calculate FW content in the Arctic Ocean
% relative to Sref
% GLBb data were remapped onto ARCc grid
% These files include GLBb reanalysis (for 1993-2012)
% and analysis (expt 90.[01])
% see: ~ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
%  /REMAP_ARCc/remap_gridGLB2ARC_archv/remap_annualNyrs_analysisGLBa.csh
% integrating to the depth of Sref
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

pthARC = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mean/';
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
%ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
fmat = sprintf('%sFWC_ArcticOcean_1993-2015.mat',pthmat);

YRPLT=[];
cnc=0;
for yr=1993:2015
  cnc=cnc+1;
  YRPLT(cnc,1)=yr;
end

ftopo = sprintf('%sdepth_ARCc0.08_%s.nc',pthtopo,TV); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

IND = sub_get_ARCgrid;
ind1=IND.i1;
ind2=IND.i2;
jnd1=IND.j1;
jnd2=IND.j2;

HH  = HH(jnd1:jnd2,ind1:ind2);
LON = LON(jnd1:jnd2,ind1:ind2);
LAT = LAT(jnd1:jnd2,ind1:ind2);
[mm,nn]=size(HH);


% Mask of region of interest:
LMSK=HH;
LMSK(LMSK>=-5)=0;
LMSK(LMSK<0)=1;
IDEEP=find(LMSK>0);
INAN =find(LMSK==0);

fprintf('FWC calculation, Arctic Ocean, save mat=%i\n',s_mat);


np=length(YRPLT);
iy1=1;
cc=iy1-1;
for iyr=iy1:np
  tic;
  yr=YRPLT(iyr);
  fprintf('Calculating FWC for yr=%i\n',yr);
  
  pthbin=pthARC;
  if yr<1995
    expt=190;
  elseif yr>=1995 & yr <2013
    expt=191;
  elseif yr==2013
    expt=910; % 
  elseif yr==2014 | yr==2015
    expt=911; % 
  end
  
  fnm=sprintf('%s%3.3i_archMN_GLBb2ARCc.%4.4i_01_%4.4i_12',pthbin,expt,yr,yr);
  fina=[fnm,'.a'];
  finb=[fnm,'.b'];
  
  
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
  
  fprintf('%i: Beaufort Gyre, FWC=%8.2f m\n',yr,Fwc(1388,334)); 
  fprintf('%i: Eurasian Bas., FWC=%8.2f m\n',yr,Fwc(1125,814)); 
%keyboard
  cc=cc+1;
  FWC(cc).Title = 'Annual mean FWC in Arctic Ocean from GLBb0.08';
  FWC(cc).Source=pthbin;
  FWC(cc).Year  = yr;
  FWC(cc).Fwc_m = Fwc;
  FWC(cc).Sref = Sref;
% saving data
  if s_mat>0
    fprintf('Saving %s\n\n',fmat);
    save(fmat,'FWC');
  end;
end; % year
  