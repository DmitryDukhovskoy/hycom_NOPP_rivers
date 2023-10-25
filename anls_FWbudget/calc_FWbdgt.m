% Caculate FW content wrt to Sref=34.8
% full depth and to isohaline = Sref (Haine et al, 2015)
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 110;
s_mat  = 2; % overwritten by s_extr; =2 - skip existing years
s_fig  = 0;
%s_extr = 1; %=0 - plot saved FWC, =1 - calculate FWC

rg=9806;  % convert pressure to depth, m
Sref=34.8;

regn = 'ARCc0.08';

%pthbin  = '/nexsan/GLBb0.08/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/GLBb0.08/fig_fwc/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

fmat    = sprintf('%sfwc_hycom_arc08_%3.3i.mat',pthmat,expt);

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);

% Mask of the region of interest
% exclude Pacific Ocean and North.Atl.
% Mask of region of interest:
LMSK=HH;
LMSK(1:200,:)   = 0;
LMSK(1935:end,:)= 0;
LMSK(LMSK>=-500)=0;
LMSK(LMSK<-500)=1;
IDEEP=find(LMSK>0);
INAN =find(LMSK==0);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% Define indices of the select regions:
% See anls_Greenland sub_define_boxes
FWBX = sub_define_boxesAO_Natl(HH,LON,LAT,0);
nbx = length(FWBX);

% Monthly values
YRPLT=[];
cc=0;
for yr=1993:2016
  for mo=1:12
    cc=cc+1;
    YRPLT(cc,1)=yr;
    YRPLT(cc,2)=mo;
  end
end
nrec=cc;


%Start time loop:
cc  = 0;
cyr = 0;
Fwc1Yr = HH*0;
Fwc2Yr = HH*0;
yre = 0;
if s_mat==2
  fprintf('Loading saved %s\n',fmat);
  load(fmat);
  yre=FWC(end).Year;
  cyr=length(FWC);
  cc=length(FWBX(1).TM);
end


for it=1:nrec
  tic;
  yr = YRPLT(it,1);
  mo = YRPLT(it,2);
  dnmb = datenum(yr,mo,15);
  iday = dnmb-datenum(yr,1,1)+1;
  
  if yr<=yre,
    fprintf('Skipping saved %i/%2.2i\n',yr,mo);
    continue;
  end
  
  
  fprintf('Extracting %i/%i\n',yr,mo);
  
  pthbin  = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',...
		    expt,yr);
  
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
% Get layer thickness to be able
% to construct depth arrays of model layers
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  DP=F./rg;
  DP(DP<0.1)=nan; % 0-m layers, vanished
% Interf. depths and depths of the middle of the layers (m)
% NOTE: sign convention: depths are negative
% dP - Pa or m
  [ZZ,ZM] = sub_thck2dpth(DP); % 

  fld='salin';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  S=F;

% Calculate FWC using 2 approaches: whole depth integrated
% and integrated from depth of Sref
%
% Follow ~Haine et al., 2015 definition of FWC
% The better way would be to interpolate
% between S values to find the exact depth
%  ZZi=[(0:-1:-20)';(-22:-2:-50)';(-60:-10:-500);(-550:-50:-1000)'];
  [ll,mm,nn]=size(S);
  Fwc1 = zeros(mm,nn);
  Fwc2 = zeros(mm,nn);
%keyboard
%  tic;
  for k=1:ll-1
    dz=abs(squeeze(DP(k,:,:)));
    dz(isnan(dz))=0;
    ss1=squeeze(S(k,:,:));
    ss2=squeeze(S(k+1,:,:));
    ss1(INAN)=nan;
    ss2(INAN)=nan;
% First method: simply depth integrate    
    fwc1 = dz.*(Sref-ss1)/Sref; % m
    Fwc1 = Fwc1+fwc1;
% 2nd method - following Haine et al    
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
    Fwc2=Fwc2+fwc;
  end
  Fwc2(Isrf)=0; 
%  Fwc(Ilnd)=nan; % exclude not needed regions

%keyboard
  Fwc1Yr = Fwc1Yr+Fwc1;
  Fwc2Yr = Fwc2Yr+Fwc2;
  if mo==12
    Fwc1Yr=Fwc1Yr/12; % annual mean
    Fwc2Yr=Fwc2Yr/12;
    
    cyr=cyr+1;
    FWC(cyr).Year=yr;
    FWC(cyr).Sref=Sref;
    FWC(cyr).FWC_DpthIntgr_anl = Fwc1Yr;
    FWC(cyr).FWC_Z_Sref_anl    = Fwc2Yr;
    Fwc1Yr = HH*0;
    Fwc2Yr = HH*0;
  end    

  cc=cc+1;
  for kb=1:nbx
    INP = FWBX(kb).IN_polygon;
    FWBX(kb).Sref    = Sref;
    FWBX(kb).TM(cc,1)= dnmb;
    FWBX(kb).Title   = sprintf('MoMean FWC, ARCc0.08-%i, %i/%i',expt,yr,mo);
    FWBX(kb).Fwc1_m(cc) = nansum(Fwc1(INP).*Acell(INP))./sum(Acell(INP));
    FWBX(kb).Fwc2_m(cc) = nansum(Fwc2(INP).*Acell(INP))./sum(Acell(INP));
  end
  
  fprintf('FWC calculation %6.2f min\n',toc/60);
  
  fprintf('%i: Beaufort Gyre, FWC1=%6.1f, FWC2=%6.1fm\n',...
	  yr,FWBX(7).Fwc1_m(cc),FWBX(7).Fwc2_m(cc)); 
  fprintf('%i: LabrdSea, FWC1=%6.1f, FWC2=%6.1fm\n',...
	  yr,FWBX(2).Fwc1_m(cc),FWBX(2).Fwc2_m(cc)); 
end

% saving data
if s_mat>0
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'FWBX','FWC');
end;
  
  




