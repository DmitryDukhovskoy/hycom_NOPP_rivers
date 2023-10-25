% Extract U,V and Tracers (Arctic only, no Greenland Tr.) to calculate
% across Davis Strait
% Do for 0-50m 
% and whole water column
% Collocate U,V with Tr points
% For S and T fluxes - see calc_FwFlx_straits.m
% extract_TS*m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2013;
YR2 = 2016;

expt = 110;
Zlvl = -50;  % integrate down to 50m and bottom
ntr  = 4;
hgg  = 1e20;
rg   = 9806;

s_mat = 1;

TrNm{1} = 'Greenland';
TrNm{2} = 'Mackenzie';
TrNm{3} = 'East Euras.R.';
TrNm{4} = 'West Euras.R.';
TrNm{5} = 'Pacific Water';

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

cBr = sub_BeringTrConc(0);% coefficient for Bering Str., tracer, kg/m3

fprintf('Calculating years: %i - %i\n',YR1,YR2);

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

SCT = sub_set_Davis_xsct(HH,LON,LAT);

TFLX = struct;
TFLX.Info = 'Fluxes of Tracers 2-5, Davis Strait';
TFLX.Z0_m   = Zlvl;
TFLX.Indx   = SCT.II;
TFLX.Jndx   = SCT.JJ;
TFLX.Dist_m = SCT.Hbottom;

mold  = 0;
yrold = YR1;
dday  = 15;

for iyr=YR1:YR2
  cc=0;
  TM = [];
  yr = iyr;
  
  for iday = 2:dday:365
    pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
    cc   = cc+1;
    j1d  = datenum(yr,1,1);
    dnmb = j1d+iday-1;
    DV   = datevec(dnmb);
    imo  = DV(2);
    TM(cc,1) = dnmb;
  

% ==================
% Save monthly means
% ==================
    if mold~=imo
      fmat = sprintf('%s%3.3i_TrFlx_Davis_%i.mat',...
		   pthmat,expt,yrold);
      nlr = 41;
      TFLX = sub_io_TrFlxFram_month(TFLX, s_mat, mold, imo, fmat, nlr);
      mold = imo;
      yrold = yr;
    end
   
    
    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);
    
    tic;
%    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
%    F(F>hgg)=nan;
%    T=F;

%    [F,n,m,nlr] = read_hycom(fina,finb,'salin');
%    F(F>hgg)=nan;
%    S=F;
%
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F(F>hgg)=0;
    U=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>hgg)=0;
    V=F;

    fld='thknss';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>1e18)=0;
    F=F/rg;
    F(F<1e-2)=0;
    dH = F;
    
% Process segments
    II = TFLX.Indx;
    JJ = TFLX.Jndx;
    np = length(II);
 %   FWflx = zeros(nlr,np)*nan;
    ZZ    = zeros(nlr,np)*nan;
    Vflx  = zeros(1,np); % volume flux for checking
    TrFlx = zeros(ntr,np);
    
%    FWF = sub_calc_sflx(S,U,V,dH,SCT,HH,Sref,DX,DY,Zlvl);
%    clear S

    TFLX(imo).nrec        = TFLX(imo).nrec+1;

% ----------------
% Tracrers
% ----------------
    for nTr = 2:5
      trnm = TrNm{nTr};
      iTr = nTr-1;
      fprintf('   Reading Tr# %i, %s\n',nTr,trnm);

      [F,nn,mm,ll] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
      F(F>1e6)=nan;
      
% For Pacific Water,
% Convert nondim concentration 1 to 
% actual FW flux-based conc.
% used Woodgate & Aagaard estimate of Bering FW flux
%
      if nTr==5, 
%	fprintf('Adjusting Bering Strait flux\n');
	F = F*cBr;
      end
      Ctr = F; % kg/m3
      
      srf=0;
      FTR = sub_calc_sflx(Ctr,U,V,dH,SCT,HH,srf,DX,DY,Zlvl);
      
      TFLX(imo).Trace_Name = trnm;
      dmm = TFLX(imo).TrFlx_kgs(iTr,:);
      TFLX(imo).TrFlx_kgs(iTr,:)   = dmm+FTR.FWFlx;
      dmm = TFLX(imo).TrFlxZ_kgs(iTr,:);
      TFLX(imo).TrFlxZ_kgs(iTr,:)  = dmm+FTR.FWFlxZ;
    end
% ----------------
    fprintf('1 day: %8.4f min\n\n',toc/60);
    
  end % iday
end   % year

fmat = sprintf('%s%3.3i_TrFlx_Davis_%i.mat',...
		pthmat,expt,yrold);
nlr = 41;
TFLX = sub_io_TrFlxFram_month(TFLX, s_mat, mold, imo, fmat, nlr);
