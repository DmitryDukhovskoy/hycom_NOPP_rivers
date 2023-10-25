% Calculate ocean heat flux to Greenland
% across specified contour -
% isobath around Greenland
% save daily flux
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

YR1 = 2005;
YR2 = 2005;
Cp = 4200; % J/kg K
%Tref= -273.15; % Ref T to calc. H flux
Tref= -1.8; % Ref T to calc. H flux

expt = 110;

s_mat = 1; % =0 - extract data, no save; =1 extract & save, =2 - load saved


rg=9806;  % convert pressure to depth, m
hgg=1e20; 

plr=0; % highlight this interface
btx = 'ocn_hflx_greenl008.m';

fprintf('Oceanic Heat Flux, Greenland Section, %i-%i\n',YR1,YR2);

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthfig=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_green_xsct/',...
		  regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section

% ================================================
% Plot Greenland map and the contour
% ================================================
f_pltgr=0;
if f_pltgr==1
  fn=10;
  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  
  bottom_text(btx,'pwd',1);
end

% 
if s_mat==1
  fprintf('Mat file will be saved %s\n',pthmat);
end



% ================================================
HFLX = struct;
HFLX.Info = 'Heat & Vol. Fluxes across 800m isobath around Greenland';
HFLX.Tref      = Tref; 
HFLX.GrCntr_II = GC.cntr_Iindx;
HFLX.GrCntr_JJ = GC.cntr_Jindx;
HFLX.Hbottom   = GC.Hbottom;
HFLX.DistCntr  = GC.Distance_m;

nlr = 41;
mold  = 0;
yrold = YR1;
dday  = 1; 
for iyr=YR1:YR2
  cc=0;
  TM = [];
  yr = iyr;
  
% Initialze
  fmat = sprintf('%s%3.3i_Greenl_HVflx_daily%i.mat',...
	       pthmat,expt,iyr);
  HFLX = sub_io_hflx_month(HFLX, 0, 0, 1, fmat, nlr);

  for iday = 1:dday:365
    pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
    cc   = cc+1;
    j1d  = datenum(yr,1,1);
    dnmb = j1d+iday-1;
    DV   = datevec(dnmb);
    imo  = DV(2);
    TM(cc,1) = dnmb;

    
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

    fld='thknss';
    [F,n,m,l] = read_hycom(fina,finb,fld);
%    [F,n,m,l] = read_hycom(fina,finb,fld,'r_layer',34);
    F(F>1e18)=0;
    F=F/rg;
    F(F<1e-2)=0;
    dH = F;
    
    
%    [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
%    ZZ(isnan(ZZ))=100;
%    ZM(isnan(ZM))=100;

    hflg = 1; % weight U correct on heat flux bias in deep layers
    VCT = sub_crrct_vflx(GC,U,V,dH,DX,DY,HH,S,T,Tref,hflg);
%    Vct = 0;
% Process segments
    FSGM = sub_vhFlx_sgm(GC,U,V,S,T,dH,DX,DY,HH,VCT,Tref);
    Vflx = FSGM.Vflx_Dintgr_m3_s;
    Hflx = FSGM.Hflx_grcell_W;
    ZZ   = FSGM.ZZ;
    
   fprintf('VFlux=%10.8d, HFlux=%10.8d\n',nansum(Vflx), nansum(nansum(Hflx)));
    
    fprintf('1 day processed %6.3f min\n\n',toc/60);

    HFLX(cc).TM           = TM;
    HFLX(cc).Vol_flux_m3s = Vflx;
    HFLX(cc).Hflux_W      = Hflx;
    HFLX(cc).ZZ           = ZZ;
    HFLX(cc).TM           = dnmb;
%keyboard
  end  % day

  if s_mat==1    
    fprintf('Saving %s\n',fmat);
    save(fmat,'HFLX');
  end
  
end   % year




