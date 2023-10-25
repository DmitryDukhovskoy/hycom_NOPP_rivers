% Extract / calculate monthly mean T/S profiles
% vertical distribution along specified contour -
% ~800 m isobath around Greenland
% Calculate monthly means for given year
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

%pfld = 'salin'; % field to extract
%pfld = 'temp';
YR1 = 2012;
YR2 = 2016;
%expt = 110; % experiment without runoff
expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
regn = 'ARCc0.08';
s_mat = 0; % =0 - extract data, no save; =1 extract & save, =2 - load saved

rg=9806;  % convert pressure to depth, m
hgg=1e20; 
btx = 'TSprofile_contour_greenl008.m';


fprintf('Extracting T&S, expt %3.3i, mean %i-%i\n',expt,YR1,YR2);

if s_mat==0
  fprintf(' \n   DATA NOT SAVED  \n');
  fprintf(' \n   DATA NOT SAVED  \n');
  fprintf(' \n   DATA NOT SAVED  \n');
end


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

%fgr=sprintf('%sARCc008_Topo11_GreenlandContour.mat',pthmat);
%save(fgr,'GC');

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
TSGR = struct;
TSGR.Info = 'T/S profiles along the Greenland Contour';
TSGR.GrCntr_II = GC.cntr_Iindx;
TSGR.GrCntr_JJ = GC.cntr_Jindx;
TSGR.Hbottom   = GC.Hbottom;
TSGR.DistCntr  = GC.Distance_m;

mold  = 0;
yrold = YR1;
dday  = 10; 
for iyr=YR1:YR2
  cc=0;
  TM = [];
  yr = iyr;

  fmat2 = sprintf('%s%3.3i_Greenl_TScontour_%i.mat',...
		   pthmat,expt,iyr);
  
%  for iday = 311:dday:365
  for iday = 1:dday:365
    pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    if expt==112
      pthbin=sprintf('/nexsan/hycom/ARCc0.08_112/data/%4.4i/',yr);
    end
    
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
      nlr = 41;
      fmat = sprintf('%s%3.3i_Greenl_TScontour_%i.mat',...
		   pthmat,expt,yrold);
      TSGR = sub_io_TS_month(TSGR, s_mat, mold, imo, fmat, nlr, Hs);
      mold = imo;
      yrold = yr;
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
%    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
%    F(F>hgg)=0;
%    U=F;
%
%    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
%    F(F>hgg)=0;
%    V=F;
%
    fld='thknss';
    [F,n,m,l] = read_hycom(fina,finb,fld);
%    [F,n,m,l] = read_hycom(fina,finb,fld,'r_layer',34);
    F(F>1e18)=0;
    F=F/rg;
    F(F<1e-2)=0;
    dH = F;
    
% Process segments
    FSGM = sub_TS_sgm(GC,S,T,dH,DX,DY,HH);
    Tsgm = FSGM.Tsgm;
    Ssgm = FSGM.Ssgm;
    ZZ   = FSGM.ZZ;
    DZ   = FSGM.DZ;
    
    fprintf(':::  Min T= %6.4f, Max T=%6.4f\n',min(min(Tsgm)),max(max(Tsgm)));
    fprintf(':::  Min S= %6.4f, Max S=%6.4f\n',min(min(Ssgm)),max(max(Ssgm)));
    fprintf('1 day processed %6.3f min\n\n',toc/60);
    
    TSGR(imo).nrec   = TSGR(imo).nrec+1;
    dmm              = TSGR(imo).T;
    TSGR(imo).T      = dmm+Tsgm.*DZ;
    dmm              = TSGR(imo).S;
    TSGR(imo).S      = dmm+Ssgm.*DZ;
    dmm              = TSGR(imo).ZZ;
    TSGR(imo).ZZ     = dmm+ZZ;
    dmm              = TSGR(imo).DZ; % layer thicknesses
    TSGR(imo).DZ     = dmm+DZ;
    TSGR(imo).TM     = dnmb;
%keyboard
  end  % day
    
end   % year

fmat = sprintf('%s%3.3i_Greenl_TScontour_%i.mat',...
		   pthmat,expt,yrold);
TSGR = sub_io_TS_month(TSGR, s_mat, mold, imo, fmat, nlr, Hs);




f_pflx = 0;
if f_pflx==1
  im=9; 

  dx = TSGR(1).DistCntr*1e-3; % m->km
  Hs = TSGR(1).Hbottom;
  plot(dx,Hs); % plot Bottom profile along contour

  ar = ZZ*0;
  [Dst,dmb] = meshgrid(dx,[1:nlr]);
  
  T  = TSGR(im).T;
  S  = TSGR(im).S;
  TM = TSGR(im).TM;
  dv = datevec(TM);
  ZZ = TSGR(im).ZZ;

  tstr=sprintf('arc008_%3.3i, T %4.4i/%2.2i',expt,dv(1:2));
% For plotting add surface layer      
  ZZp = [ZZ(1,:);ZZ];
  ZZp(1,:)=0;
  Tp  = [T(1,:); T];
  Dstp = [Dst(1,:);Dst];
  c1=-2;
  c2=2;
  fnmb = 1;
  pfld='temp';
  sub_plot_TS_Zcntr(Tp,ZZp,Dstp,fnmb,tstr,Hs,c1,c2,pfld);

  Sp  = [S(1,:); S];
  c1=31;
  c2=35;
  fnmb = 2;
  tstr=sprintf('arc008_%3.3i, S %4.4i/%2.2i',expt,dv(1:2));
  pfld='salin';
  sub_plot_TS_Zcntr(Sp,ZZp,Dstp,fnmb,tstr,Hs,c1,c2,pfld);
end




