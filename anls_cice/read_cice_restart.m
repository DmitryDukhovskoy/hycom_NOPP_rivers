addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.08';
%R = 'ARCc0.04';

%PTH.data = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_ice_restart/';
PTH.data = '/nexsan/people/ddmitry/hycom/ARCc0.08/103/ice_restart/';
%PTH.data = '/Net/mars/ddmitry/hycom/ARCc0.04/restart/';
txtb = 'read_cice_restart.m';

switch(R),
 case('ARCc0.08');
  PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
  fltopo=sprintf('%sdepth_ARCc0.08_09.nc',PTH.topo);
%  fnmcice = 'cice.restart008_093a';
  fnmcice = 'cice.restart_piomas_ARCc0.08_T11_093a';
 case('ARCc0.04');
  PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
  fltopo=sprintf('%sdepth_ARCc0.04_17DD.nc',PTH.topo);
  fnmcice = 'cice.restart004_105a';
end

HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LAT  = alat;
LON  = elon;
[m,n]= size(HH);
[mm,nn]= size(HH);

%Pick point for test:
jj0 = 1500;
ii0 = 800;
%jj0=633;
%ii0=709;
PT.ij = [ii0,jj0];

% Fields to save:
FSV.ivel     = 0; % ice vel. components
FSV.rad      = 0; % radiation fields
FSV.ostrs    = 0; % ocean str. components
FSV.int_strs = 0; % internal stress
FSV.imask    = 0; % ice mask

%    restart dumped by HYCOM/CICE
%fmat0 = 'cice_restart093j';
%fcice = sprintf('%sarc08_T11_cice.restart_093j',PTH.data); 
%   restart created in write_cice_restart
fmat0 = 'cice_restart105a';

fcice = sprintf('%s%s',PTH.data,fnmcice); % prepared restart
%fcice = sprintf('%scice.restart_piomas_ARCc0.08_T11_093a',PTH.data); % prepared restart
%fcice = sprintf('%scice.restart_in',PTH.data);
%fcice = sprintf('%scice.restart_109i',PTH.data); % restart from CICE
%fcice = sprintf('%scice.restart_093b',PTH.data); % restart from CICE

fprintf(' ==== OPENING cice.restart: %s\n',fcice);

fidc  = fopen(fcice,'r','ieee-be');

fprintf('Opening restart file to read: %s\n',fcice);
if FSV.ivel == 1
  fprintf('Saving ice velocity into %s\n',fmat0);
end
if FSV.rad == 1
  fprintf('Saving radiation into %s\n',fmat0);
end
if FSV.ostrs == 1
  fprintf('Saving ocean stress into %s\n',fmat0);
end
if FSV.int_strs == 1
  fprintf('Saving internal stress into %s\n',fmat0);
end
if FSV.imask == 1
  fprintf('Saving ice mask into %s\n',fmat0);
end

% Constants:
L0 = 3.34e5; % J/kg, latent heat of fusion of fresh ice at 0C


% Domain size, categories, layers:
% see ice_dominan_size.F90
ncat  = 5; % # of categories
nilyr = 4; % ice layers per category
ntilyr = ncat*nilyr; % # of ice layers in all categories
nslyr = 1; % # of snow layers
ntslyr = ncat*nslyr; % # of snow layers in all cat.
ntrcr = 3; % # of tracers, 1- surf. T


% Read header
%         write(nu_dump) istep1,time,time_forc
frewind(fidc);

dmm=fread(fidc,1,'int'); % length of following recrd. byte (*8=bits)
fprintf('Reading Header\n');
istep = fread(fidc,1,'int32'); % 4byte
rday  = fread(fidc,1,'float64'); % 8, run day, sec
fday  = fread(fidc,1,'float64'); % 8
dmm=fread(fidc,1,'int'); % bites written

fprintf('Header: %i %15.4f %15.4f\n',istep,rday,fday);
fprintf('istep, rday, fday: %i %15.4f %15.4f\n',...
	istep,rday/(3600*24),fday/(3600*24));



% Do not keep fields in memory if not needed
cnt=0;
% Read for thickness categories (thickness are in hin_max)
% hin_max(n-1) < Cat n < hin_max(n)
%  0.000000000000000E+000  < Cat            1  <   0.644507216819426
%  0.644507216819426       < Cat            2  <    1.39143349757630
%   1.39143349757630       < Cat            3  <    2.47017938195989
%   2.47017938195989       < Cat            4  <    4.56728791885049
%   4.56728791885049       < Cat            5  <    9.33384181586817


% state variables:
Hmean=zeros(mm,nn);
aice = zeros(mm,nn); % aggregated ice area
for n=1:ncat
  fprintf('Category: %i min/max area, vol ice, vol snow, Tsfc\n',n);
% Concentration of ice in each category: aicen  
% Note area = fraction of grid cell for cat. n
%  aicen = (fread(fidc,[nn,mm],'float64'))'; % 
  cnt=cnt+1;
  dmm   = fread(fidc,1,'int'); % = nn*mm*8bit
%  fprintf(' dmm = %i\n',dmm);
  AA    = (fread(fidc,[nn,mm],'float64'))'; % 
  dmm   = fread(fidc,1,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: cat %i, %16.6f %16.5f\n',cnt,n,mina,maxa);
%  keyboard
  aicen=AA;
  aice=aice+aicen;
  
  PT.aicen(n) = AA(jj0,ii0);
  PT.cat(n)   = n;
%volume per unit area of ice (m)
% vicen(i,j)=aicen(i,j)*hin(i,j) 
  cnt=cnt+1;
  dmm   = fread(fidc,1,'int');
  AA    = (fread(fidc,[nn,mm],'float64'))'; 
  dmm   = fread(fidc,1,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: cat %i, %16.6f %16.6f\n',cnt,n,mina,maxa);

  PT.vicen(n) = AA(jj0,ii0);
  
  vicen=AA;
  hin = vicen./aicen; % ice thicknes, cat=n
  Hmean=Hmean+vicen;  % cell-mean ice thickness
  splt=0;
  if splt>0
    fgn=1;
    ttl=sprintf('hin: Ice Thckn,m, cat%i',n);
    sub_plot_cice(hin,fgn,ttl,HH);
  end
%  keyboard;
  
  
%volume per unit area of snow (m)
% vsnon(i,j)=aicen(i,j)*hsn(i,j) 
  cnt=cnt+1;
%  fprintf('%i: Ice cat %i, Reading snow vol/area\n',cnt,n);
  dmm   = fread(fidc,1,'int');
  AA    = (fread(fidc,[nn,mm],'float64'))';
  dmm   = fread(fidc,1,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: cat %i, %16.6f %16.6f\n',cnt,n,mina,maxa);
  
  PT.vsnon(n) = AA(jj0,ii0);
  VSNON(n,:,:) = AA;  
% Ice tracer from Tsfcn - T of 
% ice/snow top surface (C)
% trcrn 
  cnt=cnt+1;
%  fprintf('%i: Ice cat %i, Reading ice tracer\n',cnt,n);
  dmm   = fread(fidc,1,'int');
  AA    = (fread(fidc,[nn,mm],'float64'))';
  dmm   = fread(fidc,1,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: cat %i, %16.6f %16.6f \n',cnt,n,mina,maxa);
%  trcrn = AA;

  PT.trcrn(n) = AA(jj0,ii0);
  
  splt=0;
  if splt>0
    fgn=1;
    ttl='trcrn: Ice Tracer, top T';
    sub_plot_cice(AA,fgn,ttl,HH);
  end
  
end


splt=0;
if splt>0
  fgn=1;
  ttl=sprintf('Ice Thkn,m, %s %s',fnmcice,R);
  Hmean(Hmean==0)=nan;
  sub_plot_cice(Hmean,fgn,ttl,HH);
  contour(HH,[0 0],'k');
  caxis([0 3]);
  colorbar
  axis('equal');
  set(gca,'xlim',[1 nn],...
	  'ylim',[1 mm]);
%  stt=sprintf('Ice Thkn, CICE restart 1993 Jan 1');
%  title(stt);
  bottom_text(txtb,'pwd',1);
  keyboard  
end

% Plot aggregated area:
spltA = 1;
if spltA>0
  fprintf('Max Aggregated Ice are=%12.9f\n',max(max(aice)));
  I=find(aice>1.0);  
  [j,i]=ind2sub(size(aice),I);
  contour(HH,[0 0],'k')    
  hold on;
%  plot(i,j,'r.');
%  keyboard  
  
  fgn=1;
  ttl=sprintf('Aggregated Ice Area, %s %s',fnmcice,R);
  aice(aice==0)=nan;
  sub_plot_cice(aice,fgn,ttl,HH);
  contour(HH,[0 0],'k');
  caxis([0 1]);
  colorbar
  axis('equal');
  set(gca,'xlim',[1 nn],...
	  'ylim',[1 mm]);
%  stt=sprintf('Ice Thkn, CICE restart 1993 Jan 1');
%  title(stt);
  bottom_text(txtb,'pwd',1);
  keyboard  
end



clear aicen vicen vsnon trcrn
%VICEN = vicen;

%keyboard
% eicn
% Energy of melting for each ice layer (J/m2) per unit area
% (it equals to the negative enthalpie)
% This is the energy required to melt a unit volume of 
% ice and raise its T to 0C
% The ice enthalpy is a function of T and S (because
% of internal melting and freezing in
% brine pockets), S in CICE is prescribed 
% with a fixed vert. profile, so enthalpie is 
% a 1-1 relationship with temp. 
% Description - see CICE manual
fprintf('min/max eicen for each layer:\n');
% Note: sea ice layers go downward
% Layer 1 - ice/atm, Layer 4 - ice/ocean
for k=1:ntilyr
  cnt=cnt+1;
%  fprintf('%i: Tot layer %i, Reading ice melt energy\n',cnt,k);
  dmm   = fread(fidc,1,'int');
  AA    = (fread(fidc,[nn,mm],'float64'))'; % eicn
  dmm   = fread(fidc,1,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: Layer %i, %16.6f %16.6f\n',cnt,k,mina,maxa);
  PT.eicen(k) = AA(jj0,ii0);
%keyboard  
end 

clear eicen

% esnon
% Energy of melting for each snow layer (J/m2)* each category
% same as sea ice but no brain pockets problem
fprintf('min/max esnon for each layer:\n');
for k=1:ntslyr
  cnt=cnt+1;
%  fprintf('%i: Tot. snow Layer %i, Reading energy snow melt\n',cnt,k);
  dmm   = fread(fidc,1,'int');
  AA    = (fread(fidc,[nn,mm],'float64'))'; % esnon
  dmm   = fread(fidc,1,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: cat %i, %16.6f %16.6f\n',cnt,k,mina,maxa);
  PT.esnon(k) = AA(jj0,ii0);
  ESNON(k,:,:) = AA;
end 


clear esnon

%keyboard

% From enthalpie eicen
% saved in restart (J/m2) get
% T in the layer 
% Using CICE approach: ice_therm_vertical.F90:
%  qsn  (ij,k) = esnon(i,j,k)*rnslyr/vsnon(i,j)
%  where rsnlyr - real(nslyr), # of snow layers
%       !-----------------------------------------------------------------
%      ! Compute snow temperatures from enthalpies.
%      ! Note: qsn <= -rhos*Lfresh, so Tsn <= 0.
%      !-----------------------------------------------------------------
%            Tsn(ij,k) = (Lfresh + qsn(ij,k)/rhos)/cp_ice
%      ! Compute ice enthalpypacket_write_wait: Connection to 140.32.32.98: Broken pipe
%        pacity = F, qin and Tin are never used.
%      !-----------------------------------------------------------------
%            ! qin, eicen < 0
%            qin(ij,k) = eicen(i,j,k)*real(nilyr,kind=dbl_kind) &
%                        /vicen(i,j)
%
%      !-----------------------------------------------------------------
%      ! Compute ice temperatures from enthalpies using quadratic formula
%      !-----------------------------------------------------------------
%
%            if (l_brine) then
%               aa1 = cp_ice
%               bb1 = (cp_ocn-cp_ice)*Tmlt(k) - qin(ij,k)/rhoi - Lfresh
%               cc1 = Lfresh * Tmlt(k)
%               Tin(ij,k) =  (-bb1 - sqrt(bb1*bb1 - c4*aa1*cc1)) /  &
%                             (c2*aa1)
%               Tmax = Tmlt(k)
%
%            else                ! fresh ice
%               Tin(ij,k) = (Lfresh + qin(ij,k)/rhoi) / cp_ice
%               Tmax = -qin(ij,k)*puny/(rhos*cp_ice*vicen(i,j))
%                         ! as above for snow
%            endif
%
c0 = 2106;    % J/(kg*Cdeg) - specific heat of fresh ice, at 0C
L0 = 3.34e5; % J/kg  latent heat of fusion of fresh ice at 0C
mu = 0.054;  % deg/ppt - liquidus, ratio betw. Tfreez and S of brine
cw = 4218;   % J/(kg*Cdeg) - specific heat of sea water
rho_snow = 330;  % kg/m3
rho_ice  = 917;  % kg/m3 
Tm       = -0.4; % ice melt T, C, depends on S

f_getT = 0;
if f_getT>0
 
  for ict=1:ncat;
    n1=(ict-1)*nilyr+1;
    n2=ict*nilyr;
    eicen = PT.eicen(n1:n2);
    vicen = PT.vicen(ict);
    qi = eicen*nilyr/vicen;

    a=c0;
    b=(cw-c0)*Tm-qi/rho_ice-L0;
    c=L0*Tm;

    Tlr = (-b-sqrt(b.^2-4*a*c))./(2*a);
    fprintf('Ice cat: %i, Surf T = %4.1f\n',ict, PT.trcrn(ict));
    fprintf('  T in layers: surf (1) - btm (4): \n');
    fprintf('  %4.1f, %4.1f, %4.1f, %4.1f\n',Tlr);

    esnon = PT.esnon(ict);
    vsnon = PT.vsnon(ict);
    qsn = esnon/vsnon; % 1 layer of snow
    Tsn = (L0 + qsn/rho_snow)/c0;
    fprintf('   T snow: %4.1f\n\n', Tsn);
  end
  keyboard
end

f_chck_esnon=1;
if f_chck_esnon>0
  for ict=1:ncat
    esnon=squeeze(ESNON(ict,:,:));
    vsnon=squeeze(VSNON(ict,:,:));
    Qsn = esnon./vsnon;
    Tsn = (L0 + Qsn./rho_snow)/c0;
    
    tmax=max(max(Tsn));
    tmin=min(min(Tsn));
    fprintf('   ### From Enthalpie Tsnow, cat, Tmin, Tmax:%i %9.4f %9.4f\n',...
	    ict,tmin,tmax);
    if tmax>1e-7
      fprintf('Tmax is too warm ...');
    end
  end
%  keyboard
end





% Velocity
fprintf('min/max velocity components\n');
cnt=cnt+1;
%fprintf('%i: Reading ice uvel\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % uvel
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if isfield(FSV,'ivel')
  flg=FSV.ivel;
  if flg>0
    U=AA;
  end
end


cnt=cnt+1;
%fprintf('%i: Reading ice vvel\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % vvel
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if isfield(FSV,'ivel')
  flg=FSV.ivel;
  if flg>0
    fmat = sprintf('%s%s_ivel.mat',PTH.data,fmat0);
    fprintf('Saving %s\n',fmat);
    V = AA;
    save(fmat,'U','V');
    clear V U
  end
end

%keyboard
splt=0;
if splt>0
  fgn=1;
  ttl='Ice vvel';
  sub_plot_cice(AA,fgn,ttl,HH);
end
 

clear uvel vvel

% Radiation fields:
% 4 radiative categories
% for calculating albedo for visible and IR wavelengths
% and penetrating sh/wave
% CICE asumes that all IR is absorbed at the surface
% only visible is absorbed in the ice interior or
% transmitted to the ocean
%
% Scale factor to change MKS units
% for shortwave components
% default = 1
fprintf('Radiation fields \n');
cnt=cnt+1;
%fprintf('%i: Reading scale factor\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % scale_factor
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.rad == 1
  scale_fct=AA;
end


% Shortwave down visible direct radiation, W/m2
cnt=cnt+1;
%fprintf('%i: Reading Sh/wv vis. dir\n',cnt);
dmm   = fread(fidc,1,'int');
AA    = (fread(fidc,[nn,mm],'float64'))'; % swvdr
dmm   = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.rad == 1
  swvdr=AA;
end

% Shortwave down visible diffusive radiation, W/m2
cnt=cnt+1;
%fprintf('%i: Reading Sh/wv vis. diff\n',cnt);
dmm   = fread(fidc,1,'int');
AA    = (fread(fidc,[nn,mm],'float64'))'; % swvdf
dmm   = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.rad == 1
  swvdf=AA;
end

% Shortwave down near IR direct radiation, W/m2
cnt=cnt+1;
%fprintf('%i: Reading Sh/wv IR dir\n',cnt);
dmm   = fread(fidc,1,'int');
AA    = (fread(fidc,[nn,mm],'float64'))'; % swidr
dmm   = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.rad == 1
  swidr=AA;
end

% Shortwave down near IR diffusive radiation, W/m2
cnt=cnt+1;
%fprintf('%i: Reading Sh/wv IR diff\n',cnt);
dmm   = fread(fidc,1,'int');
AA    = (fread(fidc,[nn,mm],'float64'))'; % swidf
dmm   = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.rad == 1
  swidf=AA;
  fmat = sprintf('%s%s_rad.mat',PTH.data,fmat0);
  fprintf('Saving %s\n',fmat);
  save(fmat,'scale_fct','swvdr','swvdf','swidr','swidf');
  clear scale_fct swvdr swvdf swidf
end

%clear swvd* swid*

splt=0;
if splt>0
  fgn=1;
  ttl='swvdf: Shortwave vis diffus, W/m2\n';
  sub_plot_cice(AA,fgn,ttl,HH);
end

% Ocean stress, N/m2
% xdirection
fprintf('Min/max ocean stress components \n');
cnt=cnt+1;
%fprintf('%i: Reading ocean str.X\n',cnt);
dmm = fread(fidc,1,'int');
AA  = (fread(fidc,[nn,mm],'float64'))'; % strocnxT
dmm = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.ostrs == 1
  strocnxT = AA;
end


% ydirection
cnt=cnt+1;
%fprintf('%i: Reading ocean str.Y\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % strocnyT
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.ostrs == 1
  strocnyT = AA;
  fmat = sprintf('%s%s_ostr.mat',PTH.data,fmat0);
  fprintf('Saving %s\n',fmat);
  save(fmat,'strocnxT','strocnyT');
  clear strocnxT strocnyT
end


%clear strocn*

% Internal stress, stress tensor, kg/s2
% (1) northeast, (2) northwest, (3) southwest, (4) southeast
fprintf('internal stress components \n');
cnt=cnt+1;
%fprintf('%i: Reading intern.stress p1\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % stressp_1
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stressp_1 = AA;
end

cnt=cnt+1;
%fprintf('%i: Reading intern.stress p3\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % stressp_3
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stressp_3 = AA;
end

cnt=cnt+1;
%fprintf('%i: Reading intern.stress p2\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % stressp_2
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stressp_2 = AA;
end

cnt=cnt+1;
%fprintf('%i: Reading intern.stress p4\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % stressp_4
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stressp_4 = AA;
  fmat = sprintf('%s%s_intstrp.mat',PTH.data,fmat0);
  fprintf('Saving %s\n',fmat);
  save(fmat,'stressp_1','stressp_2','stressp_3','stressp_4');
  clear stressp_*
end


%keyboard
splt=0;
if splt>0
  fgn=1;
  ttl='stressp 4: Internal stress, kg/s2';
  sub_plot_cice(AA,fgn,ttl,HH);
end


cnt=cnt+1;
%fprintf('%i: Reading intern.stress m1\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % stressm_1
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stressm_1 = AA;
end

cnt=cnt+1;
%fprintf('%i: Reading intern.stress m3\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % stressm_3
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stressm_3 = AA;
end


cnt=cnt+1;
%fprintf('%i: Reading intern.stress m2\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % stressm_2
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stressm_2 = AA;
end


cnt=cnt+1;
%fprintf('%i: Reading intern.stress m4\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % stressm_4
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stressm_4 = AA;
  fmat = sprintf('%s%s_intstrm.mat',PTH.data,fmat0);
  fprintf('Saving %s\n',fmat);
  save(fmat,'stressm_1','stressm_2','stressm_3','stressm_4');
  clear stressm_*
end



cnt=cnt+1;
%fprintf('%i: Reading intern.stress 12_1\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % stress 12_1
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stress12_1 = AA;
end


cnt=cnt+1;
%fprintf('%i: Reading intern.stress 12_3\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % 12_3
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stress12_3 = AA;
end


cnt=cnt+1;
%fprintf('%i: Reading intern.stress 12_2\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % 12_2
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stress12_2 = AA;
end


cnt=cnt+1;
%fprintf('%i: Reading intern.stress 12_4\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % 12_4
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.int_strs == 1
  stress12_4 = AA;
  fmat = sprintf('%s%s_intstr12.mat',PTH.data,fmat0);
  fprintf('Saving %s\n',fmat);
  save(fmat,'stress12_1','stress12_2','stress12_3','stress12_4');
  clear stress12_*
end



% Ice mask
fprintf('ice mask for dynamics \n');
cnt=cnt+1;
%fprintf('%i: Reading intern.stress ice mask\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % iceumask
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

if FSV.imask == 1
  iceumask = AA;
  fmat = sprintf('%s%s_iceumask.mat',PTH.data,fmat0);
  fprintf('Saving %s\n',fmat);
  save(fmat,'iceumask');
  clear iceumask
end


%keyboard
splt=0;
if splt>0
  fgn=1;
  ttl='Ice Mask';
  sub_plot_cice(AA,fgn,ttl,HH);
  keyboard
end


% For trully coupled HYCOM-CICE these fields
% are not needed
% if defined ocean mixed layer in CICE
fprintf('Ocean mixed layer min/max\n');
cnt=cnt+1;
%fprintf('%i: Reading ocean SST\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))';
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);

splt=0;
if splt>0
  fgn=1;
  ttl='Ocean ML: SST';
  sub_plot_cice(AA,fgn,ttl,HH);
end


if isempty(dmm); 
  fprintf('Reached EOF: Ocean ML data are not in the restart file\n');
  cnt=cnt-1;
end


% ocean mixed L: freezing melting potential, W/m2
fprintf('Ocean frzmlt\n');
cnt=cnt+1;
%fprintf('%i: Reading ocean frzmlt\n',cnt);
dmm  = fread(fidc,1,'int');
AA   = (fread(fidc,[nn,mm],'float64'))'; % frzmlt
dmm  = fread(fidc,1,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: %16.6f %16.6f\n',cnt,mina,maxa);
if isempty(dmm); 
  fprintf('Reached EOF\n');
  cnt=cnt-1;
end

splt=0;
if splt>0
  fgn=1;
  ttl='Ocean ML: frzmlt';
  sub_plot_cice(AA,fgn,ttl,HH);
end

% This should be the end
%keyboard


fclose(fidc);

fprintf('Read in %i fields\n',cnt);

f_esn=0;
if f_esn, % calculate T snow from enthalpie
  aicen=PT.aicen;
  esnon=PT.esnon; % snow enthalpie
  vsnon=PT.vsnon; % snow volume/m2
  qsn = esnon./vsnon;
  
  T = (L0+qsn/rho_snow)/c0;
  
end
