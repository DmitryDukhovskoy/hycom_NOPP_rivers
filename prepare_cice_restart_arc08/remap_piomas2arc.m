% Remap PIOMAS to ARCc0.08 grid
% Find interpolation weights
% as inverse distance weights for 
% HYCOM grid points
% Note: CICE hin(n) = vicen(n)/aicen(n) - thickness category=n = ice vol.(n)/ ice area (n)
%  vicen is ice volume per 1 m2 of area, cat. n
%  or this is the grid-cell area mean thickness 
%
% PIOMAS gives (area) mean ice thickness for the grid cell!
%
% COAPS, FSU, D. Dukhovskoy, Aug. 2016,
% Fixed bug with sea ice thickness, Jan 2018
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

s_mat_indx = 0; % =1 - derive indices/weights for interpolation, save mat
	      % =0 - load saved mat with weights
s_mat   = 0;  % =1 - save PIOMAS fields interpolated into HYCOM 
s_njobs = 1;  % divide all points into several jobs and run indendpently
              % =0 or =1 - run all points in 1 time
s_jobid = 1;  % current job id, from 1, ..., to s_njobs


if s_mat == 0
  fprintf('\n ====   No fields will be saved ... ===\n\n');
end


% Fields to remap/save:	 
%FLD.state_var = 1; % state var: thickn, area, enthalpie by laeyrs/categories
%FLD.velocity  = 1; % ice vel. comp.
%FLD.rad       = 1; % radiation
%FLD.ostress   = 1; % ocean stresses
%FLD.istress   = 1; % internal stresses
%FLD.imask     = 1; % ice mask

% Restart time stamp:
YR = 1993;
MM = 1;
fmat0 = 'cice_restart093j'; % template for saved rest. fields 
                            % from existing CICE restart
	 
PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.data = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/PIOMAS_ice_data/';
PTH.rest = '/nexsan/people/ddmitry/hycom/ARCc0.08/103/ice_restart/';
PTH.mat  = PTH.rest;

%fmat=sprintf('%sremap_piomas2hycom.mat',PTH.rest);

fgrds = sprintf('%sgrid.dat',PTH.data); % grid for scalar fields
nxl=360;
nyl=120;
nxy=nxl*nyl;

% Create template of sea ice distribution
% with mean 0
%  Hc = sg*(randn(100,1))+hi;
%xx=[-3.5:0.07:3.49];
%fprb = 1/sqrt(2*pi)*exp(-xx.^2/2);
%rng('default');
%hi_tmp = randn(100,1);
%save('hi_tmp.mat','hi_tmp');
load('hi_tmp');
% Adjust - to make mean=0;
mn=mean(hi_tmp);
hi_tmp=hi_tmp-mn;



% read lon/lat scalar fields
dmm = load(fgrds);
[a1,a2]=size(dmm);
nrw=nxy/a2;

LONp = dmm(1:nrw,:);
LONp = reshape(LONp',[nxy,1]);
LONp = reshape(LONp,[nxl,nyl])';
[mp,np]=size(LONp);

LATp = dmm(nrw+1:end,:);
LATp = reshape(LATp',[nxy,1]);
LATp = reshape(LATp,[nxl,nyl])';

% HYCOM:
fltopo=sprintf('%sdepth_ARCc0.08_11.nc',PTH.topo);
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LATh  = alat;
LONh  = elon;
[m,n]= size(HH);
[mm,nn]= size(HH);

% Grid cell spacing
[DX,DY]=sub_dx_dy(elon,alat);
Acell=DX.*DY; % Grid cell area, m2


H2P = sub_intrp_indx_piomas2hycom(s_mat_indx);

% Read PIOMAS:
% For scalar data use LONp/LATp
fld='ithkn';
yr=1993;
imo=1; 
APIOMAS = sub_read_piomas(fld,yr,imo);
FP = APIOMAS.Field;

% Interpolate PIOMAS to HYCOM
FF = HH*0;
WW = H2P.Weights;
IP = H2P.PIOMAS_LinIndx;
IH = H2P.HYCOM_IceLinIndx;
dmm = FP(IP(:,1)).*WW(:,1)+...
      FP(IP(:,2)).*WW(:,2)+...
      FP(IP(:,3)).*WW(:,3)+...
      FP(IP(:,4)).*WW(:,4)+...
      FP(IP(:,5)).*WW(:,5);
FF(IH)=dmm;      
FF(FF==0)=nan;

% PIOMAS Sea ice T:
tK = 273.15;
fld='itemp';
APIOMAS = sub_read_piomas(fld,yr,imo);
FP = APIOMAS.Field;
FP(FP==0) = tK-10; % land
FP=FP-tK;  % K -> C

% Interpolate PIOMAS to HYCOM
TI = HH*0;
WW = H2P.Weights;
IP = H2P.PIOMAS_LinIndx;
IH = H2P.HYCOM_IceLinIndx;
dmm = FP(IP(:,1)).*WW(:,1)+...
      FP(IP(:,2)).*WW(:,2)+...
      FP(IP(:,3)).*WW(:,3)+...
      FP(IP(:,4)).*WW(:,4)+...
      FP(IP(:,5)).*WW(:,5);
TI(IH)=dmm;      
TI(TI==0)=nan;

%
% PIOMAS Snow
fld='snow';
APIOMAS = sub_read_piomas(fld,yr,imo);
FP = APIOMAS.Field;

% Interpolate PIOMAS to HYCOM
SN = HH*0;
WW = H2P.Weights;
IP = H2P.PIOMAS_LinIndx;
IH = H2P.HYCOM_IceLinIndx;
dmm = FP(IP(:,1)).*WW(:,1)+...
      FP(IP(:,2)).*WW(:,2)+...
      FP(IP(:,3)).*WW(:,3)+...
      FP(IP(:,4)).*WW(:,4)+...
      FP(IP(:,5)).*WW(:,5);
SN(IH)=dmm;      
SN(SN==0)=nan;



fplt=0;
if fplt>0
  figure(1); clf;
  pcolor(FF); shading flat;
  hold on;
  contour(FF,[0.2 0.2],'m');
  contour(HH,[0 0],'k');
  caxis([0 3]);
  colorbar
  axis('equal');
  set(gca,'xlim',[1 nn],...
	  'ylim',[1 mm]);
  stt=sprintf('PIOMAS %i/%i intrp to HYCOM',imo,yr);
  title(stt);
  txtb='remap_piomas2arc.m';
  bottom_text(txtb,'pwd',1);
  keyboard
end


% CICE thickness categories
% May find in the output *log file:
% hin_max(n-1) < Cat n < hin_max(n)
%  0.000000000000000E+000  < Cat            1  <   0.644507216819426
%  0.644507216819426       < Cat            2  <    1.39143349757630
%   1.39143349757630       < Cat            3  <    2.47017938195989
%   2.47017938195989       < Cat            4  <    4.56728791885049
%   4.56728791885049       < Cat            5  <    9.33384181586817

hi_min=0.05;
icat=[hi_min;0.644507216819426; 1.39143349757630;...
      2.47017938195989; 4.56728791885049; ...
      100];

ncat = length(icat)-1; % 5 categories of ice
nlr  = 4;              % 4 layers in each cat.

% Create Sea ICe Mask based on HH and PIOMAS ice conc.
FF(FF==nan)=0;
IMSK = FF*0;
IMSK(isnan(IMSK))=0;
I=find(FF>0.1 & HH<0);
IMSK(I)=1;
if s_mat > 0;
  fmat = sprintf('%shycom_ice_mask_%4.4i%2.2i.mat',...
		 PTH.mat,YR,MM);
  fprintf('Saving %s\n',fmat);
  save(fmat,'IMSK');
end

% Save Sea Ice thicknesses:
% Assume Gaussian distribution of thickness
% with mean - from PIOMAS
xx=[-5:.1:5];
sg=2.;
%Hw = exp(-(xx).^2/(2*sg.^2));
%Hw = length(xx)*Hw./sum(Hw);


%B=D;
IN = find(IMSK==1);
nin=length(IN);
%if FLD.state_var == 1
fprintf('Remapping state variables, go over ice grid cells ...\n');

iiS = 1;
iiE = nin;

if s_njobs>1
  dNp = round(nin/s_njobs);
  iiS = (s_jobid-1)*dNp+1;
  iiE = iiS+dNp;
  iiE = min([iiE,nin]);

  fprintf('Partitioned job, n_jobs %i\n',s_njobs);
  fprintf('Job # %i, indices: %i - %i\n\n',s_jobid,iiS,iiE);
  
end

%keyboard
iir = 0;
for ii=iiS:iiE;
%for ii=522216:522216;
  if mod(ii,20000)==0
    fprintf('Ice points, ii=%i %6.2f%%\n',ii,(ii-iiS+1)/(iiE-iiS+1)*100);
  end

  iir = iir+1; % in normal case, iir=ii
  i0=IN(ii);
  INDX(iir,1) = i0; 
%keyboard  
  hi=FF(i0);
  sg=0.5*hi;
% Normally distribute with mean hi
%  Hc = sg*(randn(100,1))+hi;
  Hc = sg*hi_tmp+hi;
  nH=length(Hc);
  Hc(Hc<hi_min)=0;
% adjust mean H ice  
  dh = hi-mean(Hc);
  icc=0;
  while abs(dh)>1e-5
    Hc = Hc+dh;
    Hc(Hc<hi_min)=0;
    dh = hi-mean(Hc);
    icc=icc+1;
    if icc>1000, error('Endless loop'); end;
  end
  
  Nc0=find(Hc==0);
  iW=length(Nc0)/nH;  % fraction of open water
  area0 = Acell(i0);
%keyboard
% aicen: Ice area fraction, by categories:
  for ic=1:ncat
    h1=icat(ic);
    h2=icat(ic+1);
    Jc=find(Hc>=h1 & Hc<h2);
    if isempty(Jc); 
      AICEN(iir,ic) = 0;
      VICEN(iir,ic) = 0;
      continue; 
    end;
    aicen=length(Jc)/nH; % fraction ice in this cat. 
    AICEN(iir,ic) = aicen;
   

% Ice volume:
% volume per unit area of ice (m) per category
    hin = mean(Hc(Jc)); % ice thckn, category=ic
%    VICEN(iir,ic) = aicen*hin/ncat;
    VICEN(iir,ic) = aicen*hin; % already hin for category=ic

% Snow volume:
% volume per unit area of snow (m)
% vsnon(i,j)=aicen(i,j)*hsn(i,j)
% Need to divide by # of ice categories
    hsn=SN(i0);
%      ICE(ii).vsnon(ic)=aicen*hsn;
    VSNON(iir,ic) = aicen*hsn/ncat;

% Ice tracer: T of ice/snow top surf.
%      ICE(ii).Tsfc(ic) = TI(i0);
    TSFC(iir,ic) = TI(i0);

% eicen: % Energy of melting for each ice layer (J/m2)
% (it equals to the negative enthalpie)
% This is the energy required to melt a unit volume of 
% ice and raise its T to 0C
% Calculate enthalpy for the whole category assuming
% T cat = half of T surf and T freeze
% and distribute enthalpie equally between the layers
% Note: sea ice layers go downward
% Layer 1 - ice/atm, Layer 4 - ice/ocean
    Tocn = -3.; % ice near ocean interf - from CICE ~-3.3
    Tm   = -0.2; % ice melt T
%    Tice = 0.5*(TI(i0)-1.8);
    T0 = TI(i0);
    if T0<0;
      Tsrf = 0.7*T0;
    else
      Tsrf = Tm;
    end
% Note ice layer numbering increases downward    
% in sea ice, bottom (ocean) layer ~ -3.2C
% top is much warmer than Tsurf (in winter)
%    Tsrf = min([Tsrf,Tm]);
    Tsrf = min([Tsrf,Tm]);
    Tsrf_ice = min([0.4*Tsrf,Tm]);
    dTl = (Tocn-Tsrf_ice)/nlr; % T change in 1 layer
    for ilr = 1:nlr
      Tlr_top = (ilr-1)*dTl+Tsrf_ice;
      Tlr_btm = ilr*dTl+Tsrf_ice;
      Tlr = 0.5*(Tlr_btm+Tlr_top);
      Si = 3;
      qi = sub_enthalpie('ice',Tlr,Tm); % J/m3
%	ICE(ii).eicn(ilr,ic) = qi*aicen*hin; % J/m2
%        dmm = qi*aicen*hin/nlr
      vicen = VICEN(iir,ic)/nlr; 
      eicen = qi*vicen;
      EICEN(iir,ic,ilr) = eicen; % J/m2
    end

% esnon - snow enthalpy (J/m2) for each layer    
%         1 snow layer x ncat ice categories
    qs = sub_enthalpie('snow',Tsrf,0); % J/m3
%    vsnon = aicen*hsn/ncat;  % Need to divide by # of ice cat.
    vsnon = VSNON(iir,ic); % check if divided by ncat before
    esnon = qs*vsnon;
    ESNON(iir,ic) = esnon; % J/m2, per 1 ice cat.
% Reconstruct T snow following CICE:    
    [Tsn,Tmax] = sub_tsnow_enthalpie(esnon,vsnon,aicen); % 

    if Tsn>0 
      fprintf('Tsnow > 0: %8.2f \n',Tsn);
      keyboard;
    end
    if Tsn<-100 
      fprintf('Tsnow is too low<-100: %8.2f \n',Tsn);
      keyboard;
    end
%keyboard    
  end;  % ice categories
  
% Check:
% Grid cell area mean ice thickness should be = PIOMAS hi
  hmm = 0;
  for ic=1:ncat
    hin = VICEN(iir,ic); % grid-cell area mean ice thkns = VICEN (m3/1 m2 of area = m, mean h)
    hmm = hmm+hin;
  end;
  
%keyboard
end;   % all sea ice points

% Save fields
if s_mat > 0 & s_njobs<=1;
%    fmat = sprintf('%srest_ice_state_%4.4i%2.2i.mat',...
%		   PTH.mat,YR,MM);
%    fprintf('Saving %s\n\n',fmat);
%    save(fmat,'ICE');
%    clear ICE
  fmat = sprintf('%srest_ice_aicen_%4.4i%2.2i.mat',...
		 PTH.mat,YR,MM);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'AICEN');

  fmat = sprintf('%srest_ice_vicen_%4.4i%2.2i.mat',...
		 PTH.mat,YR,MM);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'VICEN');

  fmat = sprintf('%srest_ice_vsnon_%4.4i%2.2i.mat',...
		 PTH.mat,YR,MM);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'VSNON');

  fmat = sprintf('%srest_ice_trcrn_%4.4i%2.2i.mat',...
		 PTH.mat,YR,MM);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'TSFC');

  fmat = sprintf('%srest_ice_eicen_%4.4i%2.2i.mat',...
		 PTH.mat,YR,MM);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'EICEN');

  fmat = sprintf('%srest_ice_esnon_%4.4i%2.2i.mat',...
		 PTH.mat,YR,MM);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'ESNON');

end
%end;   % state variables

% for parallel jobs - save all fields
% for each job id but TSFC - ice T
if s_mat > 0 & s_njobs>1
  fmat = sprintf('%sTMP_aicen_%4.4i%2.2i_%2.2i.mat',...
		 PTH.mat,YR,MM,s_jobid);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'AICEN','iiS','iiE','INDX');

  fmat = sprintf('%sTMP_vicen_%4.4i%2.2i_%2.2i.mat',...
		 PTH.mat,YR,MM,s_jobid);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'VICEN','iiS','iiE','INDX');

  fmat = sprintf('%sTMP_vsnon_%4.4i%2.2i_%2.2i.mat',...
		 PTH.mat,YR,MM,s_jobid);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'VSNON','iiS','iiE','INDX');

  fmat = sprintf('%sTMP_eicen_%4.4i%2.2i_%2.2i.mat',...
		 PTH.mat,YR,MM,s_jobid);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'EICEN','iiS','iiE','INDX');

  fmat = sprintf('%sTMP_esnon_%4.4i%2.2i_%2.2i.mat',...
		 PTH.mat,YR,MM,s_jobid);
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'ESNON','iiS','iiE','INDX');
end

exit





  



