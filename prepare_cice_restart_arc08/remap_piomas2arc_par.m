% Parallel code
%
% Remap PIOMAS to ARCc0.08 grid
% Find interpolation weights
% as inverse distance weights for 
% HYCOM grid points
% 
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

s = matlabpool('size');
if s==0
  matlabpool open 12;
end



s_mat_indx=0; % =1 - derive indices/weights for interpolation, save mat
	      % =0 - load saved mat with weights
s_mat = 0;  % =1 - save PIOMAS fields interpolated into HYCOM 
	 
if s_mat == 0
  fprintf('No fields will be saved ...\n\n');
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
  txtb='hycom_NOPP_rivers/prepare_cice_restart/remap_piomas2arc.m';
  bottom_text(txtb);
end


% CICE thickness categories
% May find in the output *log file:
% hin_max(n-1) < Cat n < hin_max(n)
%  0.000000000000000E+000  < Cat            1  <   0.644507216819426
%  0.644507216819426       < Cat            2  <    1.39143349757630
%   1.39143349757630       < Cat            3  <    2.47017938195989
%   2.47017938195989       < Cat            4  <    4.56728791885049
%   4.56728791885049       < Cat            5  <    9.33384181586817

hi_min=0.1;
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

AICEN = zeros(nin,ncat);
VICEN = zeros(nin,ncat);

parfor ii=1:nin
  dmm = mean(AICEN,2);
  II = find(dmm>0);
  npsd = length(II);
    if mod(npsd,10000)==0
      fprintf('Ice points, ii=%i %6.2f%%\n',npsd,npsd/nin*100);
    end

    i0=IN(ii);
    hi=FF(i0);
    sg=0.5*hi;
  % Normally distribute with mean hi
    Hc = sg*(randn(100,1))+hi;
    nH=length(Hc);
    Hc(Hc<hi_min)=0;
    Nc0=find(Hc==0);
    iW=length(Nc0)/nH;  % fraction of open water
    area0 = Acell(i0);

  % aicen: Ice area fraction, by categories:
    for ic=1:ncat
      h1=icat(ic);
      h2=icat(ic+1);
      Jc=find(Hc>=h1 & Hc<h2);
      if isempty(Jc); 
        AICEN(ii,ic) = 0;
	VICEN(ii,ic) = 0;
	continue; 
      end;
      aicen=length(Jc)/nH; % fraction ice in this cat. 
%      ICE(ii).aicen(ic)=aicen;
      AICEN(ii,ic) = aicen;

  % Ice volume:
  %volume per unit area of ice (m)
  % vicen(i,j)=aicen(i,j)*hin(i,j)
      hin = mean(Hc(Jc));
%      ICE(ii).vicen(ic)=aicen*hin;
      VICEN(ii,ic) = aicen*hin;
      
  % Snow volume:
  % volume per unit area of snow (m)
  % vsnon(i,j)=aicen(i,j)*hsn(i,j) 
      hsn=SN(i0);
%      ICE(ii).vsnon(ic)=aicen*hsn;
      VSNON(ii,ic) = aicen*hsn;
      
  % Ice tracer: T of ice/snow top surf.
%      ICE(ii).Tsfc(ic) = TI(i0);
      TSFC(ii,ic) = TI(i0);
      
  % eicen: % Energy of melting for each ice layer (J/m2)
  % (it equals to the negative enthalpie)
  % This is the energy required to melt a unit volume of 
  % ice and raise its T to 0C
  % Calculate enthalpy for the whole category assuming
  % T cat = half of T surf and T freeze
  % and distribute enthalpie equally between the layers
      Tocn = -1.8;
  %    Tice = 0.5*(TI(i0)-1.8);
      Tsrf = min([TI(i0),-0.1]);
      dT = Tsrf-Tocn;
      dTl= dT/nlr; % T change in 1 layer
      for ilr = 1:nlr
	Tlr_btm = (ilr-1)*dTl+Tocn;
	Tlr_top = ilr*dTl+Tocn;
	Tlr = 0.5*(Tlr_btm+Tlr_top);
	Si = 3;
	qi = sub_enthalpie('ice',Si,Tlr); % J/m3
%	ICE(ii).eicn(ilr,ic) = qi*aicen*hin; % J/m2
	EICEN(ii,ic,ilr) = qi*aicen*hin; % J/m2
      end

  % esnon - snow enthalpy (J/m2) for each layer    
  %         1 snow layer x ncat ice categories
      qs = sub_enthalpie('snow',0,Tsrf);
%      ICE(ii).esnon(ic) = qs*aicen*hsn; % J/m2
      ESNON(ii,ic) = qs*aicen*hsn/ncat; % J/m2, per 1 ice cat.

    end;  % ice categories
end;   % parfor all sea ice points

  
s = matlabpool('size');
if s==0
  matlabpool close;
end
  
  % Save fields
  if s_mat > 0;
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


  



