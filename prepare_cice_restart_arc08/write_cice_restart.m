% Prepare cice restart for HYCOM-CICE 0.08
% to match sea ice distribution in the Arctic
% on the day of the model simulation start
% using existing old restart file
% and PIOMAS sea ice fields
% http://psc.apl.uw.edu/research/projects/...
%      arctic-sea-ice-volume-anomaly/data/model_grid
% citation: Zhang, Jinlun and D.A. Rothrock: 
% Modeling global sea ice with a thickness and 
% enthalpy distribution model in generalized 
% curvilinear coordinates, Mon. Wea. Rev. 131(5), 681-697, 2003.
%
% Temporary fields for restart files were created
% in remap_piomas2arc.m & read_cice_restart.m
%

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

% Dates of the restart file;
YY = 1993;
MM = 1; 
DD = 1;

fmat0 = 'cice_restart093j'; % cice template, any year/mo


PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.data = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/PIOMAS_ice_data/';
PTH.rest = '/nexsan/people/ddmitry/hycom/ARCc0.08/103/ice_restart/';
PTH.mat  = PTH.rest;
fltopo=sprintf('%sdepth_ARCc0.08_11.nc',PTH.topo);
HH   = nc_varget(fltopo,'Bathymetry');
%alat = nc_varget(fltopo,'Latitude');
%elon = nc_varget(fltopo,'Longitude');
%LAT  = alat;
%LON  = elon;
%[m,n]= size(HH);
[mm,nn]= size(HH);

fcice = sprintf('%scice.restart_piomas_ARCc0.08_T11_093a',PTH.rest);
fid = fopen(fcice,'w','ieee-be');

if fid>0 
  fprintf('Open file %s\n',fcice);
else
  fprintf('ERR: cannot open file %s\n',fcice);
end

% Domain size, categories, layers:
% see ice_dominan_size.F90
ncat  = 5; % # of categories
nilyr = 4; % ice layers per category
ntilyr = ncat*nilyr; % # of ice layers in all categories
nslyr = 1; % # of snow layers
ntslyr = ncat*nslyr; % # of snow layers in all cat.
ntrcr = 3; % # of tracers, 1- surf. T


c0 = 2106;    % J/(kg*Cdeg) - specific heat of fresh ice, at 0C
L0 = 3.34e5; % J/kg  latent heat of fusion of fresh ice at 0C
mu = 0.054;  % deg/ppt - liquidus, ratio betw. Tfreez and S of brine
cw = 4218;   % J/(kg*Cdeg) - specific heat of sea water
rho_snow = 330;  % kg/m3
rho_ice  = 917;  % kg/m3 
Tm       = -0.4; % ice melt T, C, depends on S



% Ice mask
A = sub_restart_fields('iceumask',[]);
IMSK = A.IMSK;
clear A

% Write header:
[d1,d2,d3] = get_dtime(datenum(YY,MM,DD));
istep=(datenum(YY,MM,DD)-datenum(1900,12,31))*24*3600/450; % cice t.step since Jan.1,1901
rday = d1*86400; % day in HYCOM day convention, secs
fday=1*86400;   % seconds
fprintf('Writing header: istep= %i, time= %i, time_forc= %i\n',...
	istep,rday/86400,fday/86400);
nbt=20;
fwrite(fid,nbt,'int');
fwrite(fid,istep,'int');
fwrite(fid,rday,'float64');
fwrite(fid,fday,'float64');
fwrite(fid,nbt,'int');
%keyboard

% ------------------------
% State variables
% ------------------------
% 
%     area, vol ice, vol snow T surf by cat & layers
%fmat = sprintf('%srest_ice_state_%4.4i%2.2i.mat',...
%		 PTH.rest,YY,MM);
%fprintf('Loading %s\n',fmat);
%load(fmat);
%nin = length(ICE); % # of sea ice points

IN = find(IMSK==1);
nin = length(IN);
%if length(IN) ~= nin
%  fprintf('Ice mask IMSK does not match ice points in ICE\n');
%  fprintf(' Check IMSK, IMSK can be created from ICE(i).Index');
%  error('  IMSK /= Ice points in ICE');
%end
%if IMSK(1) ~= ICE(1).Indx
%  fprintf('Ice mask IMSK does not match ice points in ICE\n');
%  fprintf(' Check IMSK, IMSK can be created from ICE(i).Index');
%  error('  IMSK(1) differe from ICE(1).Indx');
%end
% Create IMSK:
%for ii =1:nin
%  IMSK(ii,1) = ICE(ii).Indx;
%end

fmat = sprintf('%srest_ice_aicen_%4.4i%2.2i.mat',...
		PTH.mat,YY,MM);
fprintf('Loading %s\n',fmat);
load(fmat);

fmat = sprintf('%srest_ice_vicen_%4.4i%2.2i.mat',...
		PTH.mat,YY,MM);
fprintf('Loading %s\n',fmat);
load(fmat);

fmat = sprintf('%srest_ice_vsnon_%4.4i%2.2i.mat',...
		PTH.mat,YY,MM);
fprintf('Loading %s\n',fmat);
load(fmat);

fmat = sprintf('%srest_ice_trcrn_%4.4i%2.2i.mat',...
		PTH.mat,YY,MM);
fprintf('Loading %s\n',fmat);
load(fmat);

% Check sea ice mask if it matches saved fields:
nif = length(AICEN);

if nif ~= nin
  error('Check IMASK, does not match saved ice fields');
end


cnt=0;
nbt = nn*mm*8;
for ic = 1:ncat
  fprintf('Category: %i min/max area, vol ice, vol snow, Tsfc\n',ic);
  cnt=cnt+1;
  
  AA = zeros(mm,nn);
  AA(IN) = AICEN(:,ic);
  AICENc(ic,:,:)=AA;
  AA=AA';
  fwrite(fid,nbt,'int');
  fwrite(fid,AA,'float64');
  fwrite(fid,nbt,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: cat %i, %16.6f %16.5f\n',cnt,ic,mina,maxa);

  cnt=cnt+1;
  AA = zeros(mm,nn);
  AA(IN) = VICEN(:,ic);
  AA=AA';
  fwrite(fid,nbt,'int');
  fwrite(fid,AA,'float64');
  fwrite(fid,nbt,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: cat %i, %16.6f %16.5f\n',cnt,ic,mina,maxa);
  
  cnt=cnt+1;
  AA = zeros(mm,nn);
  AA(IN) = VSNON(:,ic);
  VSNONc(ic,:,:) = AA;
  AA=AA';
  fwrite(fid,nbt,'int');
  fwrite(fid,AA,'float64');
  fwrite(fid,nbt,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: cat %i, %16.6f %16.5f\n',cnt,ic,mina,maxa);
  
  cnt=cnt+1;
  AA = zeros(mm,nn);
  AA(IN) = TSFC(:,ic);
  AA=AA';
  fwrite(fid,nbt,'int');
  fwrite(fid,AA,'float64');
  fwrite(fid,nbt,'int');
  mina=min(min(AA));
  maxa=max(max(AA));
  fprintf(' %i: cat %i, %16.6f %16.5f\n',cnt,ic,mina,maxa);
  
end  % ncat - write state variables by categories

clear AICEN VICEN VSNON TSFC
%clear AICEN VICEN TSFC

% ---------------
% Ice Enthalpie
% ---------------
fmat = sprintf('%srest_ice_eicen_%4.4i%2.2i.mat',...
		PTH.mat,YY,MM);
fprintf('Loading %s\n',fmat);
load(fmat);
fprintf('min/max eicen for each layer:\n');
for ic=1:ncat
  for ilr=1:nilyr
    cnt=cnt+1;
    AA = zeros(mm,nn);
    AA(IN) = EICEN(:,ic,ilr); %  
%    AA(IN) = EICEN(:,ic,ilr)/nilyr; %
% Check for errors: 
    ierr = find(EICEN(:,ic,ilr)>0);
    if ~isempty(ierr)
      fprintf('Found >0 enthalpie: ic=%i, ilr=%i\n',ic,ilr);
      keyboard
%      ee = squeeze(EICEN(:,ic,ilr));
%      ee = sub_fix_enthalpie(ee,ierr,ic,ilr,PTH,YY,MM);
%      AA(IN) = ee; 
%      keyboard;
    end
    
    AA=AA';
    fwrite(fid,nbt,'int');
    fwrite(fid,AA,'float64');
    fwrite(fid,nbt,'int');
    mina=min(min(AA));
    maxa=max(max(AA));
    fprintf(' %i: cat=%i layer=%i, %16.6f %16.5f\n',cnt,ic,ilr,mina,maxa);
  end  
end 
clear EICEN

% ----------------
% Snow Enthalpie
% ----------------
fmat = sprintf('%srest_ice_esnon_%4.4i%2.2i.mat',...
		PTH.mat,YY,MM);
fprintf('Loading %s\n',fmat);
load(fmat);
fprintf('min/max esnon for each layer:\n');
for ic=1:ncat
    ict=ic;
    cnt=cnt+1;
    AA = zeros(mm,nn);
%    IE=find(abs(ESNON<1e-7));
%    ESNON(IE)=-1e-9;
    AA(IN) = ESNON(:,ic);  
    esnon=AA;
    AA=AA';
    fwrite(fid,nbt,'int');
    fwrite(fid,AA,'float64');
    fwrite(fid,nbt,'int');
    mina=min(min(AA));
    maxa=max(max(AA));
    fprintf(' %i: cat=%i, %16.6f %16.5f\n',cnt,ic,mina,maxa);

    aicen=squeeze(AICENc(ict,:,:));
    vsnon=squeeze(VSNONc(ict,:,:));
    hsn=vsnon./aicen;
    II=find(hsn<1e-7);
    
    Qsn = esnon./vsnon;
%    IE = find(abs(esnon)<=1e-7);
%    Qsn(IE) = nan;
    Tsn = (L0 + Qsn./rho_snow)/c0;
    Tsn(II) = nan;
    
    tmax=max(max(Tsn));
    tmin=min(min(Tsn));
    fprintf('### Check Tsnow, cat, Tmin, Tmax:%i %9.4f %9.4f\n',...
	    ict,tmin,tmax);
    if tmax>1e-7
      fprintf('Tmax is too warm ...\n');
      keyboard
    end
    if tmin<-100  % min T in CICE
      fprintf('Tmin is too low ...\n');
      keyboard
    end
    
end 
clear ESNON

k=0;
% ------------------
% Ice velocity
% ------------------
A = sub_restart_fields('ivel',IMSK);
fprintf('min/max velocity components\n');
cnt=cnt+1;
AA=A.U';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i:    %16.6f %16.5f\n',cnt,mina,maxa);

cnt=cnt+1;
AA=A.V';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i:    %16.6f %16.5f\n',cnt,mina,maxa);

clear AA
% ------------------
% Radiation fields
% ------------------
A = sub_restart_fields('rad',IMSK);
fprintf('min/max Radiation Fields\n');

cnt=cnt+1;
AA=A.scale_factor';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i:    %16.6f %16.5f\n',cnt,mina,maxa);

cnt=cnt+1;
AA=A.swvdr';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i:   %16.6f %16.5f\n',cnt,mina,maxa);

cnt=cnt+1;
AA=A.swvdf';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i:   %16.6f %16.5f\n',cnt,mina,maxa);

cnt=cnt+1;
AA=A.swidr';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i:   %16.6f %16.5f\n',cnt,mina,maxa);

cnt=cnt+1;
AA=A.swidf';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i:   %16.6f %16.5f\n',cnt,mina,maxa);

clear A 

% -----------------
% ocean stress
% -----------------
A = sub_restart_fields('ostress',IMSK);
fprintf('min/max ocean stress components\n');

cnt=cnt+1;
AA=A.strocnxT';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i:   %16.6f %16.5f\n',cnt,mina,maxa);

cnt=cnt+1;
AA=A.strocnyT';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i:   %16.6f %16.5f\n',cnt,mina,maxa);

k=0;
% -----------------
% Internal Stress
% -----------------
A = sub_restart_fields('instress',IMSK);
fprintf('min/max internal stresses\n');

cnt=cnt+1;
AA=A.stressp_1';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stressp_3';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stressp_2';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stressp_4';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stressm_1';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stressm_3';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stressm_2;
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stressm_4';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stress12_1';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stress12_3';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stress12_2';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);

cnt=cnt+1;
AA=A.stress12_4';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: cat/layer %i, %16.6f %16.5f\n',cnt,k,mina,maxa);


clear AA

% ----------------
% Sea ice umask 
% ----------------
%AA = zeros(mm,nn);
%AA(IN) = 1;
AA = IMSK;
cnt=cnt+1;
AA=AA';
fwrite(fid,nbt,'int');
fwrite(fid,AA,'float64');
fwrite(fid,nbt,'int');
mina=min(min(AA));
maxa=max(max(AA));
fprintf(' %i: Ice Mask %16.6f %16.5f\n',cnt,mina,maxa);

fprintf('Written all fields in %s\n',fcice);
fclose(fid);




