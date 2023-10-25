% Save depth-averaged fluxes 
% Extracted T, S, normal U at the specified straits
% for processing in python
%  Corrected flux calculation with T,S,U allocation is used
% see anls_fluxes/extr_TSVdaily_straits04.m
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear

f_save = 1;
snm = 'FramStr';

expt=023;
TV=17;
YR1=2017;
YR2=2019;
dday=7;


pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.04/data_straits/';
%pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.04/data_theresa/';
pthmat = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/023/strait_fluxes/';

Cp = 4200; % J/kg K
Tref1 = -1.8; % Ref T to calc. H flux
Tref2 = 0;    % Ref T to calc. H flux
Sref1 = 34.8;
Sref2 = 34.9;
hgg=1e20;

ixx    = 9; % experiment name and dir - check with EXPT - expt 023
%ixx    = 6;  % expt 022 original 
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;



% Process data:
TM   = [];
Vflx = [];
FWflx1 = [];
FWflx2 = [];
Hflx1  = [];
Hflx2  = [];
for YR = YR1:YR2
	yr=YR;
	dE=datenum(yr,12,31);
	dJ1=datenum(yr,1,1);
	ndays=dE-dJ1+1;

%    tmy = [dJ1:dday:dE]';

%	fmatout=sprintf('%shycom004_%3.3i_StraitFluxesDay_%4.4i.mat',...
%									pthmat,expt,YR);
  fmatout=sprintf('%shycom004_%3.3i_%s_StraitFluxes_%4.4i.mat',...
                    pthmat,expt,texpt,YR);
  
	fprintf('Loading %s\n',fmatout);
	load(fmatout);

	nsc = length(SCT);
	if ~exist('ii0','var');
		for ik=1:nsc
			nm=SCT(ik).Name;
			if strncmp(nm,snm,4); ii0=ik; break; end;
		end
	end
%
%  tmy
	tmy = SCT(ii0).Time;
%
	if ~exist('ZZ','var')
		ZZ = SCT(ii0).ZZintrp;
		dZ = abs(diff(ZZ));
		dz1 = 0.5*abs(ZZ(2)-ZZ(1));
		dZ = [dz1;dZ];
		dZ(end) = dZ(end)-dz1;
		dL = SCT(ii0).segm_dL;
		[DL,DZ]=meshgrid(dL,dZ);
	end

% Vol Flux:
  uu = SCT(ii0).Unrm;
  tt = SCT(ii0).T;
	ss = SCT(ii0).S;
	ndays = size(uu,1);

	Vflx = [];
	Sflx = [];
	Tflx = [];
	for iday = 1:ndays
		u1   = squeeze(uu(iday,:,:));
		t1   = squeeze(tt(iday,:,:));
		s1   = squeeze(ss(iday,:,:));
    rhow = sw_dens0(s1,t1);
		vf1  = nansum(u1.*DL.*DZ);
		hf1  = nansum(u1.*Cp.*rhow.*(t1-Tref1).*DL.*DZ);
		fwf1 = nansum(u1.*(Sref1-s1)./Sref1.*DL.*DZ);
		sf1  = nansum(u1.*s1.*DL.*DZ); 

		if isempty(Vflx)
			Vflx = vf1;
			Sflx = sf1;
			Tflx = hf1;
			FWflx= fwf1;
		else
			Vflx = Vflx+vf1;
			Sflx = Sflx+sf1;
			Tflx = Tflx+hf1;
			FWflx= FWflx+fwf1;
		end

	end
	Vflx = Vflx/ndays;
	Sflx = Sflx/ndays;
	Tflx = Tflx/ndays;
	FWflx = FWflx/ndays;
%
% Filter out high-freq oscillation due to +/- fluxes at zigzaging section
  if ~exist('Wn','var')
  	Wn = 1/5;
	  [Bf,Af] = butter(9,Wn,'low');
  end
	VflxF = filtfilt(Bf,Af,Vflx); % 
	SflxF = filtfilt(Bf,Af,Sflx); % 
  TflxF = filtfilt(Bf,Af,Tflx); % 
	FWflxF = filtfilt(Bf,Af,FWflx); % 

	Dist = cumsum(dL);

	if f_save == 1
		nm = SCT(ii0).Name;
    fflux = sprintf('%shycom004_023_Fluxes_%s_%i.dat',pthout,nm,YR);
		fprintf('Writing fluxes --> %s\n',fflux);
		fid = fopen(fflux,'w','ieee-be');
		nfields = 4;
		npnts = length(VflxF);
		fwrite(fid,nfields,'int');
		fwrite(fid,npnts,'int');
		fwrite(fid,dL,'float32');
		fwrite(fid,VflxF,'float32');
		fwrite(fid,TflxF,'float32');
		fwrite(fid,FWflxF,'float32');
		fclose(fid);
  end
	
end








