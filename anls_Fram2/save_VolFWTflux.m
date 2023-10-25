% Save depth-averaged fluxes 
% Extracted T, S, normal U at the specified straits
% for processing in python
%  Corrected flux calculation with T,S,U allocation is used
% see anls_fluxes/extr_TSVdaily_straits08.m
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear

f_save = 1;
snm = 'FramStr';

expt=112;
TV=11;
YR1=2005;
YR2=2019;
dday=7;


pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';


%tmm=[datenum(2008,1,1):7:datenum(2008,12,31)];
%tmm=tmm(:);
%for isc=1:5
%  SCT(isc).Time=tmm;
%end;

Cp = 4200; % J/kg K
Tref1 = -1.8; % Ref T to calc. H flux
Tref2 = 0;    % Ref T to calc. H flux
Sref1 = 34.8;
Sref2 = 34.9;
hgg=1e20;




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

	fmatout=sprintf('%shycom008_%3.3i_StraitFluxesDay_%4.4i.mat',...
									pthmat,expt,YR);
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
	Wn = 1/5;
	[Bf,Af] = butter(9,Wn,'low');
	VflxF = filtfilt(Bf,Af,Vflx); % 
	SflxF = filtfilt(Bf,Af,Sflx); % 
  TflxF = filtfilt(Bf,Af,Tflx); % 
	FWflxF = filtfilt(Bf,Af,FWflx); % 

	Dist = cumsum(dL);

	if f_save == 1
		nm = SCT(ii0).Name;
%		fflux = sprintf('%sFluxes_%s_%i.dat',pthout,nm,YR);
		fflux = sprintf('%shycom008_%3.3i_Fluxes_%s_%i.dat',pthout,expt,nm,YR);
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








