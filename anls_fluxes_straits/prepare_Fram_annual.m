% For Laura prepare monthly climatology
% of Fram Strait flux
% prepared in extr_TSVdaily_straits08.m
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=112;
TV=11;
YR1=2015;
YR2=2019;
dday=7;

snm = 'FramStr';
s_mat=1;  % ==2 - start from last saved rec # in YEAR

Cp = 4200; % J/kg K
Tref1 = -1.8; % Ref T to calc. H flux
Tref2 = 0;    % Ref T to calc. H flux
Sref1 = 34.8;
Sref2 = 34.9;
hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation

pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='extr_TSVdaily_straits08.m';

fprintf('arc08-%3.3i Heat and Vol fluxes Greenland Shelf gates %i-%i, save=%i\n',...
        expt,YR1,YR2,s_mat);


icc=0;
YR1=2005;
YR2=2019;
for iyr=YR1:YR2
  Usm=[];
  Tsm=[];
  Ssm=[];

  yr=iyr;
  YR=iyr;
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

% Net Vol flux:
  icc=icc+1;
		dmm = SCT(ii0).VolFlx_m3s;
  Unrm = SCT(ii0).Unrm;
  U=squeeze(mean(Unrm,1));
  SS = SCT(ii0).S;
  S=squeeze(mean(SS,1));

  II = SCT(ii0).I;
  JJ = SCT(ii0).J;

  FRAM(icc).Long = SCT(ii0).long;
  FRAM(icc).Lat = SCT(ii0).latd;
  FRAM(icc).Iindx= II;
  FRAM(icc).Jindx= JJ;
  FRAM(icc).Z_depth = SCT(ii0).ZZintrp;
  FRAM(icc).Year = iyr;
		FRAM(icc).VolFlux_m3s = mean(dmm);
  FRAM(icc).VolFlux_err_m3s = std(dmm)/sqrt(length(dmm));
  FRAM(icc).NormalU_ms = U;
  FRAM(icc).Saln = S;

end

if s_mat==1
  fgnm = sprintf('%sHYCOM_CICE_008_%3.3i_%s_%i-%i.mat',pthout,expt,snm,YR1,YR2);
  fprintf('Saving %s\n',fgnm);
  save(fgnm,'FRAM');
end





