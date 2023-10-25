% For Laura prepare daily snapshots
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

btx='prepare_Fram_daily.m';

fprintf('arc08-%3.3i Heat and Vol fluxes Greenland Shelf gates %i-%i, save=%i\n',...
        expt,YR1,YR2,s_mat);

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);



icc=0;
YR1=2019;
YR2=2019;
for iyr=YR1:YR2
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
  TM = SCT(ii0).Time;
  DV = datevec(TM);

  Unrm = SCT(ii0).Unrm;
  SS = SCT(ii0).S;
  VolFlx=SCT(ii0).VolFlx_m3s;

  II = SCT(ii0).I;
  JJ = SCT(ii0).J;

  FRAM.NormalU_ms = Unrm;
  FRAM.Saln = SS;
  FRAM.VolFlux_m3s = VolFlx;
  FRAM.Time_mat = TM;


		FRAM.Long  = SCT(ii0).long;
		FRAM.Lat   = SCT(ii0).latd;
		FRAM.Iindx = II;
		FRAM.Jindx = JJ;
		FRAM.Z_depth = SCT(ii0).ZZintrp;
		FRAM.Years = [YR1:YR2];
		for j=1:length(JJ)
				i0=II(j);
				j0=JJ(j);
				FRAM.Hbottom(j,1) = HH(j0,i0);
		end


		if s_mat==1
				fgnm = sprintf('%sHYCOM_CICE_008_%3.3i_%s_daily%i.mat',pthout,expt,snm,YR);
				fprintf('Saving %s\n',fgnm);
				save(fgnm,'FRAM');
		end

end


f_chck=0;
if f_chck
  Z   = FRAM.Z_depth;
  Lon = FRAM.Long;
  id = 5;
  S = squeeze(FRAM.Saln(id,:,:));
  U  = squeeze(FRAM.NormalU_ms(id,:,:));
  Hb = FRAM.Hbottom;
  dnmb = TM(id);

  figure(1); clf;
  axes('Position',[0.09 0.3 0.86 0.62]);
  hold on;
  pcolor(Lon,Z,S); shading flat;
%  pcolor(Lon,Z,U); shading flat;
  caxis([34 35]);
  plot(Lon,Hb,'k-');
  set(gca,'xlim',[-20 14],...
          'ylim',[-3500 0]);
  colorbar('SouthOutside');
  stl=sprintf('S %s',datestr(dnmb));
  title(stl);

  bottom_text(btx,'pwd',1);
end





