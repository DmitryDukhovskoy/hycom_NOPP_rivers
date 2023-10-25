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

btx='prepare_Fram_clim.m';

fprintf('arc08-%3.3i Heat and Vol fluxes Greenland Shelf gates %i-%i, save=%i\n',...
        expt,YR1,YR2,s_mat);

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);


for im=1:12
  FRAM(im).nrecs=0;
  FRAM(im).NormalU_ms = [];
  FRAM(im).Saln = [];
  FRAM(im).VolFlux_m3s = [];
end;


icc=0;
YR1=2005;
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

  for im=1:12
    Im=find(DV(:,2)==im);
    cnt = FRAM(im).nrecs;   
    if cnt==0
						FRAM(im).Saln = squeeze(mean(SS(Im,:,:)));
						FRAM(im).NormalU_ms = squeeze(mean(Unrm(Im,:,:)));
						FRAM(im).VolFlux_m3s = +mean(VolFlx(Im));
    else
						um = FRAM(im).NormalU_ms;
						sm = FRAM(im).Saln;
						vfl= FRAM(im).VolFlux_m3s;
						FRAM(im).Saln = sm+squeeze(mean(SS(Im,:,:)));
						FRAM(im).NormalU_ms = um+squeeze(mean(Unrm(Im,:,:)));
						FRAM(im).VolFlux_m3s = vfl+mean(VolFlx(Im));
    end

    FRAM(im).nrecs=cnt+1;
  end;

end

FRAM(1).Long  = SCT(ii0).long;
FRAM(1).Lat   = SCT(ii0).latd;
FRAM(1).Iindx = II;
FRAM(1).Jindx = JJ;
FRAM(1).Z_depth = SCT(ii0).ZZintrp;
FRAM(1).Years = [YR1:YR2];
for j=1:length(JJ)
  i0=II(j);
  j0=JJ(j);
  FRAM(1).Hbottom(j,1) = HH(j0,i0);
end

for im=1:12
  cnt = FRAM(im).nrecs;
  uu  = FRAM(im).NormalU_ms;
  FRAM(im).NormalU_ms = uu./cnt;

  ss = FRAM(im).Saln;
  FRAM(im).Saln = ss./cnt;

  vfl = FRAM(im).VolFlux_m3s;
  FRAM(im).VolFlux_m3s = vfl./cnt;
end;


if s_mat==1
  fgnm = sprintf('%sHYCOM_CICE_008_%3.3i_%s_mnth%i-%i.mat',pthout,expt,snm,YR1,YR2);
  fprintf('Saving %s\n',fgnm);
  save(fgnm,'FRAM');
end

f_chck=0;
if f_chck
  Z   = FRAM(1).Z_depth;
  Lon = FRAM(1).Long;
  im = 5;
  S = FRAM(im).Saln;
  U  = FRAM(im).NormalU_ms;
  Hb = FRAM(1).Hbottom;

  figure(1); clf;
  axes('Position',[0.09 0.3 0.86 0.62]);
  hold on;
  pcolor(Lon,Z,S); shading flat;
%  pcolor(Lon,Z,U); shading flat;
  caxis([33 35]);
  plot(Lon,Hb,'k-');
  set(gca,'xlim',[-20 14],...
          'ylim',[-3500 0]);
  colorbar('SouthOutside');
  stl=sprintf('S mo=%i',im);
  title(stl);

end





