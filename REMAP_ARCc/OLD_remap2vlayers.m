% Test code for remapping
% 32 v layers from GLBb0.08 T07 experiment 19.0
% onto 41 v layers ARCc0.07 T11
% The first 2 steps (GLB-> ARC and ARC T07 ->T11) 
% should be already done
% This code uses ARCc T11 32 layer archive files
% created from GLBb0.08 T07
% and interpolates into 41 layers used in GOFS3.1
%
% Use layer depths from PHC climatology file
% for 41 layers (see: RELAX_PHC/*)
%
% Targ. densities of these 41 layers are such that
% added layers are at the top of the ocean
% above sigma2=35.50 (L41=24) all layers beneath it
% match 32-layer grid
% The added layers are anticipated to be fixed depth layers 
% Inspection of hybrid layers generated for HPC climatology
% indicates that the 2 last layers above sigma2=35.50 may be 
% not fixed-depth
% 
% This algorithm assumes:
% 1) All layers L41<24 (sigma2=35.50) are fixed depth
% 2) Depth of the bottom interface of L41(24) = L32(14) 
%    i.e. sum[i=1,24]L41(i)=sum[i=1,14]L32(i), however 
%    at some locations 32-layer 14 is shallower than
%    L41(24) and even L41(23) - will need to adjust this
% 
% 3) The bottom layer L32(32), sigma2=37.48 does not exist in L41
%    will need to merge the last two L32 layers in 1 for L41

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

%pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
pthrlx = '/Net/kronos/ddmitry/hycom/ARCc0.08/relax_41layers/output/';
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
pthclim = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/relax/110/';
pthglb  = pthbin;
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

rg=9806;  % convert pressure to depth, m
yr=1993;
iday=1;
hr=0;
TV=11;
nsigma=14; % # of sigma levels
%mo=7; % month, 1 or 7

% Topo:
fltopo=sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn] = size(HH);
%HH(HH>0)=nan;


ID=n;
JD=m;
% Read dP41 layer thicknesses from
% climatology:
finRa = sprintf('%s41layers_T%2.2i/relax.m%2.2i.a',pthbin,TV,mo);
finRb = sprintf('%s41layers_T%2.2i/relax.m%2.2i.b',pthbin,TV,mo);
fld='thknss';
[Fg,nr,mr,lr] = read_hycom(finRa,finRb,fld);
lr=size(Fg,1);

Fg(Fg>1e10)=nan;
Fg(Fg<0.1)=nan;

Fg=Fg./rg;
DP=Fg;

Lxx=24; % where 41 and 32-lr grids should match

% Find points not in sigma-region, i.e.
% where depth> sum over (nsigma dP)
sgmDP=squeeze(nansum(DP(1:nsigma,:,:),1));

% Fix last layer above Lxx
sgmDP(sgmDP==0)=nan;
HH(HH>0)=nan;
sgmH  = abs(1-abs(sgmDP)./abs(HH)); % sigma-coord. points
Isgm=find(sgmH<1e-3); % sigma points
%[js,is]=ind2sub(size(HH),Isgm);
sgmH(Isgm)=nan;
Insgm = find(~isnan(sgmH)); % not-sigma coord. points

dP=squeeze(DP(Lxx-1,:,:));
dP(Isgm)=nan;

figure(1); clf;
contour(HH,[0 0],'k');
hold on;
pcolor(dP); shading flat;



% Input:
fina = sprintf('%sarchv_arcT%2.2i.%4.4i_%3.3i_%2.2i.a',pthbin,TV,yr,iday,hr);
finb = sprintf('%sarchv_arcT%2.2i.%4.4i_%3.3i_%2.2i.b',pthbin,TV,yr,iday,hr);








