% Test runs with 0.08 HYCOM-CICE5 GOFS3.5
% to debug sea ice melt problem
%
% Greenland runoff, no passive tracers
% analyze monthly mean sea ice fields
% extract from daily instantenous output
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

%expt = 122;  % GOFS3.5 HYCOM-CICEv5
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.08/datamat/cice_mnth/';
pth1 = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/';
pthmat = pthout;

%texpt = 'dltEddTfrz';
%texpt = 'BL99Tmlt';   % CCSM3 sh/wave radiation, Tmlt > Tfrz ice heat potent.
%texpt = 'BL99Tfrz';   % CSM3 sh/wave radiation, Tmlt=Tfrz (similar to 0.04 original run)
%texpt = 'dltEddTmlt';

iex = 3; % expt # in EXPT - run sub_cice_experiments to see all expts
EXPT = sub_cice_experiments;

expt   = EXPT(iex).Nmb;
pthmat = EXPT(iex).pthmat;
texpt  = EXPT(iex).cice_opt;
mmean  = EXPT(iex).day_mean; % = 0 - instanteneous or daily mean fields
res    = EXPT(iex).res;
vcice  = EXPT(iex).cice;
pthout = EXPT(iex).dir_outp;

mmean = 1; % = 0 - instanteneous or daily mean fields


% start :
yr1   = 2017;
mo1   = 6;
mday1 = 1;
dnmb1 = datenum(yr1,mo1,mday1);

% end:
yr2   = 2017;
mo2   = 12;
mday2 = 31;
dnmb2 = datenum(yr2,mo2,mday2);
dskp  = 4;

fprintf('Extracting time series of CICE %s: %s-%s\n\n',texpt,datestr(dnmb1),datestr(dnmb2));
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
%Lmsk = HH*0;
%Lmsk(HH<0)=1;
Ioc = find(LAT>70. & HH<-50);   % for surf fluxes, Arctic region


DNMB = [dnmb1:dskp:dnmb2];
nrc  = length(DNMB);

for irc=1:nrc
  dnmb = DNMB(irc);
  DV   = datevec(dnmb);
  yr   = DV(1);
  mo   = DV(2);
  mday = DV(3);
  YR   = yr;
  iday = dnmb-datenum(YR,1,1)+1;


  pthbin = sprintf('%sARCc0.08_122_%s/%i_cice/',pth1,texpt,YR);


		if mmean 
				flnm = sprintf('%s122_cice.%i-%2.2i-%2.2i.nc',pthbin,YR,mo,mday);
		else
				flnm = sprintf('%s%3.3i_cice_inst.%i-%2.2i-%2.2i-00000.nc',pthbin,expt,YR,mo,mday);
		end

  ICE(irc).dnmb = dnmb;

  if ~exist(flnm,'file')
    fprintf('Not found %s\n',flnm);
    continue;
  end
  fprintf('Reading %s \n',flnm);

  ai = squeeze(nc_varget(flnm,'aice')); % ice area  (aggregate)
%  I = find(ai>0.15 & LAT>80);
  I = find(ai>0.15 & LAT>75);  % for ice characteristics that make sense over sea ice
  if isempty(I); continue; end;
  ai = ai(I);
  hi = squeeze(nc_varget(flnm,'hi'));   % grid cell mean ice thickness
  hi = hi(I);
  si = squeeze(nc_varget(flnm,'hs'));   % grid cell mean snow thickness
  si = si(I);
%  tsfc = squeeze(nc_varget(flnm,'Tsfc')); % snow/ice surface temperature
%  tsfc = tsfc(I);
%  sice = squeeze(nc_varget(flnm,'sice'));  % bulk ice S
%  sice = sice(I);
  fswdn = squeeze(nc_varget(flnm,'fswdn')); % down solar flux, W/m2
  fswdn = fswdn(Ioc);
  flwdn = squeeze(nc_varget(flnm,'flwdn')); % down longwave flux, W/m2
  flwdn = flwdn(Ioc);
  frzmlt = squeeze(nc_varget(flnm,'frzmlt')); % freeze-melt potential, >0 - ice forms
  frzmlt = frzmlt(I);
% fswint - all zeros?
%  fswint = squeeze(nc_varget(flnm,'fswint_ai')); % sh/wave absorbed in ice interior, W/m2
%  fswint = fswint(I)./ai;  % denomralize - normalized by ice area
  fsnow = squeeze(nc_varget(flnm,'snow_ai')); % snowfall rate
  fsnow = fsnow(I)./ai;
  fswabs = squeeze(nc_varget(flnm,'fswabs_ai')); % snow/ice/ocean absorbed solar flux
  fswabs = fswabs(I)./ai; % denomralize - normalized by ice area
  albsni = squeeze(nc_varget(flnm,'albsni')); % snow/ice broad band albedo
  albsni = albsni(I);
  flat = squeeze(nc_varget(flnm,'flat_ai')); % latent heat flux, weight by ice area W/m2
  flat = flat(I)./ai; % denomralize - normalized by ice area
  fsens = squeeze(nc_varget(flnm,'fsens_ai')); % sensible heat flux, weighted by ice area
  fsens = fsens(I)./ai; % denomralize - normalized by ice area
%  congel = squeeze(nc_varget(flnm,'congel')); % congelation ice grwoth, cn/day
%  congel = congel(I);
  meltt = squeeze(nc_varget(flnm,'meltt')); % top ice melt, cm/day
  meltt = meltt(I);
  melts = squeeze(nc_varget(flnm,'melts')); % top snow melt, cm/dat
  melts = melts(I);
  meltb = squeeze(nc_varget(flnm,'meltb')); % basal ice melt
  meltb = meltb(I);
%  meltl = squeeze(nc_varget(flnm,'meltl')); % lateral ice melt
%  meltl = meltl(I);
  fhocn = squeeze(nc_varget(flnm,'fhocn_ai')); % heat flux ice to ocean, W/m2
  fhocn = fhocn(I)./ai; % denomralize - normalized by ice area
  fswthru = squeeze(nc_varget(flnm,'fswthru_ai')); % SW flux thru ice to ocean
  fswthru = fswthru(I)./ai; % denomralize - normalized by ice area
  sst   = squeeze(nc_varget(flnm,'sst')); % sea surf T
  sst   = sst(I);


  ICE(irc).ai = ai;
  ICE(irc).hi = hi;
  ICE(irc).si = si;
%  ICE(irc).tsfc = tsfc;
%  ICE(irc).sice = sice;
  ICE(irc).fswdn = fswdn;
  ICE(irc).flwdn = flwdn;
  ICE(irc).frzmlt = frzmlt;
%  ICE(irc).fswint = fswint;
  ICE(irc).fsnow = fsnow;
  ICE(irc).fswabs = fswabs;
  ICE(irc).albsni = albsni;
  ICE(irc).flat = flat;
  ICE(irc).fsens = fsens;
%  ICE(irc).congel = congel;
  ICE(irc).meltt = meltt;
  ICE(irc).melts = melts;
  ICE(irc).meltb = meltb;
%  ICE(irc).meltl = meltl;
  ICE(irc).fhocn = fhocn;
  ICE(irc).fswthru = fswthru;
  ICE(irc).sst   = sst;
  ICE(irc).Indx = I;

end;

fmatout = sprintf('%s%3.3icice_outp_tser_%s.mat',pthmat,expt,texpt);
fprintf('Saving %s\n',fmatout);
save(fmatout,'ICE','-v7.3');


  















