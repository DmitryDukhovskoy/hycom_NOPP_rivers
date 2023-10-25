% Experiments with CICEv5
%
% Extract time series of sea ice area & volume
% from specified experiments
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

pthout   = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat   = pthout;
pthtopo  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo4 = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';


f_add = 8;  % what experiment to add to existing saved, = 0 - all experiments

EXPT{1} = '0.08_122_dltEddTmlt';
EXPT{2} = '0.08_122_dltEddTfrz';
EXPT{3} = '0.08_122_BL99Tmlt';
EXPT{4} = '0.08_122_BL99Tfrz';
EXPT{5} = '0.04_BL99Tfrz';               % Old run 0.04 used BL99 Tfrz
EXPT{6} = '0.04_dltEddTmlt';             % New run 0.04 with delta Eddington sh/w and Tmlt  
EXPT{7} = '0.04_dltEddTmltEAP';          % dlt Edd, Tmlt, EAP rhoeology
EXPT{8} = '0.04_dltEddTmltEAPJRA';       % JRA55 forcing expt_023
EXPT{9} = '0.08_dltEddTmltEAPJRA';       % 0.08 JRA forc, similar to 0.04 expt 02.3

ftopo = sprintf('%s/depth_ARCc0.04_17DD.nc',pthtopo4); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell04=DX.*DY; % Grid cell area, m2


ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;

[DX,DY]=sub_dx_dy(LON,LAT);
Acell08=DX.*DY; % Grid cell area, m2

% Extract 1 year at a time
dd = 5; 
YR = 2019;
d1 = datenum(YR,1,1);
d2 = datenum(YR,12,31);

% Output:
%fmatout = sprintf('%sseaice_VolArea_tests.mat',pthmat);
fmatout = sprintf('%sseaice_VolArea_tests%i.mat',pthmat,YR);

if f_add >0; % add to existing saved experiments

% If file does not exist - create empy array and add only 1 experiment
  if exist(fmatout,'file')
    fprintf('Loading %s\n',fmatout);
    load(fmatout);
  else
    fprintf('Does not exist %s\n',fmatout);
    fprintf('Only %s will be saved\n',EXPT{f_add});
    ixx = f_add;
    ICE(ixx).TM       = [];
    ICE(ixx).Area_km2 = [];
    ICE(ixx).IExt_km2 = [];
    ICE(ixx).Vol_km3  = [];
  end
 
end


nexpt = length(EXPT);
for ixx = 1:nexpt

  if f_add>0 & ixx~=f_add, continue; end;

  nme = EXPT{ixx};
  if contains(nme,'0.04')
    if strncmp(nme,'0.04_BL99Tfrz',13)
      pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i_cice_BL99Tfrz/',YR);
      fl0 = '022_cice';
    elseif contains(nme,'JRA')
      pthbin = ...
       sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_023/data/%i_cice_dltEddTmltEAPJRA/',YR);
      fl0 = '023_cice';
    elseif strncmp(nme,'0.04_dltEddTmlt',15)
      pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i_cice_dltEddTmlt/',YR);
      fl0 = '022_cice';
    elseif strncmp(nme,'0.04_dltEddTmltEAP',18);
      pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i_dltEddTmltEAP/',YR);
      fl0 = '022_cice';
    end
    Acell = Acell04;
  else
    if contains(nme,'BL99Tfrz')
      pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.08_122/data/%i_cice_BL99Tfrz/',YR);
      fl0 = '122_cice';
    elseif contains(nme,'BL99Tmlt')
      pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.08_122/data/%i_cice_BL99Tmlt/',YR);
      fl0 = '122_cice';
    elseif contains(nme,'dltEddTfrz')
      pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.08_122/data/%i_cice_dltEddTfrz/',YR);
      fl0 = '122_cice';
    elseif contains(nme,'dltEddTmlt') & ~contains(nme,'JRA')
      pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.08_122/data/%i_cice_dltEddTmlt/',YR);
      fl0 = '122_cice';    
    elseif contains(nme,'JRA')
      pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.08_123/data/%i_cice_dltEddTmltEAPJRA/',YR);
      fl0 = '123_cice';
    end
    Acell = Acell08;
  end

  ICE(ixx).Name = nme;
%keyboard

  icc = 0;
  for dnmb = d1:dd:d2
    dv   = datevec(dnmb);
    mo   = dv(2);
    mday = dv(3);

    fprintf('Reading %s %s\n',nme,datestr(dnmb));
    flnm = sprintf('%s%s.%i-%2.2i-%2.2i.nc',pthbin,fl0,YR,mo,mday);

    if ~exist(flnm,'file');
      fprintf('Not found %s\n',flnm);
      continue;
    end

    fprintf('Reading %s \n',flnm);
    Hi = squeeze(nc_varget(flnm,'hi'));
    Ci = squeeze(nc_varget(flnm,'aice'));

		cice0=1.e-20;  % for ice extent - should be 0.15, for ice area - all ice >0% is considered
%    Ic = find(Ci>=0.15); % do not need this for the area, only for extent!
% also for extent - make conc > cice0 conc=1 - area of grid cell with ice
		Ic = find(Ci>=cice0);
    cci = Ci;
    cci(Ci<cice0)=0;
    cciE = Ci;     % something similar to ice area extent to compare to NOAA
    cciE(Ci<cice0)=0;
    cciE(Ci>=cice0)=1.;
    Aice = nansum(nansum(Acell.*cci*1e-6)); % km2 - ice area
    Iext = nansum(nansum(Acell.*cciE*1e-6)); % km3 - ice Ext. 
    Vice = nansum(nansum(Acell.*cci.*Hi*1e-9));  % km3

    fprintf('  A=%4.1fe6 km2, IExt=%4.1fe6 km2, Vol=%6.4d km3\n',...
            Aice*1e-6',Iext*1e-6,Vice);
    icc = icc+1;
    ICE(ixx).TM(icc)       = dnmb;
    ICE(ixx).Area_km2(icc) = Aice;
    ICE(ixx).IExt_km2(icc) = Iext;
    ICE(ixx).Vol_km3(icc)  = Vice;


  end
end

fprintf('Saving %s\n',fmatout);
save(fmatout,'ICE');


