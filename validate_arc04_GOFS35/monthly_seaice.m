% Arctic Ocean 0.04 HYCOM-CICEv5 GOFS3.5
% Greenland runoff, no passive tracers
% analyze monthly mean sea ice fields
% extract from daily instantenous output
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 22;
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat = pthout;

yr1=2017;
yr2=2017;

for YR=yr1:yr2
  pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%i_cice/',expt,YR);

  ICEMO = struct;
  ICEMO.Info = 'Monthly Sea Ice Fields from 0.04 HYCOM-CICEv5 GOFS3.5';
  ICEMO.hycom_expt = sprintf('%3.3i',expt);

  for im=1:12
    fmatout = sprintf('%scice_monthly_%i%2.2i.mat',pthmat,YR,im);

    d1=datenum(YR,im,1);
    d2=d1+32;
    dv=datevec(d2);
    dmm=datenum(dv(1),dv(2),1);
    d2=dmm-1;
    id1=1;
    id2=d2-d1+1;

    Hice = [];
    Cice = [];
    Sice = [];
    Uice = [];
    Vice = [];
    icc  = 0;
    for iday=id1:id2
      flnm = sprintf('%s022_cice_inst.%i-%2.2i-%2.2i-00000.nc',pthbin,YR,im,iday);

      if ~exist(flnm,'file')
        fprintf('Not found %s\n',flnm);
        continue;
      end
  
      tic;

      fprintf('Reading %s \n',flnm);
      Hi = squeeze(nc_varget(flnm,'hi'));
      Si = squeeze(nc_varget(flnm,'hs'));
      Ci = squeeze(nc_varget(flnm,'aice'));
      Ui = squeeze(nc_varget(flnm,'uvel'));
      Vi = squeeze(nc_varget(flnm,'vvel'));
      Hi(isnan(Hi))=0;
      Si(isnan(Si))=0;
      Ci(isnan(Ci))=0;
      Ui(isnan(Ui))=0;
      Vi(isnan(Vi))=0;

      icc=icc+1;
      if isempty(Hice)
        Hice = Hi;
        Sice = Si;
        Cice = Ci;
        Uice = Ui;
        Vice = Vi;
      else
        Hice = Hice+Hi;
        Sice = Sice+Si;
        Cice = Cice+Ci;
        Uice = Uice+Ui;
        Vice = Vice+Vi;
      end;
%keyboard
      Hmin = min(min(Hice(Hice>0)))/icc;
      Hmax = max(max(Hice))/icc;
      Smin = min(min(Sice(Sice>0)))/icc;
      Smax = max(max(Sice))/icc;

      fprintf('H ice  min/max = %8.4f %8.4f m\n',Hmin,Hmax);
      fprintf('H snow min/max = %8.4f %8.4f m\n',Smin,Smax);
      fprintf('1 day processed %6.4f min\n',toc/60);
    end

    Hice = Hice/icc;
    Sice = Sice/icc;
    Cice = Cice/icc;
    Uice = Uice/icc;
    Vice = Vice/icc;

    ICEMO.Ice_thck  = Hice;
    ICEMO.Snow_thck = Sice;
    ICEMO.Ice_conc  = Cice;
    ICEMO.Ice_u     = Uice;
    ICEMO.Ice_v     = Vice;

    if icc>1
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'ICEMO');
    end
  end
end




