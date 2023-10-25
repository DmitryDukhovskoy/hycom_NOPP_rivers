% Calculate monthly mean conc/h CICE output - 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt=110;
s_mat = 1;

yr1=1993;
yr2=2016;
mo=9;
dd=7; 

%pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthbin  = '/nexsan/people/ddmitry/hycom/ARCc0.08/110/cice/';
%pthbin = sprintf('/nexsan/archive/ARCc0.08_%s/data/%i_cice/',expt,yr);
%pthbin = '/nexsan/hycom/ARCc0.08_112/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_icemat/';

fmat    = sprintf('%sarc08_%3.3i_cice_mean_mnth%2.2i.mat',pthmat,expt,mo);

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);

% Domain: AO and SPNA
i1=35;
j1=246;
i2=nn;
j2=1916;


cc=0;
for yr=yr1:yr2
  pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%i_cice/',expt,yr);

  cdd=0;
  h0=[];
  a0=[];
  u0=[];
  v0=[];
  for mday=2:dd:31
    tic;
    
    fin = sprintf('%s%3.3i_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',...
     	         pthbin,expt,yr,mo,mday);
    if ~exist(fin,'file')
      fprintf('%s not found \n',fin);
      continue;
    end
    
    dnmb = datenum(yr,mo,mday);
    iday = dnmb-datenum(yr,1,1)+1;
    fprintf('Extracting %i/%i/%i\n',yr,mo,mday);
    
    hi=squeeze(nc_varget(fin,'hi'));
    ai=squeeze(nc_varget(fin,'aice'));
    ui=squeeze(nc_varget(fin,'uvel'));
    vi=squeeze(nc_varget(fin,'vvel'));
    hi=hi(j1:j2,i1:i2);
    ai=ai(j1:j2,i1:i2);
    ui=ui(j1:j2,i1:i2);
    vi=vi(j1:j2,i1:i2);
    
    cdd=cdd+1;
    if cdd==1
      h0=hi*0;
      a0=ai*0;
      u0=ui*0;
      v0=vi*0;
    end
    h0=h0+hi;
    a0=a0+ai;
    u0=u0+ui;
    v0=v0+vi;
    
    fprintf('1 record %4.1f sec\n',toc);
    
  end  % days for given month
  
  cc=cc+1;
  ICE(cc).Year=yr;
  ICE(cc).hice=h0./cdd;
  ICE(cc).aice=a0./cdd;
  ICE(cc).uice=u0./cdd;
  ICE(cc).vice=v0./cdd;
  
  if s_mat==1
    fprintf('Saving %s\n',fmat);
    save(fmat,'ICE');
  end
  
  
end    % years
