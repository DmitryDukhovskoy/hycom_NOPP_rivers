% CICE5 is mapd from GLBc0.04 GOFS3.5 CICE5
% using D. Hebert code
% on Gordon
%  022rstrt_cice5_preprocess.com
%
% Change restart date/time in global attributes
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear


pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthdat = '/Net/kronos/ddmitry/hycom/ARCc0.04/cice/';

%frst = sprintf('%scice.restart.2017010109_interp.nc',pthdat);
frst = sprintf('%scice.restart.2017090109_interp.nc',pthdat);

%ftopo = sprintf('%s/depth_ARCc0.04_17DD.nc',pthtopo); % 
%HH  = nc_varget(ftopo,'Bathymetry');
%LON = nc_varget(ftopo,'Longitude');
%LAT = nc_varget(ftopo,'Latitude');
%[mm,nn]=size(LON);

%
%uu = nc_varget(frst,'uvel');

%myrstd = datenum(2016,1,1);
myrstd = datenum(2016,9,1);
hycomd = get_dtime(myrstd);


% Check the restart date
% in the existing file:
time_in = ncreadatt(frst,'/','time');
tfrc_in = ncreadatt(frst,'/','time_forc');
%
% Convert to HYCOM days:
hday_in = time_in/86400+1;

hnyr_in = ncreadatt(frst,'/','nyr');
mnth_in = ncreadatt(frst,'/','month');
mday_in = ncreadatt(frst,'/','mday');
dsec_in = ncreadatt(frst,'/','sec');
hrst_in = dsec_in/3600; 

fprintf('Global info %s\n',frst);
fprintf('Restart HYCOMtime:     %12.5e\n',time_in);
fprintf('Restart forcetime:     %12.5e\n',tfrc_in);
fprintf('Restart HYCOM day:     %10.3f\n',hday_in);
fprintf('Restart HYCOM year:    %i\n',hnyr_in);
fprintf('Restart month/day/HR : %2.2i/%2.2i/%5.2f\n\n',mnth_in,mday_in,hrst_in);


% ===================
fprintf('Changing to %s\n',datestr(myrstd));
dv = datevec(myrstd);
hnyr_out = dv(1)-1900;
mnth_out = dv(2);
mday_out = dv(3);
hrst_out = 0;
dsec_out = hrst_out*3600;
hday_out = get_dtime(myrstd);
time_out = (hday_out-1)*86400;
tfrc_out = time_out;

fprintf('Global attributes will be changed to:\n');
fprintf('==> Restart HYCOMtime:    %12.5e\n',time_out);
fprintf('==> Restart forcetime:    %12.5e\n',tfrc_out);
fprintf('==> Restart HYCOM day:     %10.3f\n',hday_out);
fprintf('==> Restart HYCOM year:    %i\n',hnyr_out);
fprintf('==> Restart month/day/HR : %2.2i/%2.2i/%5.2f\n\n',mnth_out,mday_out,hrst_out);

keyboard

if time_out == time_in
  fprintf('Time attributes are the same, nothing to change\n');
  return
end

fprintf('Writing time: %i\n',time_out);
tic;
ncwriteatt(frst,'/','time',time_out);
fprintf(' Done %6.2f min\n',toc/60);


fprintf('Writing time_forc: %i\n',tfrc_out);
tic;
ncwriteatt(frst,'/','time_forc',tfrc_out);
fprintf(' Done %6.2f min\n',toc/60);

if hnyr_out ~= hnyr_in
		fprintf('Writing nyr: %i\n',hnyr_out);
		tic;
		ncwriteatt(frst,'/','nyr',hnyr_out);
		fprintf(' Done %6.2f min\n',toc/60);
end

if mnth_out ~= mnth_in
  fprintf('Writing month: %i\n',mnth_out);
  tic;
  ncwriteatt(frst,'/','month',mnth_out);
  fprintf(' Done %6.2f min\n',toc/60);
end

if mday_out ~= mday_in
  fprintf('Writing mday: %i\n',mday_out);
  tic;
  ncwriteatt(frst,'/','mday',mday_out);
  fprintf(' Done %6.2f min\n',toc/60);

end

if dsec_out ~=dsec_in
  fprintf('Writing sec: %i\n',dsec_out);
  tic;
  ncwriteatt(frst,'/','sec',dsec_out);
  fprintf(' Done %6.2f min\n',toc/60);
end

fprintf('All done \n');





