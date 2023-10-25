% Check cice atm. forcing fields
% Climatology precipitation
% cice.rhoa_ncar file
%  mm/mo

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.04';
E = '010';
%ntopo1=09;
%ntopo1=11; % topo for ARCc0.08
ntopo2=17; % topo for ARCc0.04
TV = sprintf('%2.2iDD',ntopo2);

ptharc    = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
pthin     = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/force/';
pthout    = sprintf('/Net/mars/ddmitry/hycom/%s/force/',R);
%pthtopo   = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo04 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);


% New topo
fltopo_new=sprintf('%sdepth_%s_%s.nc',pthtopo04,R,TV);
HH   = nc_varget(fltopo_new,'Bathymetry');
LAT = nc_varget(fltopo_new,'Latitude');
LON = nc_varget(fltopo_new,'Longitude');
[mm,nn]= size(HH);

JDM = mm;
IDM = nn;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

fprintf('%s domain, Topo=%s, ID=%i JD=%i\n',R,TV,IDM,JDM);

% -------------------
%   Read precip data
% -------------------
fin = sprintf('%scice.rhoa_ncar85-88_12_arc04.r',pthout);
fid = fopen(fin,'r','ieee-be');

frewind(fid);

imp=12; % month to plot
% For CICE forcing fields are not fixed-length records
% as for HYCOM
for ii=1:imp
  fprintf('Reading , ii = %i\n',ii);
  dmm=fread(fid,[IDM,JDM],'float64');
  if (isempty(dmm)); 
    fprintf('E-o-F, month = %i\n',ii);
    keyboard
  end
end
A04 = dmm';

figure(1); clf;
contour(HH,[0 0],'k');
hold on;
pcolor(A04); shading flat;
colorbar
title('Rho Air, NCAR');

fclose(fid);










