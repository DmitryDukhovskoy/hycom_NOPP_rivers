% Plot wndewd wndnwd *.r files focring cice
% C  hycom2raw8 - Usage:  hycom2raw8 fhycom.a idm jdm [i1 j1 idms jdms] [spval] fraw.a
%C
%C  Outputs a raw8 (no control words, no padding) subarray.
%C
%C  fhycom.a is assumed to contain idm*jdm 32-bit IEEE real values for
%C   each array, in standard f77 element order, followed by padding
%C   to a multiple of 4096 32-bit words, but otherwise with no control
%C   bytes/words, and input values of 2.0**100 indicating a data void.
%C
%C  The output is the sub-array (i1:i1+idms-1,j1:j1+jdms-1), with data
%C   voids where this extends outside (1:idm,1:jdm).
%C  fraw.a will contain idms*jdms 64-bit IEEE real values for each array,
%C   in standard f77 element order, with no control words, no padding,
%C   and data voids indicated by spval (default 2.0**100).
%C
%C  this version for "serial" Unix systems.


addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

regn = 'ARCc0.04';

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

%mm  = 5040;
%nn  = 3200;
JDM = mm;
IDM = nn;
IJDM=IDM*JDM;
%npad=4096-mod(IJDM,4096);
%toto=ones(npad,1);

pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/cice/';
fin = sprintf('%swndewd_116b.r',pthbin);

fid = fopen(fin,'r','ieee-be');

frewind(fid);
% Hourly fields: max # records 
iday=17;
ihr=0;
imp=(iday-1)*24+ihr+1; % record to plot
% For CICE forcing fields are not fixed-length records
% as for HYCOM
%k0 = imp-1;
%stat=fseek(fid,k0*(IJDM+npad)*8,-1);
%if abs(stat)>0
%  fprintf('Error finding record %i\n',imp);
%  error('Stopping');
%end

% Estimate # of records in the file:
% get size of the *r file:
nrec = 89929728000/(IDM*JDM*8);


% For CICE forcing fields are not fixed-length records
% as for HYCOM
for ii=1:imp
  fprintf('Reading record %i, need %i\n',ii,imp);
%  imm=fread(fid,1,'uint32');
  dmm=fread(fid,[IDM,JDM],'float64');
%  imm=fread(fid,1,'uint32');
  if (isempty(dmm)); 
    fprintf('E-o-F, month = %i\n',ii);
  end
end
A04 = dmm';


ikl = max(findstr(fin,'/'));
fnm = fin(ikl+1:end);

figure(1); clf;
contour(HH,[0 0],'k');
hold on;
pcolor(A04); shading flat;
colorbar
title(sprintf('%s record %i, %i:%2.2ih',fnm,imp,iday,ihr),'Interpreter','None');


fclose(fid);

btx = 'plot_wind_cice.m';
bottom_text(btx,'pwd',1);



