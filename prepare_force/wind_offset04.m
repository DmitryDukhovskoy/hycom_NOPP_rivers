% Create 0 offset 
% for wind fields
% ARCc0.04

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.04';
E = '010';
%ntopo1=09;
ntopo1=11; % topo for ARCc0.08
ntopo2=17; % topo for ARCc0.04
TV = sprintf('%2.2iDD',ntopo2);

ptharc    = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
pthin     = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/force/';
pthout    = sprintf('/Net/mars/ddmitry/hycom/%s/force/',R);
pthtopo   = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo04 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);

fltopo_new=sprintf('%sdepth_%s_%s.nc',pthtopo04,R,TV);
fprintf('Reading topo: %s\n',fltopo);
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

fouta = sprintf('%stauewd_zero_arc04.a',pthout);
foutb = sprintf('%stauewd_zero_arc04.b',pthout);
fida  = fopen(fouta,'w','ieee-be');
fidb  = fopen(foutb,'wt');
fprintf('Writing files %s\n',fouta);

A = zeros(nn*mm,1);
nbt = nn*mm*4; % bytes, record lenght

fwrite(fida,nbt,'int');
fwrite(fida,A,'float32');
fwrite(fida,nbt,'int');

aa = 'Zero Stress Offset ARCc0.04                  ';
fprintf(fidb,[aa,'\n']);
aa = '                                             ';
fprintf(fidb,[aa,'\n']);
aa = '                                             ';
fprintf(fidb,[aa,'\n']);
aa = '                                             ';
fprintf(fidb,[aa,'\n']);
aa = sprintf('i/jdm = %i %i',IDM, JDM);
fprintf(fidb,[aa,'\n']);
aa = ' tau_ewd: month,range = 01   0.0000000E+00   0.0000000E+00';
fprintf(fidb,[aa,'\n']);

fclose(fida);
fclose(fidb);


fouta = sprintf('%staunwd_zero_arc04.a',pthout);
foutb = sprintf('%staunwd_zero_arc04.b',pthout);
fida  = fopen(fouta,'w','ieee-be');
fidb  = fopen(foutb,'wt');
fprintf('Writing files %s\n',fouta);

fwrite(fida,nbt,'int');
fwrite(fida,A,'float32');
fwrite(fida,nbt,'int');

aa = 'Zero Stress Offset ARCc0.04                  ';
fprintf(fidb,[aa,'\n']);
aa = '                                             ';
fprintf(fidb,[aa,'\n']);
aa = '                                             ';
fprintf(fidb,[aa,'\n']);
aa = '                                             ';
fprintf(fidb,[aa,'\n']);
aa = sprintf('i/jdm = %i %i',IDM, JDM);
fprintf(fidb,[aa,'\n']);
aa = ' tau_nwd: month,range = 01   0.0000000E+00   0.0000000E+00';
fprintf(fidb,[aa,'\n']);

fclose(fida);
fclose(fidb);






