% Check created archv files
% Nest files
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

%R = 'ARCc0.08';
R = 'ARCc0.04';
f_fort = 1; % =0 - check interp. fields created remap32to41lrs_fixedZ.m
            % =1 - created in remap_32to41FORT/remap32to41lrs.F90 

pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthout  = '/Net/kronos/ddmitry/hycom/ARCc0.08/archv_41layers/';
pthout  = '/Net/mars/ddmitry/hycom/ARCc0.04/nest_files/';
%pthrlx = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/110/';
%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);

hg=2^100;
rg=9806;  % convert pressure to depth, m

yr   = 2004;
iday = 362;
hr   = 0;
switch(R)
 case('ARCc0.08'),
  TV  = 11; % arc08
  TVn = '11';
  flnmbs = 'archv_arcT11L41';
 case('ARCc0.04'),
  TV  = 17; % arc04
  TVn = '17DD';
%  flnmbs = 'archv_arc04T17L41'
  flnmbs = 'archv_arc04T17L41'
end

% Topo:
fltopo=sprintf('%sdepth_%s_%s.nc',pthtopo,R,TVn);
%fltopo=sprintf('%sdepth_ARCc0.04_%2.2iDD.nc',pthtopo,TV);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn] = size(HH);
%HH(HH>0)=nan;

%  files:
finaNEW = sprintf('%sarchv41lr_arcT%2.2i.%4.4i_%3.3i_%2.2i.a',...
		 pthout,TV,yr,iday,hr);
finbNEW = sprintf('%sarchv41lr_arcT%2.2i.%4.4i_%3.3i_%2.2i.b',...
		 pthout,TV,yr,iday,hr);

if f_fort==1
%  finaNEW = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00_FORTerr.a',...
%		 pthout,TV,yr,iday);
%  finbNEW = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00_FORTerr.b',...
%		 pthout,TV,yr,iday);
%  finaNEW = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00.a',...
%		 pthout,TV,yr,iday);
%  finbNEW = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00.b',...
%		 pthout,TV,yr,iday);
% ARCc0.04 nest files created on Gordon:
  finaNEW = sprintf('%sarchv_arc04T%2.2iL41.%4.4i_%3.3i_00.a',...
		 pthout,TV,yr,iday);
  finbNEW = sprintf('%sarchv_arc04T%2.2iL41.%4.4i_%3.3i_00.b',...
		 pthout,TV,yr,iday);
%  finaNEW = sprintf('%s%s.%4.4i_%3.3i_00.a',pthout,flnmbs,yr,iday);
%  finbNEW = sprintf('%s%s.%4.4i_%3.3i_00.b',pthout,flnmbs,yr,iday);
%  finaNEW = sprintf('%s%s.%4.4i_%3.3ifxd.a',pthout,flnmbs,yr,iday);
%  finbNEW = sprintf('%s%s.%4.4i_%3.3ifxd.b',pthout,flnmbs,yr,iday);
end

% Open files for reading
disp('Openning files ...');
fida_NEW = fopen(finaNEW,'r');
fidb_NEW = fopen(finbNEW,'rt');

for nl=1:10
  aa=fgetl(fidb_NEW);

  if nl==8
    dmm=aa(2:8);
    IDM=str2num(dmm);
  end
  if nl==9
    dmm=aa(2:8);
    JDM=str2num(dmm);
  end

  if exist('JD','var')
  if (ID~=nn | JD~=mm)
    fprintf('ERR: Check dimensions in new grid  and old grid:\n');
    fprintf('New: IDM=%i, JDM=%i; Old: IDM=%i, JDM=%i\n',nn,mm,IDM,JDM);
    error('*** STOPPING ***');
  end
  end

end

IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

% Start reading fields:
while ~isempty(aa),
  aa=fgetl(fidb_NEW);
  if isempty(aa) | aa==-1; fprintf('--- EoF -----\n'); break; end;
%  if ~isstring(aa); fprintf('--- EoF -----\n'); break; end;
  Is=strfind(aa,'=');
  fld=aa(1:Is-2);
  S = sscanf(aa(Is+1:end),'%g');
  smin=S(5);
  smax=S(6);
  Lr=S(3);
  
  [A,counta]=fread(fida_NEW,IJDM,'float32','ieee-be');
  toto = fread(fida_NEW,npad,'float32','ieee-be');
  I=find(A<0.01*hg);
  amin=min(A(I));
  amax=max(A(I));
  i99=find(A < -999.89999 & A > -999.900001);
  
  fprintf('\n %s  Layr %i, .a, .b min = %16.7e %16.7e\n',fld,Lr,amin,smin);
  fprintf(' %s  Layer %i, .a, .b max = %16.7e %16.7e\n',fld,Lr,amax,smax);
  
  dmn=abs(1-amin/smin);
  dmx=abs(1-amax/smax);
  
  if dmn>1e-4
    fprintf('.a and .b not consistent dmn=%9.7f\n',dmn);
    keyboard;
  end
  if dmx>0.0001
    fprintf('.a and .b not consistent dmx=%9.7f\n',dmx);
    keyboard;
  end
  if ~isempty(i99)
    fprintf('Found -999.9 in the field %s\n',fld);
    keyboard
  end

  f_plt=0;
%  if strncmp(fld,'u-vel',5),
  if strncmp(fld,'v-vel',5),
%  if strncmp(fld,'temp',4),
%   if strncmp(fld,'u_btrop',5);
    f_plt = 1;
%    keyboard;
  end
  if f_plt>0
    fprintf('Plotting %s\n',fld);
    B=reshape(A,IDM,JDM)';
    B(B>0.1*hg)=nan;
    [jj,ii]=find(B<=-990);
    figure(1); clf;
    pcolor(B); shading flat;
    hold on;
    contour(HH,[0 0],'k');
    plot(ii,jj,'ro');
    stt=sprintf('%s %s',R,fld);
    title(stt,'interpreter','none');
    keyboard
%  [jj,ii]=find(B>0);
%
%  Iland=find(B>0 & HH>=0);
  end
  
%  if strncmp(fld,'u-vel',5)
%    keyboard
%  end
  
end

fprintf(' No errors found \n');


  
  
  


