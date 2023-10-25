% Plot XY view of 
% final archive fields (nest fields)
% created from GLBb0.08 T07 32 layers ->  to ARCc0.08 on T11 41 layers
% or any other archv file
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.08';
f_fort = 1; % =0 - check interp. fields created remap32to41lrs_fixedZ.m
            % =1 - created in remap_32to41FORT/remap32to41lrs.F90 

lplt = 1; % layer to plot
%fld  = 'v-vel.'; 
fld  = 'u-vel.'; 
%fld  = 'v_btrop'; 

switch(fld),
 case('salin')
  c1 = 32; 
  c2 = 34;
 case ({'u-vel.','v-vel.'})
  c1 = -0.4;
  c2 = 0.4;
 case ({'u_btrop','v_btrop'})
  c1 = -0.2;
  c2 = 0.2;
end


pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
pthout  = '/Net/kronos/ddmitry/hycom/ARCc0.08/archv_41layers/';
pthrlx = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/relax/110/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/topo_grid/',R);

hg=2^100;
rg=9806;  % convert pressure to depth, m

%yr   = 1993;
yr   = 2005; % nest files corrected on Gordon
iday = 248;
hr   = 0;
TV   = 11; % arc08
%TV   = 17; % arc04


% Topo:
fltopo=sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV);
%fltopo=sprintf('%sdepth_ARCc0.04_%2.2iDD.nc',pthtopo,TV);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn] = size(HH);
%HH(HH>0)=nan;

%  files:
% created in matlab - there can be error - U,V not rotated from GLBb
finaNEW = sprintf('%sarchv41lr_arcT%2.2i.%4.4i_%3.3i_%2.2i.a',...
		 pthout,TV,yr,iday,hr);
finbNEW = sprintf('%sarchv41lr_arcT%2.2i.%4.4i_%3.3i_%2.2i.b',...
		 pthout,TV,yr,iday,hr);

if f_fort==1
  fnm = sprintf('archv_arcT%2.2iL41.%4.4i_%3.3i_00',TV,yr,iday);
% fnm = sprintf('archv_arcT%2.2iL41.%4.4i_%3.3i_00_FORTerr',TV,yr,iday);
% fnm = sprintf('archv_arcT%2.2iL41.%4.4i_%3.3i_00_FORT',TV,yr,iday); %corrcted right way
% fnm = sprintf('archv_arcT%2.2iL41.%4.4i_%3.3i_00_crct',TV,yr,iday); % corrcted fast way
 finaNEW = sprintf('%s%s.a',pthout,fnm); 
 finbNEW = sprintf('%s%s.b',pthout,fnm); 
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

[F,n,m,l] = read_hycom(finaNEW,finbNEW,fld,'r_layer',lplt);
F(F>1e20)=nan;
A=squeeze(F);


figure(2); clf;
pcolor(A); shading flat;
hold on;
contour(HH,[0 0],'k');
caxis([c1 c2]);
colorbar;
stt{1}=sprintf('%s, Layer: %i, ARCc0.08T11L41, from GLBb0.08T07L32',fld,lplt);
stt{2}=sprintf('%s',finaNEW);

axis('equal');
set(gca,'xlim',[1 IDM],'ylim',[1 JDM]);
title(stt,'Fontsize',14,'Interpreter','none');

  
  
  


