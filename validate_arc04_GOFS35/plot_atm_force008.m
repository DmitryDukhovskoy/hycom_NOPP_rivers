% Atmospheric forcing *.[ab] files
% already interpolated onto HYCOM grid
% from *d or nc atm files
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

regn = 'ARCc0.08';
expt = 122;
s_fig = 0;

f_plt = 0; % plot instant 1-hr fields
f_avrg = 1; % plot average fields

%fld = 'glbrad';
%fld = 'lwdflx';
fld = 'airtmp';
hYR = 117;
AB  = 'f';
% Average over nrec
% assuming 24 records per day
dS=7;   % day to start
hrS=12;  % hr to start
dt = 1; % 1hr data
irc1 = (dS-1)*24/dt+hrS+1;
nrec = 1; % # of recrods to read/average


btx = 'plot_cice_forcing008.m';

%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

ID = nn;
JD = mm;
IJDM=ID*JD;
lrec = ID*JD*8;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

%pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/incoming/';
pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/atm_force008/';

asum = zeros(JD,ID);

fina = sprintf('%s%s_%3.3i%s.a',pthbin,fld,hYR,AB);
finb = sprintf('%s%s_%3.3i%s.b',pthbin,fld,hYR,AB);
fid1 = fopen(fina,'r');  %

% Skip to start record:
recS = (dS-1)*24/dt+hrS+1; % record to start, hr=0, day=dS
stat=fseek(fid1,(recS-1)*(IJDM+npad)*4,-1);
fprintf('Skip %i records \n',recS-1);
ccL=0;
for irc=1:nrec
  fprintf('Reading record %i\n',irc);
  dmm=fread(fid1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid1,npad,'float32','ieee-be');  % read npad 
%  fprintf('ii=%i, dmm=%ix%i\n',ii,size(dmm));
%keyboard
  dmm=reshape(dmm,ID,JD);
  ccL=ccL+1;
  A=dmm';
  asum = asum+A;
end

AA = asum./ccL;


switch(fld)
 case('lwdflx');
  c1=220;
  c2=360;
 case('glbrad');
  c1=100;
  c2=400;
 case('airtmp');
  c1=-10;
  c2=20;
end

CMP = create_colormap6(400,c1,c2);
cmp = CMP.colormap;

figure(10); clf;
set(gcf,'Position',[1750 519 777 810]);
pcolor(AA); shading flat;
colormap(cmp);
caxis([c1 c2]);

dt=1;
nhrs = ccL*dt;
maxAav = max(max(AA));

SM ='abcdefghijkl';
yr = hYR+1900;
mo = strfind(SM,AB);

dn1=datenum(yr,mo,dS,hrS,0,0);
dv1=datevec(dn1);
dn2=dn1+(ccL-1)*dt/24;
dv2=datevec(dn2);

stl = sprintf('%s_%i%s, Avrg: %i/%2.2i/%2.2i:%2.2i-%2.2i/%2.2i:%2.2i, nrec=%i, max=%5.4g W/m2',...
														fld,hYR,AB,dv1(1:4),dv2(2:4),ccL,maxAav);



hold on;
contour(HH,[0 0],'k');
axis('equal');
title(stl,'Interpreter','none');

set(gca,'xlim',[50 1575],...
								'ylim',[250 2000],...
								'fontsize',12);

clb = colorbar;
set(clb,'Position',[0.89 0.11 0.025 0.81],...
								'Fontsize',12);

btx = 'plot_atm_force008.m';
bottom_text(btx,'pwd',1);




