% Monthly mean T/S
% Plot only 1 layer at a time
% if need to layer-average - modify the code
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.04';
%expt = 011;  
expt = 012;  

pfld  = 'salin';
s_fig = 0;

iyr = 2005;
%mplt = [12,1,2]; % average over N months, if not, mplt=1month
mplt = [6,7,8]; % average over N months, if not, mplt=1month
nav = length(mplt);

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/data_mat/',expt);
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

fmat = sprintf('%sarc04_%3.3i_mnthly_%s_%i.mat',pthmat,expt,pfld,iyr);

fprintf('Loading %s\n',fmat);
load(fmat);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

SM={'J','F','M','A','M','J','J','A','S','O','N','D'};

smm=HH*0;
smo=[];
ccn=0;
for imm=1:nav
  im=mplt(imm);
  A = meanS(im).Favrg;
  if max(max(A))<0.1, 
    fprintf(' !!!!!  Month %i is missing\n',im);
    continue; 
  end;
  ccn=ccn+1;
  smm = smm+A;
  smo=[smo,SM{im}];
end
S = smm./ccn;

% Subp. N. Atl
xlim1 = 740;
xlim2 = 2500;
ylim1 = 300;
ylim2 = 2200;

% Greenland
%xlim1 = 900;
%xlim2 = 2200;
%ylim1 = 600;
%ylim2 = 2200;
%xlim1 = 
% SE Greenland:
%xlim1 = 1200;
%xlim2 = 1720;
%ylim1 = 700;
%ylim2 = 1280;

fprintf('Plotting ...\n');
nf = 1;
sb=33;
stl = sprintf('ARCc0.04-%3.3i, Mean S (cntr0=%3.1f) %i, %s',expt,sb,iyr,smo);
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
pfld='salin';
hps = [0.93 0.1 0.035 0.8];
sub_plot_scalar(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
		'c1',30,'c2',35,'cmp',2,'clbpos',hps);
contour(S,[30:0.5:36],'k');
contour(S,[sb sb],'k','Linewidth',1.8);

txtb = 'plot_mnthly_arc04_S.m';
bottom_text(txtb,'pwd',1);
