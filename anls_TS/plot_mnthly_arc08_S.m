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

regn = 'ARCc0.08';
%expt = 110;  
expt = 112;  

pfld  = 'salin';
s_fig = 0;

iyr1 = 2017;
iyr2 = 2017;
%iyr1 = 2012;
%iyr2 = 2014;
%mplt = [12,1,2]; % average over N months, if not, mplt=1month
mplt = [6,7,8]; % average over N months, if not, mplt=1month
%mplt = [6,7,8,9]; % average over N months, if not, mplt=1month
nav = length(mplt);

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

SM={'J','F','M','A','M','J','J','A','S','O','N','D'};

smm=HH*0;
smo=[];

ncc=0;
for iyr=iyr1:iyr2
  fmat = sprintf('%sarc08_%3.3i_mnthly_%s_%i.mat',pthmat,expt,pfld,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);

  for imm=1:nav
    im=mplt(imm);
    A = meanS(im).Favrg;
    smm = smm+A;
    if iyr==iyr1
      smo=[smo,SM{im}];
    end
    ncc=ncc+1;
  end

end
S = smm./ncc;

% subpolar NA
xlim1 = 370;
xlim2 = 1250;
ylim1 = 150;
ylim2 = 1100;
% Greenland
%xlim1 = 450;
%xlim2 = 1100;
%ylim1 = 300;
%ylim2 = 1100;
% SE Greenland:
%xlim1 = 600;
%xlim2 = 860;
%ylim1 = 350;
%ylim2 = 640;

nf = 1;
sb=33;
stl = sprintf('ARCc0.08-%3.3i, Mean Surf S (cntr0=%3.1f) %i-%i, %s',...
	      expt,sb,iyr1,iyr2,smo);
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
pfld='salin';
%hps = [0.9 0.1 0.035 0.8];
hps = [0.93 0.1 0.02 0.8];
Fpos = [1350  255  1034  1085]; % Figure gcf position
sub_plot_scalar(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
		'c1',30,'c2',35,'cmp',2,'clbpos',hps,'figpos',Fpos);
contour(S,[30:0.5:36],'k');
contour(S,[33 33],'k','Linewidth',1.6);


txtb = 'plot_mnthly_arc08_S.m';
bottom_text(txtb,'pwd',1);
