% Extract T and Layer thicknes of Atl. Water
% in specified regions (Canada Basin, BG, ...)
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps/
startup

format long g
clear all
close

f_mat=1;

regn = 'ARCc0.08';
expt = 110;

pfld  = 'temp';
%f_extr = 1;  % =0 - load in extracted depth of Atl. Water, =1 -extract
s_fig  = 0;

sfig=0;

rg = 9806;
txb = 'plot_atlw_TdH.m';


pthfig  = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/%3.3i/fig_AtlLayer/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);

fmat = sprintf('%sarc08_110_atlw_TdH.mat',pthmat);
load(fmat);

TM = ATLW.dnmb;
DV = datevec(TM);
nrc= length(TM);
tx = [DV(1,1):1/12:DV(end,1)+0.99];

tmx = ATLW.Tmax;
zmx = ATLW.Z_Tmax;
z0c = ATLW.Z_T0;
dha = ATLW.dHatl;
tav = ATLW.Tatl_av;

yr1=DV(1,1);
yr2=DV(end,1);

figure(1); clf;
axes('Position',[0.08 0.6 0.85 0.35]);
plot(tx,tmx);
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+1],...
	'xtick',[yr1:yr2+1],...
	'xgrid','on',...
	'ygrid','on');
stl=sprintf('%s_%3.3i, Beauf.Sea, AtlL: Tmax, dgr',regn,expt);
title(stl,'Interpreter','none');

axes('Position',[0.08 0.1 0.85 0.35]);
plot(tx,zmx);
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+1],...
	'xtick',[yr1:yr2+1],...
	'xgrid','on',...
	'ygrid','on');
stl=sprintf('%s_%3.3i, Beauf.Sea, AtlL: Depth Tmax, m',regn,expt);
title(stl,'Interpreter','none');
bottom_text(txb,'pwd',1);

figure(2); clf;
axes('Position',[0.08 0.6 0.85 0.35]);
plot(tx,tav);
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+1],...
	'xtick',[yr1:yr2+1],...
	'xgrid','on',...
	'ygrid','on');
stl=sprintf('%s_%3.3i, Beauf.Sea, AtlL: Tav, dgr',regn,expt);
title(stl,'Interpreter','none');

axes('Position',[0.08 0.1 0.85 0.35]);
plot(tx,z0c);
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+1],...
	'xtick',[yr1:yr2+1],...
	'xgrid','on',...
	'ygrid','on');
stl=sprintf('%s_%3.3i, Beauf.Sea, AtlL: Depth T0C, m',regn,expt);
title(stl,'Interpreter','none');
bottom_text(txb,'pwd',1);

figure(3); clf;
axes('Position',[0.08 0.6 0.85 0.35]);
plot(tx,dha);
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+1],...
	'xtick',[yr1:yr2+1],...
	'xgrid','on',...
	'ygrid','on');
stl=sprintf('%s_%3.3i, Beauf.Sea, AtlL Thickness, m',regn,expt);
title(stl,'Interpreter','none');
bottom_text(txb,'pwd',1,'position',[0.08 0.45 0.6 0.1]);

    








