% Plot diffusion fields from global hycom GOFS3.1
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

%pfld = 'thkdf4';
%pfld = 'veldf4';
pfld = 'veldf2';
expt = 'GLBb0.08_GOFS3.1_expt73.7';
pthdf = '/nexsan/people/ddmitry/hycom/GLBb0.08_GOFS3.1_expt73.7/';

idm = 4500;
jdm = 3298;
kdm = 41;
ijdm = idm*jdm;
npad = 4096-mod(ijdm,4096);


fgrid = sprintf('%sregional.grid',pthdf);
ftopo = sprintf('%sdepth_GLBb0.08_09m11.a',pthdf);
GRD = read_grid_bath(fgrid,ftopo);
LAT = GRD.PLAT;
LON = GRD.PLON;
HH  = GRD.Topo;
[mm,nn]= size(HH);
I=find(HH>1.e20);
HH=-HH;
HH(I)=100;

fina = sprintf('%s%s.a',pthdf,pfld);
finb = sprintf('%s%s.b',pthdf,pfld);

fprintf('Opening %s\n',fina);
fprintf('Opening %s\n',finb);
fida = fopen(fina,'r','ieee-be');

A=fread(fida,ijdm,'float32'); % read 2D field
%dm1=fread(fida,npad,'float32');  % Padding = size(toto)
minh=min(A);
maxh=max(A);

fprintf('Reading %s, minh/maxh: %9.4f %9.4f\n',fina,minh,maxh);

AA=reshape(A,idm,jdm)';
AA(HH>=0)=nan;

figure(10); clf;
set(gcf,'Position',[1670         585         865         744]);
axes('Position',[0.1 0.1 0.8 0.8]);
hold on;
pcolor(AA); shading flat
axis('equal');
set(gca,'Color',[0.8 0.8 0.8],...
        'xlim',[0 nn],...
        'ylim',[0 mm]);


btx = 'plot_GLB_gofs31_diffusion.m';
title(sprintf('%s %s, min=%7.5f max=%7.5f',expt,pfld,minh,maxh),'interpreter','none');
caxis([0 maxh]);
%contour(HH,[0 0],'k');
bottom_text(btx,'pwd',1);
hb = colorbar;




