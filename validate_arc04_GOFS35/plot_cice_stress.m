% Plot surface stress emposed by sea ice on ocean
% 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.04';
expt = 22;
s_fig = 0;

pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/cice/';

iyr  = 2016;
imo  = 3;
mday = 27;
fout = sprintf('%s%3.3i_cice_inst.%i-%2.2i-%2.2i-00000.nc',pthbin,expt,iyr,imo,mday);
btx = 'plot_seaice_stress.m';

fprintf('Reading %s\n',fout);

% Plot atm ice stress:
ataux = squeeze(nc_varget(fout,'strairx'));
atauy = squeeze(nc_varget(fout,'strairy'));


figure(1); clf;
pcolor(ataux); shading flat;
caxis([-0.25 0.25]);
colorbar
stl = sprintf('CICE atm/ice stress, ataux, %i/%2.2i/%2.2i',iyr,imo,mday);
title(stl);
axis('equal');
set(gca,'xlim',[350 3200],...
        'ylim',[1000 5020]);
bottom_text(btx,'pwd',1);


figure(2); clf;
pcolor(atauy); shading flat;
colorbar
caxis([-0.25 0.25]);
stl = sprintf('CICE atm/ice stress, atauy, %i/%2.2i/%2.2i',iyr,imo,mday);
title(stl);
axis('equal');
set(gca,'xlim',[350 3200],...
        'ylim',[1000 5020]);
bottom_text(btx,'pwd',1);



