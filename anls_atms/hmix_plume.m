% Calculate depth of mixing of a plume
% driven by winds
% based on
% "Impact of Offshore Winds on a Buoyant River Plume System"
% JOSEPH T. JURISA* AND ROBERT J. CHANT
% JPO, 2013
% Eq. (2)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f=2*7.29e-5*sind(65);
spd=[15:0.1:45];  % downwelling wind speed
Cd = 1.e-3*(2.7./spd+0.142+0.0764.*spd);
Cd(spd<1.e-9)=1.e-12;
%! Capping the drag for strong winds:              
Cd(Cd>1.4e-3); Cd=1.4e-3; 
%tauw=1.2*Cd*spd.^2; % wind stress

rhoa=sw_dens0(31.5, 0);
rhop=sw_dens0(31.0,0); % plume density
drho=rhoa-rhop;
drho=[0.01:0.01:0.3];;

[S,DR]=meshgrid(spd,drho);

tauw=1.2*Cd*S.^2; % wind stress
Ric=3.0;
Ue=tauw./(rhop*f);
g=9.8;
hp=(4*Ric*rhoa*Ue./(g*DR)).^(1/3);

cmp=colormap('parula');
figure(1); clf;
pcolor(S,DR,hp); shading interp;
caxis([30 100]);
hold on;
hb=colorbar;
set(hb,'Fontsize',18);
contour(S,DR,hp,[30:10:120],'k-');
colormap(cmp);
set(gca,'tickdir','out',...
	'xtick',[15:5:45],...
	'xlim',[15 45],...
	'ytick',[0:0.05:0.3],...
	'ylim',[0.01 0.3],...
	'Fontsize',18);
xlabel('Downwelling Wind Speed, m/s');
ylabel('Delta Rho, kg/m3');
title('Greenl. Plume Wind-Driven Mixing Depth (m), cntr. [30:10:120]');
btx='hmix_plume.m';
bottom_text(btx,'pwd',1);



