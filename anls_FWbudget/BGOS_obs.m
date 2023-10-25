% Plot BGOS - FWContent of B. Gyre 
% from A. Proshutinsky data
% FWC in m, relative to Sref=34.8
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/110/fig_trac/');

arBG = 1.023e12; % area of BG region in OBGS, m2
YR=[2003:2017];
% FWC x1e3 km3
fwc=[16.9; 17.2; 18.2; 18.7; 19.8;,...
     21.8; 21.8; 21.9; 21.8; 22.1; ...
     20.6; 21.4; 22.4; 23.4; 23.5];

% Convert to FWC in m
fwc = fwc*1e12/arBG; % km3 -> m
sgm=[0.6; 0.8; 0.9; 0.8; 0.8;...
     0.6; 0.9; 0.8; 0.9; 0.7;...
     0.8; 0.7; 0.8; 0.8; 0.9];

figure(1); clf;
axes('Position',[0.09 0.65 0.85 0.25]);
bb=bar(YR,fwc);
hold
for ik=1:length(YR)
  yr=YR(ik);
  ds=sgm(ik);
  ff=fwc(ik);
  f1=ff-2*ds;
  f2=ff+2*ds;
  plot([yr yr],[f1 f2],'b-','linewidth',2);
  plot([yr-0.1 yr+0.1],[f1 f1],'b-','linewidth',2);
  plot([yr-0.1 yr+0.1],[f2 f2],'b-','linewidth',2);
end

set(bb,'FaceColor',[0.2 0.7 1]);
set(gca,'tickdir','out',...
	'xtick',[2003:2017],...
	'xlim',[2002.5 2017.5],...
        'xgrid','on','ygrid','on',...
	'Fontsize',16);
sll=sprintf('FWC (m) BGOS, Sref=34.8');
title(sll);
set(gcf,'Position',[1123 612 1359 673]);
btx = 'BGOS_obs.m';
bottom_text(btx,'pwd',1,'position',[0.08 0.4 0.8 0.1]);


% Plot anomaly
mfwc = fwc(1);
dfw = diff(fwc); % change, km3/yr
dfw = [0;dfw];

figure(2); clf;
axes('Position',[0.09 0.65 0.82 0.25]);
hold on
% FWC
%bb=bar(YR,fwc-mfwc);
dfwc=fwc-mfwc;
plot([1990 2017],[0 0],'k-','Color',[0.7 0.7 0.7]);
plot(YR+0.5,dfwc,'-','Color',[0 0.4 .8],'Linewidth',2);
plot(YR+0.5,dfwc,'.','Color',[0 0.6 1],'Markersize',18);
sll=sprintf('FWC anom  wrt to FWC(2003)=%3.1fxe3km3,Sref=34.8',mfwc);
title(sll);
% Uncertainty
%for ik=1:length(YR)
%  yr=YR(ik)+0.5;
%  ds=sgm(ik);
%  ff=fwc(ik)-mfwc;
%  f1=ff-2*ds;
%  f2=ff+2*ds;
%  plot([yr yr],[f1 f2],'-','Color',[0.4 0.4 0.4],'linewidth',2);
%  plot([yr-0.1 yr+0.1],[f1 f1],'-','Color',[0.4 0.4 0.4],'linewidth',2);
%  plot([yr-0.1 yr+0.1],[f2 f2],'-','Color',[0.4 0.4 0.4],'linewidth',2);
%end
set(gca,'tickdir','out',...
	'ylim',[-0.5 7],...
	'xtick',[1990:2017],...
	'xlim',[1990 2017]);

% FWC rate
axes('Position',[0.09 0.2 0.82 0.25]);
hold on
plot(YR+0.5,dfw,'-','Color',[0.7 0.3 0.7],'Linewidth',2);
plot(YR+0.5,dfw,'.','Color',[0.9 0.5 .9],'Markersize',18);
plot([1990 2017],[0 0],'k-','Color',[0.7 0.7 0.7]);
sll=sprintf('dFWC/dt x1e3 km3/yr, Sref=34.8');


%set(bb,'FaceColor',[0.2 0.7 1]);
set(gca,'tickdir','out',...
	'xtick',[1990:2017],...
	'xlim',[1990 2017]);
%        'xgrid','on','ygrid','on');

title(sll);

set(gcf,'Position',[603 31 1482 1311]);

btx = 'BGOS_obs.m';
bottom_text(btx,'pwd',1,'position',[0.01 0.05 0.4 0.05]);

if s_fig==1
  fgnm = sprintf('%sBGOS_obs_FWCanom',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end

