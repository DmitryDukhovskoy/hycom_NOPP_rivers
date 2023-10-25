% Add also HYCOM and RASM fluxes:
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

Sref=34.9;
YR1=2002;
YR2=2016;


pthdat='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_Myers/FRAM/';

btx = 'plot_FW_vol_fluxes_myers.m';

fout=sprintf('%smonthly_FW_vol_fluxes.mat',pthdat);
fprintf('Loading %s\n',fout);
load(fout);

FWF=FWV.FW_flx*1e-3;  % mSv
VF=FWV.Vol_flx*1e-6;  % Sv

% FWFlux:
[a1,a2]=size(FWF);
FWFm=mean(FWF);
prcL=prctile(FWF,25);
prcU=prctile(FWF,75);
dmm=reshape(FWF,a1*a2,1);
mnF=mean(dmm);
stdv=std(dmm);

% HYCOM FWF flux - saved in calc_FWFlux_straits.m
fh='/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/hycom_fwf.mat';
FH=load(fh);
fwfh=FH.FWFm*1e-3;
pLh=FH.prcL*1e-3;
pUh=FH.prcU*1e-3;
vx=[1:12,12:-1:1]+0.5;
vy=[pLh,fliplr(pUh)];
mnFh=FH.mnF;
sth=FH.stdv;

% RASM fluxes: 2002-2016, mSv
frasm='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_RASM/RASMfw.dat';
FR=load(frasm);
FR=FR';
Frm=mean(FR);
mnFr=mean(mean(FR));
stFr=std(reshape(FR,15*12,1));
pLr=prctile(FR,25);
pUr=prctile(FR,75);



% Vol Flux:
[a1,a2]=size(VF);
VFm=mean(VF);
pvL=prctile(VF,25);
pvU=prctile(VF,75);
dmm=reshape(VF,a1*a2,1);
mnV=mean(dmm);
stV=std(dmm);

% HYCOM vol fluxes:
% saved in plot_vol_transp.m
vh='/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/hycom_vol.mat';
VH=load(vh);
Vhm=VH.Fmn;
phL=VH.prcL;
phU=VH.prcU;
hx=[1:12,12:-1:1]+0.5;
hy=[phL,fliplr(phU)];
mnVh=VH.mnF;
stVh=VH.stdv;

% RASM fluxes: 2002-2016, Sv
vrasm='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_RASM/RASMvol.dat';
VR=load(vrasm);
VR=VR';
Vrm=mean(VR);
mnVr=mean(mean(VR));
stVr=std(reshape(VR,15*12,1));
prL=prctile(VR,25);
prU=prctile(VR,75);


mnth=[1:12]+0.5;
clr=[0.8 0.3 0.6];
clrh=[0 0.4 0.8];
clrr=[0 0.9 0.5];


figure(1); clf;
axes('Position',[0.1 0.55 0.83 0.38]);
%axes('Position',POS(2,:));
hold on
% HYCOM:
pp=patch(vx,vy,[0.9 0.95 1]);
set(pp,'edgecolor','none');
plot(mnth,fwfh,'Linewidth',2.5,'Color',clrh); % mSv, + to the Arctic

plot(mnth,FWFm,'Linewidth',2.5,'Color',clr); % mSv, + to the Arctic
for ik=1:12
  plot([ik+0.5 ik+0.5],[prcL(ik) prcU(ik)],'--','Color',clr);
end

plot(mnth,Frm,'Linewidth',2.5,'Color',clrr);
for ik=1:12
  plot([ik+0.5 ik+0.5],[pLr(ik) pUr(ik)],'--','Color',clrr);
end


stxt=sprintf('NEMO: %4.1f+/-%4.1f',mnF,stdv);
stht=sprintf('HYCOM: %4.1f+/-%4.1f',mnFh,sth);
strt=sprintf('RASM: %4.1f+/-%4.1f',mnFr,stFr);
%text(1,0.9*min(prcL*1e-3),stxt,'Fontsize',16);
text(5,-140,stht,'Fontsize',15,'Color',clrh);
text(5,-125,stxt,'Fontsize',15,'Color',clr);
text(5,-155,strt,'Fontsize',15,'Color',clrr);

stl = sprintf('Fram FWFlux (+ North), HYCOM & NEMO Sref=%3.1f',Sref);
title(stl);
set(gca,'tickdir','out',...
	'xlim',[0.9 12.9],...
	'xtick',[1:12],...
	'ylim',[-165 -20],...
	'ytick',[-200:20:0],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);


bottom_text(btx,'pwd',1,'pos',[0.04 0.4 0.4 0.04]);


% 
% Volume transports:
figure(2); clf;
axes('Position',[0.1 0.55 0.83 0.38]);
hold on;

pp=patch(hx,hy,[0.9 0.95 1]);
set(pp,'edgecolor','none');
plot(mnth,Vhm,'Linewidth',2.5,'Color',clrh); % mSv, + to the Arctic

plot(mnth,VFm,'Linewidth',2.5,'Color',clr);
for ik=1:12
  plot([ik+0.5 ik+0.5],[pvL(ik) pvU(ik)],'k--','Color',clr);
end

plot(mnth,Vrm,'Linewidth',2.5,'Color',clrr);
for ik=1:12
  plot([ik+0.5 ik+0.5],[prL(ik) prU(ik)],'k--','Color',clrr);
end


stxt=sprintf('NEMO: %4.1f+/-%4.1f',mnV,stV);
stht=sprintf('HYCOM: %4.1f+/-%4.1f',mnVh,stVh);
strt=sprintf('RASM: %4.1f+/-%4.1f',mnVr,stVr);
%text(1,0.9*min(prcL*1e-3),stxt,'Fontsize',16);
text(5,-4.2,strt,'Fontsize',15,'Color',clrr);
text(5,-4.8,stht,'Fontsize',15,'Color',clrh);
text(5,-5.4,stxt,'Fontsize',15,'Color',clr);


stl = sprintf('Fram VolFlux (+ North), HYCOM (93-16) & NEMO, 2002-2016');
title(stl);
set(gca,'tickdir','out',...
	'xlim',[0.9 12.9],...
	'xtick',[1:12],...
	'ylim',[-6 0.1],...
	'ytick',[-10:0.5:0],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);

bottom_text(btx,'pwd',1,'pos',[0.04 0.4 0.4 0.04]);



