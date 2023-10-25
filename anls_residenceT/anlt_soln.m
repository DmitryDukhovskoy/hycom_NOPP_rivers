% Estimate lifetime of Greenland freshwater
% in the SPNA
% Analytical solution to the 
% ODE describing a
% dynamical system
% v'(y)+k*v(t)-F(t) = 0
% v - volume of Greenland FW in the Subpolar Gyre
% 1/k=tau - lifetime of Gr. FW in the Subpolar Gyre
% F(t) - gain, Gr. FWFlux (anomaly or total)
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

%
% First assumption that F(t)=const=F0
%F0=208.6; % mean GFWF anom, km3/yr
%F0=1130; % mean GFWF anom, km3/yr
F0=208.6; % mean GFWF anom, km3/yr
t=[0:24];
%t=[0:100];

KK=1./[2:2:20];
V=[];
for ik=1:length(KK)
  k=KK(ik);
  V(ik,:)=F0/k*(1-exp(-k*t));
  SS{ik}=sprintf('%i',1/k);
end

YY=[1993:2017];
%YY=t+1993;
figure(1); clf;
axes('Position',[0.1 0.5 0.85 0.42]);
plot(YY,V);
title('Analytical GFWA Volume in SPNA, km3, Const Gr Flux');
set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'ylim',[-100 1.05*max(max(V))],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
psc=legend(SS);
set(psc,'Position',[0.88,0.16 0.06 0.22]);
btx='anlt_soln.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.4 0.05]);


% Analytical solution for 
% accelerating GFWA flow rate assuming
% linear acceleration, see fitFn_GreenlRunoff.m
%
% Linear regression to GFWA wrt to mean 1958-1992
% Regr coeff: intrcp 24.1 km3/yr, slope=15.2 km3/yr2
%
% F(t)=F0+p*t, F0 - initial GFWA flux
% linear increase in GFWA ~15.2 km3/yr^2
% First assume F0=0
% p = 15.2
p=15.2;
V2=[];
for ik=1:length(KK)
  k=KK(ik);
  C = p/k^2;
  V2(ik,:)=p/k^2*((k*t-1)+exp(-k*t));
  SS2{ik}=sprintf('%i',1/k);
end

figure(2); clf;
axes('Position',[0.1 0.5 0.85 0.42]);
plot(YY,V2);
title('Analyt GFWA Volume in SPNA, km3, Linearly accelerating Gr Flux, F0=0');
set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'ylim',[-100 1.05*max(max(V2))],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
psc=legend(SS2);
set(psc,'Position',[0.88,0.16 0.06 0.22]);
btx='anlt_soln.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.4 0.05]);

% Same as above but F0=intercept
F0=24.1
p=15.2;
V3=[];
for ik=1:length(KK)
  k=KK(ik);
  V3(ik,:)=F0/k+p/k^2*(k*t-1)+(p/k^2-F0/k)*exp(-k*t);
%  SS2{ik}=sprintf('%i',1/k);
end

figure(3); clf;
axes('Position',[0.1 0.5 0.85 0.42]);
plot(YY,V3);
title(sprintf('GFWA Vol in SPNA, km3, Linearly increasing Gr Flux, F0=%4.1f',F0));
set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'ylim',[-100 1.05*max(max(V3))],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
psc=legend(SS2);
set(psc,'Position',[0.88,0.16 0.06 0.22]);
btx='anlt_soln.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.4 0.05]);

