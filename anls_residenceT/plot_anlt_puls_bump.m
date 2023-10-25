% Plot solutions for pulse  and bymp forcing functions
%
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

Ftot0=208.6; % mean GFWF anom, km3/yr - total integrated over Greenland
% From Davis/Denmark Str Tracer flux analysis - see plot_greel_Tracr_POPgates.m
% and Greenland FW flux analysis - by segments
% Estimated mean annual GFWA influx into the SPG /SPNA
% is ~138 km3/yr 
%
%fr=0.80; % fraction of total GFWA fluxed into SPG/SPNA
fr=1.918;     % make bump bigger to have higher FW content accumuated in the region
%F0=fr*Ftot0;  % Hakkinen 1993 estimated FW increase ~400-500km3/yr through Fram Strait during 1963-68
F0=2000;       % Curry and Muaritzen, 2005: GSA - 2000 km3/yr of FW anomaly flux
%t=[0:24];
t=[0:100];

V0=0; % initial volume of GFWA in the SPG

KK=1./[5:5:50];
V=[];
for ik=1:length(KK)
  k=KK(ik);
  SS{ik}=sprintf('%i',1/k);
end



xl1 = 1990;
xl2 = 2050;
%YY=[1993:2017];
YY=t+1993;
iy  = max(find(YY<=xl2));



% Pulse function: F=delta(t-a)*F0;
a=10;
b=15;
KK=1./[5:5:50];
ita=max(find(t<a));
Vp=[];
for ik=1:length(KK)
  k=KK(ik);
  v1=F0*exp(-k*(t-a));
  v1(1:ita)=0;     % step function, t<a = 0
  Vp(ik,:)=v1+V0*exp(-k*t);
end


% Bump is addition of 2 step functions u(t-a) and u(t-b)
% u(t-a) = 0 fot t<0 and =1 for t>=a, 
% Bump B = u(t-a) - u(t-b)
Vb = [];
itb = max(find(t<b));
for ik=1:length(KK)
  k=KK(ik);
  v1=F0/k*(1-exp(-k*(t-a)));
  v1(1:ita)=0;     % step function, t<a = 0
  v2=F0/k*(1-exp(-k*(t-b)));
  v2(1:itb)=0;
  Vb(ik,:)=v1-v2;
end


figure(1); clf;
set(gcf,'Position',[1190 608 1336 726]);
axes('Position',[0.08 0.5 0.82 0.42]);
plot(YY,Vb);
hold on;
title('Analytical GFWA Volume in SPNA, km3, Forcing Bump Fn');
set(gca,'tickdir','out',...
 'xlim',[xl1 xl2],...      % 'xlim',[YY(1) YY(end)],...
 'xtick',[1990:5:YY(end)],...
 'ytick',[0:500:20000],...
 'ylim',[-100 1.05*max(max(Vb(:,1:iy)))],...
 'xgrid','on',...
 'ygrid','on',...
 'Fontsize',14);
psc=legend(SS);
set(psc,'Position',[0.91,0.52 0.07 0.27]);

%plot([xl1 xl2],[Fgr2016 Fgr2016],'k--','Linewidth',1.4); % GFWA vol in 2016
%plot([xl1 xl2],[Fspna2016 Fspna2016],'k-','Linewidth',1.4); % GFWA vol from tracer

% Plot forcing function:
axes('Position',[0.08 0.1 0.82 0.27]);
hold on;
iYs=YY(1)+a;
iYe=YY(1)+b;
plot([1990 iYs],[0 0],'k','linewidth',2);
plot([iYs iYs],[0 F0],'k','linewidth',2);
plot([iYs iYe],[F0 F0],'k','linewidth',2);
plot([iYe iYe],[0 F0],'k','linewidth',2);
plot([iYe xl2],[0 0],'k','linewidth',2);

set(gca,'tickdir','out',...
        'xlim',[xl1 xl2],...
        'xtick',[1990:5:YY(end)],...
        'ylim',[-1 1.1*F0],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
stl=sprintf('Forcing function, bump F0=%4.1f km3/yr',F0);
title(stl);


btx='anlt_soln_anls.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.03 0.4 0.05]);












