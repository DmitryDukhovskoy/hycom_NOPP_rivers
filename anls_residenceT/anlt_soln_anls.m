% Estimate lifetime of Greenland freshwater
% in the SPNA
% Analytical solution to the 
% ODE describing a
% dynamical system
% v'(y)+k*v(t)-F(t) = 0
% v - volume of Greenland FW in the Subpolar Gyre
% 1/k=tau - lifetime of Gr. FW in the Subpolar Gyre
% F(t) - gain, Gr. FWFlux (anomaly or total)
%
% To track GFWA in the Subpolar Gyre, consider GFWA flux
% only into the Subp. Gyre (exclluding Baffin Bay) - in fact this is SPNA region
% According to Dukhovskoy et al., 2016, ~85% of GFWA
% leaves Gr Shelf into the Labr Sea
% 15% goes to Baffin Bay and return to the Labr Sea with the Labrador Current
% only ~2% stays in Baffin Bay ~8% 
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

%
% First assumption that F(t)=const=F0
%F0=208.6; % mean GFWF anom, km3/yr
%F0=1130; % mean GFWF anom, km3/yr
Ftot0=208.6; % mean GFWF anom, km3/yr - total integrated over Greenland
% From Davis/Denmark Str Tracer flux analysis - see plot_greel_Tracr_POPgates.m
% and Greenland FW flux analysis - by segments
% Estimated mean annual GFWA influx into the SPG /SPNA
% is ~167 km3/yr 
%
fr=0.80; % fraction of total GFWA fluxed into SPG/SPNA
F0=fr*Ftot0; 
%t=[0:24];
t=[0:500];

Fgr2016=Ftot0*24; % cumulative GFW flux anomaly 1993-2016, km3
Fspg2016=fr*Fgr2016;  % estimated Fraction of the GFWA fluxed into the SPG 
Vgfwa_step_2016 = 2075; % ccumulated GFWA in SPNA in 2016 from HYCOM, constant Gr.discharge
Vgfwa_lnr_2016  = 2240; % accumulated GFWA in SPNA 2016 HYCOM, linearly incr.
                % estimated from tracers - see dS_FWC_timeseries_SubpolarGyre.m
V0=0; % initial volume of GFWA in the SPG
		
KK=1./[5:5:50];
KK2=1./[1:1:100];
V=[];
for ik=1:length(KK)
  k=KK(ik);
  V(ik,:)=V0*exp(-k*t)+F0/k*(1-exp(-k*t));
%  B=V0*exp(-k*t)+F0/k*(1-exp(-k*t));
  SS{ik}=sprintf('%i',1/k);
end
% To check V at fixed time for different k=1/time
V11=[];
for ik=1:length(KK2)
  k=KK2(ik);
  V11(ik,:)=V0*exp(-k*t)+F0/k*(1-exp(-k*t));
end

xl1 = 1990;
xl2 = 2050;
%YY=[1993:2017];
YY=t+1993;
iy  = max(find(YY<=xl2));

figure(1); clf;
axes('Position',[0.08 0.5 0.82 0.42]);
plot(YY,V);
hold on;
title('Analytical GFWA Volume in SPNA, km3, Const Gr Flux');
set(gca,'tickdir','out',...
 'xlim',[xl1 xl2],...      %	'xlim',[YY(1) YY(end)],...
	'xtick',[1990:5:YY(end)],...
	'ytick',[0:1000:20000],...
	'ylim',[-100 1.05*max(V(:,iy))],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
psc=legend(SS);
set(psc,'Position',[0.91,0.52 0.07 0.27]);

plot([xl1 xl2],[Fgr2016 Fgr2016],'k--','Linewidth',1.4); % GFWA vol in 2016
plot([xl1 xl2],[Vgfwa_step_2016 Vgfwa_step_2016],'k-','Linewidth',1.4); % GFWA vol from tracer


f_pcum=0;
if f_pcum==1
% Plot Cumulative Greenland FW flux anomaly
% constant rate = F0
sFgr=t*F0;
axes('Position',[0.07 0.1 0.82 0.27]);
hb=bar(YY,sFgr,0.8);
set(hb,'Facecolor',[0.7 0.7 0.7]);
hold on
plot([YY(1) YY(end)],[Fgr2016 Fgr2016],'k--','Linewidth',1.4);
set(gca,'tickdir','out',...
        'xlim',[YY(1) YY(end)],...
	'xtick',[1990:5:YY(end)],...
        'ylim',[0 1.05*max(max(sFgr))],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
stl=sprintf('Cumulative GFWA, km3, F0=const=%4.1f km3/yr',F0);
title(stl);
end

% Plot forcing function:
axes('Position',[0.07 0.1 0.82 0.27]);
hold on;
plot([1990 1993],[0 0],'k','linewidth',2);
plot([1993 1993],[0 F0],'k','linewidth',2);
plot([1993 xl2],[F0 F0],'k','linewidth',2);

set(gca,'tickdir','out',...
        'xlim',[xl1 xl2],...
        'xtick',[1990:5:YY(end)],...
        'ylim',[-1 1.1*F0],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
stl=sprintf('Forcing function, F0=const=%4.1f km3/yr',F0);
title(stl);


btx='anlt_soln_anls.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.03 0.4 0.05]);

%
% Find tau from given GFWA accumulated in the SPNA
% estimated from the HYCOM tracer
fprintf('================================\n');
fn = 1;
a0 = 0;
b0 = 200;
t0 = 24; % years since 1993, 2016 - when GFWA is estimated
tau0 = sub_RegulaFalsi(fn,a0,b0,t0);
fprintf('Steady-state solution, function %i, tau = %10.4f yrs\n\n',fn,tau0);

%
% Estimate tau - time scale for linear-increasing GFWA
fn = 2;
tau0 = sub_RegulaFalsi(fn,a0,b0,t0);
fprintf('Lin-incr solution, function %i, tau = %10.4f yrs\n',fn,tau0);
fprintf('================================\n\n');

% Find tau from given GFWA accumulated in the SPNA
% estimated from the HYCOM tracer - change some 
% parameters to estimate uncertainty, use step-function case
fprintf('================================\n');
fn = 3;
a0 = 0;
b0 = 200;
t0 = 24; % years since 1993, 2016 - when GFWA is estimated
tau0 = sub_RegulaFalsi(fn,a0,b0,t0);
fprintf('Steady-state solution, function %i, tau = %10.4f yrs\n\n',fn,tau0);

keyboard


% Analytical solution for 
% accelerating GFWA flow rate assuming
% linear acceleration, see fitFn_GreenlRunoff.m
%
% Linear regression to GFWA wrt to mean 1958-1992
% Regr coeff: intrcp 24.1 km3/yr, slope=15.2 km3/yr2
%
% F(t)=iF0+p*t, F0 - initial GFWA flux
% linear increase in GFWA ~15.2 km3/yr^2
% First assume F0=0
% p = 15.2
%iF0=0;    % 0 intercept

% Regression for total GFWA:
%iF0=24.0;  % intercept
%p=15.2;

% Regression for fraction of GFWA
% fluxed into SPG
% Make sure same values are in function2.m
iF0=17.44;
p=12.15;


V2=[];
clear aa SS2
for ik=1:length(KK)
  k=KK(ik);
  V2(ik,:)=iF0/k+p/(k^2)*(k*t-1)+(V0+p/(k^2)-iF0/k)*exp(-k*t);
  SS2{ik}=sprintf('%i',1/k);
  aa(ik,:)=p/k^2*exp(-k*t);  % exponent part
end

% To check V at fixed time for different k=1/time
V22=[];
for ik=1:length(KK2)
  k=KK2(ik);
  V22(ik,:)=iF0/k+p/(k^2)*(k*t-1)+(V0+p/(k^2)-iF0/k)*exp(-k*t);
end


% find when non-steady soln > stead soln:
for ik=1:length(KK2);
  ii=min(find(V22(ik,:)>=V11(ik,:) & V22(ik,:)>0));
  Tv2(ik)=1/KK2(ii);
end


% GFW flux is not constant - linear
figure(2); clf;
%axes('Position',[0.1 0.5 0.85 0.42]);
axes('Position',[0.08 0.5 0.82 0.42]);
plot(YY,V2);
hold on;

stl=sprintf('Analyt GFWA Vol in SPNA, km3, Linearly acceler. Gr Flux=%4.1f+%4.1f*t',iF0,p);
title(stl);
set(gca,'tickdir','out',...
	'xlim',[xl1 xl2],...
	'xtick',[1990:5:YY(end)],...
	'ylim',[-100 1.05*max(V2(:,iy))],...
	'ytick',[0:2000:50000],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
psc=legend(SS2);
%set(psc,'Position',[0.88,0.16 0.06 0.22]);
set(psc,'Position',[0.91,0.52 0.07 0.27]);

plot([xl1 xl2],[Fgr2016 Fgr2016],'k--','Linewidth',1.4); % GFWA vol in 2016
plot([xl1 xl2],[Vgfwa_lnr_2016 Vgfwa_lnr_2016],'k-','Linewidth',1.4); % GFWA vol from tracer


Fgr=iF0+t*p;
sFgr=cumsum(Fgr);
if f_pcum==1
axes('Position',[0.07 0.1 0.82 0.27]);
hb=bar(YY,sFgr,0.8);
set(hb,'Facecolor',[0.7 0.7 0.7]);
hold on
plot([YY(1) YY(end)],[Fgr2016 Fgr2016],'k--','Linewidth',1.4);
set(gca,'tickdir','out',...
        'xlim',[YY(1) YY(end)],...
	'xtick',[1990:5:YY(end)],...
        'ylim',[0 1.05*max(max(sFgr))],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
stl=sprintf('Cumulative GFWA, km3, %4.1f+%4.1f*t km3/yr',iF0,p);
title(stl);
end;

% Plot forcing function:
Ftrend=iF0+p*t;
axes('Position',[0.07 0.1 0.82 0.27]);
hold on;
plot([1990 1992.6],[0 0],'k','linewidth',2);
plot([1992.6 1993],[0 iF0],'k','linewidth',2);
plot(YY,Ftrend,'k','linewidth',2);

set(gca,'tickdir','out',...
        'xlim',[xl1 xl2],...
        'xtick',[1990:5:YY(end)],...
        'ylim',[-1 1.1*Ftrend(iy)],...
        'ytick',[0:100:1000],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
stl=sprintf('Forcing Function: %4.1f+%4.1f*t km3/yr',iF0,p);
title(stl);

%btx='anlt_soln.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.02 0.4 0.05]);


% Plot V for different k=1/t on year = 2016
TT=1./KK2;
f_plt3=0;
if f_plt3==1
figure(3); clf
axes('Position',[0.08 0.5 0.9 0.45]);
ip=24; % 2016
hold on;
plot(TT,V11(:,ip),'Color',[0 0.5 0.7],'Linewidth',2);
plot(TT,V22(:,ip),'Color',[0.7 0.3 0],'Linewidth',2);
plot([0 TT(end)],[Fgr2016 Fgr2016],'k--','Linewidth',1.4); % GFWA vol in 2016
plot([0 TT(end)],[Fspg2016 Fspg2016],'k--',...
     'Linewidth',1.4,'Color',[0 0.8 0.4]); % Fraction of GFWA in SPG in 2016
plot([0 TT(end)],[Vgfwa_lnr_2016 Vgfwa_lnr_2016],'k-','Linewidth',1.4); % GFWA vol from tracer
set(gca,'tickdir','out',...
        'xlim',[0 TT(end)],...
	'xtick',[0:5:TT(end)],...
        'ylim',[0 5200],...
	'ytick',[0:500:6000],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
stl=sprintf('GFWA (km^3) from 2 solutions at 2016, for different k=1/t');
title(stl);
xlabel('Time scales (yrs), t');

axes('Position',[0.08 0.28 0.3 0.1]); hold on;
plot([0.5 1],[0.5 0.5],'-','Color',[0 0.5 0.7],'Linewidth',2);
plot([0.5 1],[0.25 0.25],'-','Color',[0.7 0.3 0.],'Linewidth',2);
text(1.2,0.5,'F=const','Fontsize',14);
text(1.2,0.25,'F=F_0+p*t','Fontsize',14);
set(gca,'xtick',[],...
	'ytick',[],...
	'xlim',[0 2],...
	'ylim',[0.1 0.6],...
	'visible','off');

axes('Position',[0.45 0.28 0.3 0.1]); hold on;
plot([0.5 1],[0.75 0.75],'-','Color',[0 0 0],'Linewidth',1.4);
plot([0.5 1],[0.5 0.5],'--','Color',[0 0 0],'Linewidth',1.4);
plot([0.5 1],[0.25 0.25],'--','Color',[0 0.8 0.4],'Linewidth',1.4);
text(1.2,0.75,'GFWA in SPG in 2016','Fontsize',14);
text(1.2,0.5,'Cumulative GFWA flux 2016','Fontsize',14);
text(1.2,0.25,'Fraction of GFWA fluxed toSPG integr by 2016','Fontsize',14);
set(gca,'xtick',[],...
	'ytick',[],...
	'xlim',[0 2],...
	'ylim',[0.1 0.8],...
	'visible','off');



bottom_text(btx,'pwd',1,'Position',[0.08 0.2 0.4 0.05]);

end


% 
% Schwartz's paper on non-steady state residence time, 1973
% For the linearly increasing influx of GFWA:
% F(t)=iF0+p*t
% Checking the scales
% Turn-over time: 
% 1st definition - wrt to the influx
kk0 = 1/13; % fix 1/tau=k - estimate from the linearly increasing soln
ik=find(KK2==kk0); % pick V for k=1/14 yr-1
kk=KK2(ik);
Vk=V22(ik,:);
nt=length(t);
%Tto1=(iF0*t+p/2*t.^2)./(iF0+p*t);
Tto1=Vk./(iF0+p*t);

Vk(Vk==0)=1e-16;
% Mean age time for non-steady state case
% Mid-point quadrature
for it=1:nt
  Tt=t(it);
  a1=Tt-t(1:it);
  b1=iF0+p*t(1:it);
  c1=exp(-kk*(Tt-t(1:it)));
%  c1=exp(-1/12.8*(Tt-t(1:it)));
%
  Ttr(it)=sum(a1.*b1.*c1)/Vk(it);
end
%for it=1:nt
%  Tt=t(it);
%  a1=Tt-(t(1:it)-0.5);
%  b1=iF0+p*(t(1:it)-0.5);
%  c1=exp(-kk*(Tt-(t(1:it)-0.5)));
%
%  Ttr(it)=sum(a1.*b1.*c1)/Vk(it);
%end

figure(4); clf;
axes('position',[0.1 0.6 0.85 0.32]);
plot(t,Tto1,'Linewidth',3,'Color',[0 0 0]);
hold on;
plot(t,Ttr,'Linewidth',3,'Color',[0.6 0.6 0.6]);
lgd=legend('Tto(1)','Tma');
set(lgd,'Position',[0.08 0.4 0.16 0.09]);
set(gca,'tickdir','out',...
        'xtick',[0:50:500],...
        'xlim',[0 max(t)],...
        'ytick',[0:15],...
        'ylim',[0 ceil(max(Ttr))],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
xlabel('Years, t');
ylabel('Tau(t), years');
stl=sprintf('Turn-over T (1) and MeanAge T non-linear forcing, tau0=%i',round(1/kk0));
title(stl);
bottom_text(btx,'pwd',1,'Position',[0.08 0.3 0.4 0.05]);













