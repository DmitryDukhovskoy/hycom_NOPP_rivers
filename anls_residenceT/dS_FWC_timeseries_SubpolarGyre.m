% Time series of FWC and dltS 
% in Subpolar Gyre region
%
% Uses corrected code extracting tracer
% that integrates tracer over the specified layers
% see: extr_MassTrcr_month
% Estimate S change in the upper layer due to surplus FW flux
% from Greenland 
% Note estimate of Greenland runoff fraction
% for specified boxes may use wrong
% code - tracer concentration averaged over layers
% should be integrated giving a total mass within
% a layer

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

%fmap='SPG_noBaf'; % Subpolar with North Sea and no Baffin
fmap='SPG_noNorth'; % SPG no North Sea and no Baffin

s_fig  = 0;
f_extr = 1; % rederive FWC within specified domain

regn = 'ARCc0.08';
expt = 112;  % 110 - no Gr runoff  or 112 - with Gr runoff
YR1 = 1993;
if expt==112
  YR1=1995;
end

f_plat = 0;   % add Platov's time series
f_myr  = 0;   % add Myer's time series 
nTr = 1; 
IntSrf = 1; % = 1 - integrates over all layers from 0 - ilv, for ilv=1 is the same
%ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % 150-300m
%ilv = 4; % 300-500 m
ilv = 5; % whole depth 

% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS); 


zz1 = LRS(ilv,1);
zz2 = LRS(ilv,2);

if IntSrf==1; zz1=0; end;

dz=abs(zz2-zz1);


fprintf('FW Volume SPG, expt=%i, ilv=%i, %i - %i, nTr=%i, %i-2016\n',expt,ilv,zz1,zz2,nTr,YR1);


pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthsav  = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt);
pthmatLP= '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
pthmat2 = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/data_mat/';


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2
[XX,YY] = meshgrid((1:nn),(1:mm));

hmsk=HH;
hmsk(HH<0)=nan;

%
hZ = abs(LRS(ilv,2)-LRS(ilv,1));

if IntSrf==1
  hZ = dz;
end


%dSregions_map_v2


% Subpolar Gyre:
% Labr - Irm - eastern & central N. Atlantic
switch(fmap) 
 case('SPG_noBaf') 

IGR = [  430         729
         461         671
         559         660
         774         676
         781         608
         839         533
         855         530
         930         501
        1031         445
        1099         410
        1192         403
        1228         321
        1218         185
        1080         145
         444         145
         416         246
         367         492
         374         591
         379         729];

 case('SPG_noNorth');

IGR = [  430         729
         461         671
         559         660
         774         676
         781         608
         839         533
         855         530
         930         501
        1031         445
        1099         410
        1065         350
        1087         288
        1106         204
        1118         148
        1080         145
         444         145
         409         188
         424         218
         416         246
         367         303
         367         492
         374         591
         367         729];

end

fspg='SPG_noNorth_indx.mat';
save(fspg,'IGR');



INP = inpolygon(XX,YY,IGR(:,1),IGR(:,2));
INGr= find(INP==1 & HH<0);
II = find(HH<0);



frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat);

f_griv = 0; %=1 - rederive cumulative Greenland FWFlux from Bamber
if f_griv==1
  [TMg,Fgr]=sub_read_Greenland_v3; % Greenland total FWF, km3/mo
  % Calculate annual anomalies
  % and cumulative up to date
  DVg = datevec(TMg);
  ngr = length(TMg);
  nyr =ngr/12;
  dmm = reshape(Fgr,[12,nyr]);

  Fyr = sum(dmm); % km3/yr
  Ygr = [DVg(1,1):DVg(end,1)];
  ii=find(Ygr==1990);
  Fmn = mean(Fyr(1:ii));
%  ism = find(Ygr==dv0(1,1));
  cFWF = cumsum(Fyr-Fmn);
%  fwf0 = cFWF(ism); % km3


% Greenland runoff with step-function increase of surplus runoff
% increase = mean (1993-2016) = 209 km3/yr
  Fmn = 818.3;
  dFmn = 209;
  Fstep = Fyr;
  iyr = find(Ygr==1993);
  Fstep(iyr:end) = Fmn+dFmn; % Greenland runoff with step-function increase rate start on 1993
  cFWF_step = cumsum(Fstep-Fmn);


  save(frv,'cFWF','Fyr','Ygr','Fstep','cFWF_step');
else
  fprintf('f_griv %i, Loading %s\n',f_griv,frv);
  load(frv);
%  ii=find(Ygr==1990);
end  


%

% Tracer fraction in grid cells
fextr=sprintf('%sFWC_TimeSer_%s_ilv%i.mat',pthmat,fmap,ilv);
if f_extr==1
  nTr = 1;
  iyy=0;
  for iyr=YR1:2016
    iyy=iyy+1;
    for im=1:12
      dnmb = datenum(iyr,im,15);
      dv0  = datevec(dnmb);
      yr   = dv0(1);
      iday = dnmb-datenum(yr,1,1)+1;

      ism = find(Ygr==dv0(1,1));
      fwf0 = cFWF(ism); % km3
%
% FWFlux from step-function increase (constant Greenland anomaly after 1993)
      fwf0_stp = cFWF_step(ism);  
      rr = sub_fraction_tracerMass(expt,regn,Acell,HH,...
                                   dnmb,nTr,ilv,hZ,IntSrf);

  % Estimate volume of Greenland surplus FW in grid cell
  % = ratio of FW tracer*total FWF (integrated over time 1990-date requested)
      Vfw = fwf0*rr*1e9; % m3 Gr FWF anomaly integrated in time
      VfwS = fwf0_stp*rr*1e9; % Gr FWF anomaly from step-function increase
  % Integrate vol of FW in the region
      vfw=nansum(Vfw(INGr))*1e-9; % m3->km3
      VFW(im,iyy)=vfw;

      vfws=nansum(VfwS(INGr))*1e-9;
      VFWstp(im,iyy)=vfws;

    end
  end
  fprintf('Saving %s\n',fextr);
  save(fextr,'VFW','VFWstp');
else
  fprintf('Loading %s\n',fextr);
  load(fextr);
end


[m,n]=size(VFW);
VFWy=mean(VFW);
VF=reshape(VFW,[m*n,1]);

% Step-function increase:
VFWSy=mean(VFWstp);
VFS=reshape(VFWstp,[m*n,1]);

% Platov
if f_plat==1
  [YRp,VFWp]=sub_read_platov;
  VFWyp=mean(VFWp);
end

% Myers:
if f_myr==1
		VFWym=[0.0824    0.2685    0.4654    0.6391    0.8471    1.0602    1.3174    1.5453, ...
									1.7227    1.8617    2.0157    2.1823    2.3531    2.5044    2.6081    2.7133, ...
									2.8187    2.9311    3.0190    3.1842    3.2929    3.3805    3.5057    3.6082]'*1e3;
		YRm=[1993:2016]';


end
stl=[];
ii=0;
YR=[YR1:2016];

figure(1); clf;
axes('Position',[0.1 0.5 0.85 0.42]);
hold on;
if f_plat==1
  plot(YRp,VFWyp,'k-','Linewidth',2.5,'Color',[0.9 0.4 0]);
  stl=sprintf('Tracer Vol GrFW anomaly SPNA, HYCOM max=%6.1f km^3, SibCIOM=%6.1f',...
      max(VFWy),max(VFWyp));
  ii=ii+1;
  strlg{ii}='SibCIOM';
end

if f_myr==1
  plot(YRm,VFWym,'k-','Linewidth',2.5,'Color',[0. 0.9 0.4]);
  stl=sprintf('Tracer Vol GFWA SPNA, HYCOM max=%6.1f km^3, SibCIOM=%6.1f, NEMO=%6.1f',...
      max(VFWy),max(VFWyp),max(VFWym));
  ii=ii+1;
  strlg{ii}='NEMO';

end

VFWSy=VFWSy-VFWSy(1); % make 0 in 1993
plot(YR,VFWy,'k-','Linewidth',2.5,'Color',[0.3 0.3 0.3]);  % lin incr
%plot(YR,VFWy,'k-','Linewidth',2.5,'Color',[0 0.6 0.8]);  % lin incr
%plot(YR,VFWSy,'k-','Linewidth',2.5,'Color',[0.9 0.5 0]); % Step fn

ii=ii+1;
strlg{ii}='HYCOM real';
ii=ii+1;
strlg{ii}='HYCOM step';
%lgd=legend(strlg);
%set(lgd,'Position',[0.2 0.2 0.14 0.13]);



stl1=sprintf('expt%i, Estimated FW Volume anomaly caused by GrFWFlux, km3',expt);
title(stl1);
%if isempty(stl)
%  stl=sprintf('Tracer-based Vol GFWA in SPNA, Bamber runoff and step-function , max=%6.1f, %6.1f km^3',...
%      max(VFWy), max(VFWSy));
%end


title(stl);
set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'ylim',[-100 3000],...
	'ytick',[0:500:4000],...
	'xtick',[1993:2016],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
btx='dS_FWC_timeseries_SubpolarGyre.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.4 0.05]);

    
axes('Position',[0.6 0.05 0.3 0.3]);
hold on
plot(XX(INGr),YY(INGr),'g.');
contour(HH,[0 0],'k');
axis('equal');
set(gca,'xlim',[300 1300],...
	'ylim',[50 1300],...
	'xtick',[],...
	'ytick',[]);


fprintf('\n ========================================= \n');
fprintf('To plot map of the region - dbcont\n');
keyboard


% Particles release location:
%fmatLP = sprintf('%sGrSh_part-lr10_1993-01.mat',pthmatLP);
%Ip = PRTCL.TRACK(1).I;
%Jp = PRTCL.TRACK(1).J;
slr=[];
fina=[];
finb=[];
Np0=1e4;
nsim=1;
[bmm,dmm,amm]=sub_seed_GrSh_prt(slr,HH,LON,LAT,fina,finb,Np0,nsim);
Ip = bmm.TRACK.I;
Jp = bmm.TRACK.J;

%
% Labrador Sea:
[bmm,dmm,amm]=sub_seed_WLabr_prt(slr,HH,LON,LAT,fina,finb,Np0,nsim);
Ilb = bmm.TRACK.I;
Jlb = bmm.TRACK.J;

% 
% Plot Greenland contour and study domain
%GC = sub_greenl_isobath(HH,LON,LAT); 
xyl=[300 10; ...
     1300  1200];
fn=2;
figure(2); clf;
set(gcf,'Position',[1078 377 937 955]);
axes('Position',[0.08 0.2 0.85 0.7]);

%sub_plot_bath3(HH,LON,LAT,fn,xyl);

xl1 = xyl(1,1);
yl1 = xyl(1,2);
xl2 = xyl(2,1);
yl2 = xyl(2,2);

LMSK = HH*0;
LMSK(HH<0)=1;
lcmp = [0.1 0.1 0.1; 1 1 1];
pcolor(LMSK); shading flat;
colormap(lcmp);
freezeColors;

hold on;

% Study Domain:
rcmp=[1 1 1; 0.8 1 0.8];
RG=HH*nan;
RG(INGr)=1;
%pcolor(RG); shading flat;
%colormap(rcmp);
%freezeColors;
IGR(end+1,:)=IGR(1,:);
%plot(IGR(:,1),IGR(:,2),'-','Linewidth',2,'Color',[0.9 0.5 0]);
plot(IGR(:,1),IGR(:,2),'-','Linewidth',2,'Color',[0.8 0.8 0.8]);

% Plot all gates:
%SCT = sub_define_SPNA_sections(HH,LON,LAT);
%nsct = length(SCT);
%for ip=1:nsct
%	clr=[0.9 0.5 0];
%		IIs=SCT(ip).I;
%		JJs=SCT(ip).J;
%		plot(IIs,JJs,'-',...
%							'Linewidth',2.5,'Color',clr);
%end



contour(HH,[-8000:1000:-100],'Linewidth',1,'Color',[0.8 0.8 0.8]);
contour(HH,[-1000 -1000],'Linewidth',1.6,'Color',[0.5 0.5 0.5]);
contour(HH,[-500 -500],'Linewidth',1,'Color',[0.8 0.8 0.8]);


%contour(HH,[-500:100:-1],'Color',[0.6 0.6 0.6]);
%contour(HH,[-50:10:-1],'Color',[0.4 0.4 0.4]);
axis('equal');
set(gca,'xlim',[xl1 xl2],'ylim',[yl1 yl2]);
%set(gcf,'Position',[666 98 1507 1241]);

% Plot Lagr Prt WLabr Shelf Particles
% release locations
%plot(Ip,Jp,'.','Color',[0.8 0.9 1]);
%plot(Ilb,Jlb,'.','Color',[0.8 0.9 1]);

% Plot Lon/Lat;
dmm=LON;
dmm(dmm<-170)=nan;
dmm(dmm>170)=nan;
dmm(LAT>88) = nan;
contour(dmm,[-180:30:180],'Color',[0.6 0.6 0.6]);
dmm=LON;
I=find(dmm<0);
dmm(I)=dmm(I)+360;
dmm(dmm>200)=nan;
dmm(dmm<160)=nan;
dmm(LAT>88)=nan;
contour(dmm,[150:30:190],'Color',[0.6 0.6 0.6]);
contour(LAT,[40:10:88],'Color',[0.6 0.6 0.6]);

bottom_text(btx,'pwd',1);







