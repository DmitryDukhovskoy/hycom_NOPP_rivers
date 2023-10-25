% Plot ccmp winds at select locations
% on Gr. Shelf
%
% Wind stat prepared in wind_stat_ccmp_Greenland.m
% CCMP winds for Greenland coast
% speed - mean scalar,
% direction - mean vector, i.e. predominant direction
% for different regions
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08/;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/;
startup

clear
close all

% Dates to plot:
ir=5;  % 5 - SE Gr Shelf Box
thtD=78; % coastline direction wrt to 
         % East where tht=0, 
	 % negative projection = downwelling

dstr=datenum(2005,1,1);  % 1st year
dend=datenum(2008,12,31);  % last year (year of Oct-Dec segment)

dv1=datevec(dstr);
dv2=datevec(dend);

fprintf('Time range: %i/%i/%i - %i/%i/%i\n',dv1(1:3),dv2(1:3));

pthdat = '/Net/data/ccmp/v02.0/';
pthdt  = '/Net/tholia/ddmitry/ovwst/data_mat/';
pth8   = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig = '/Net/tholia/ddmitry/ovwst/fig_cycl/';
pthtopo= '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

fmat = sprintf('%swrose_Greenland_ccmp020rss_1997_2016.mat',pthdt);
%fmat = sprintf('%swrose_Greenland_ccmp020rss_%i_%iv1.mat',pthdt,ys,ye);

btx = 'plot_ccmp_GrShelf.m';

dnmb=datenum(2005,1,1);
d0vct=datevec(dnmb);
year=d0vct(1);
mnth=d0vct(2);
mday=d0vct(3);
fp=sprintf('%sY%i/M%2.2i/CCMP_Wind_Analysis_%4.4i%2.2i%2.2i_V02.0_L3.0_RSS.nc',...
	       pthdat,year,mnth,year,mnth,mday);

%x1=210;
x1=270;  % Westernmost coord
x2=360;  % Easternmost coord.
y1=50;
y2=78.375;


% CCMP Ocean Wind vector data:
alat=nc_varget(fp,'latitude');
elon=nc_varget(fp,'longitude');
n=length(elon);
m=length(alat);


dx=abs(elon(2)-elon(1));
if x1<0
  x1a=x1+360;
else
  x1a=x1;
end

dst=sqrt((elon-x1a).^2);
i1=min(find(dst==min(dst)));
dst=sqrt((alat-y1).^2);
j1=min(find(dst==min(dst)));
dst=sqrt((elon-x2).^2);
i2=min(find(dst==min(dst)));
dst=sqrt((alat-y2).^2);
j2=min(find(dst==min(dst)));

[X,Y] = meshgrid(elon(i1:i2),alat(j1:j2));
subalat=alat(j1:j2);
if i2<i1
  e1=elon(i1:end);
  e2=elon(1:i2);
  subelon=[e1;e2];
else
 subelon=elon(i1:i2);
end


% HYCOM grid/topo for plotting
ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);



fprintf('Loading %s\n',fmat);
load(fmat); % WCLIM WRS

%keyboard

% Regions of interest: SW Greenland shelf
% winds are in math convention U + east, V + north
TM=WRS(1).TM;

I1=find(TM==dstr);
I2=find(TM==dend);

%keyboard

U=WRS(ir).U(I1:I2);
V=WRS(ir).V(I1:I2);
TM=TM(I1:I2);
DV=datevec(TM);
yr1=DV(1,1);
yr2=DV(end,1);

% Upwelling favorable winds are negative
% Project winds on coast direction
WDir=atan2d(V,U);
alfa=thtD-WDir;
S=sqrt(U.^2+V.^2);
%Pcst=[cosd(thtD), sind(thtD)];
%  Wcst=[U,V]*Pcst'; or:
Wcst=S.*cosd(alfa);


% Filter in time
% OUtput freq q=1/4 day => 
% Cutoff freq Wn=1 corresp to half of q
% i.e. 1/8
% 7day cutoff = (1/(7*4))/(1/8)
q=1/4;
toff_day=7;
qoff=1/toff_day*q/(q/2);
Wn = qoff;
[Bf,Af] = butter(9,Wn,'low');
WcstF = filtfilt(Bf,Af,Wcst); % W/m

% ===================
% Monthly means:
% ===================
YRS=[];
for iyr=yr1:yr2
  J=find(DV(:,1)==iyr);
  nj=length(J);
  y0=DV(J(1));
  dmm=[y0:1/nj:y0+(nj-1)/nj];
  YRS=[YRS,dmm];
end
YRS=YRS(:);

YM=[];  % time array for plotting means
cc=0;
for iyr=yr1:yr2
  for im=1:12
    I=find(DV(:,1)==iyr & DV(:,2)==im);
    wm=mean(Wcst(I));
    cc=cc+1;
    WcstM(cc)=wm;
    dmm=YRS(I);
    YM(cc,1)=dmm(1);
    YM(cc,2)=dmm(end);
  end
end


ctt=0;
ytck=[];
for yr=DV(1,1):DV(end,1)
  for im=1:12
    ii=find(DV(:,1)==yr & DV(:,2)==im,1);
    ctt=ctt+1;
    ytck(ctt)=YRS(ii);
    ylbl{ctt}=sprintf('%2.2i/%4.4i',im,yr);
  end
end

figure(2); clf;
axes('Position',[0.2 0.1 0.25 0.8]);
plot(WcstF,YRS);
hold on;
nm=length(WcstM);
for jt=1:nm
  w0=WcstM(jt);
  w1=w0;
  if jt<nm, w1=WcstM(jt+1); end;
  t1=YM(jt,1);
  t2=YM(jt,2);
  plot([w0 w0],[t1 t2],'Color',[0.8 0.3 0],'Linewidth',2);
  plot([w0 w1],[t2 t2],'Color',[0.8 0.3 0],'Linewidth',2);
end
plot([0 0],[YRS(1) YRS(end)],'k--');

set(gca,'tickdir','out',...
	'xlim',[-20 20],...
	'xtick',[-20:5:20],...
	'ylim',[YRS(1) YRS(end)],...
	'ytick',ytck,...
	'yticklabel',ylbl,...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
xlabel('Wind Proj., m/s');
title('CCMPv2.0, "-" downwl., SE Greenl. Shelf');

axes('Position',[0.6 0.8 0.3 0.1]);
stl=sprintf('Filtered cutoff T=%3.1f day & Monthly',toff_day);
text(0,0,stl,'Fontsize',12);
set(gca,'Visible','off');


bottom_text(btx,'pwd',1);
