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


ys=1997;  % 1st year
ye=2016;  % last year (year of Oct-Dec segment)
%% To Plot wind roses for select years and regions:
% make f_roses2=1 and select years when call subfunction


ndate0=datenum(ys,1,1,0,0,0);
ndate1=datenum(ye,12,31,0,0,0);     % last date

dv1=datevec(ndate0);
dv2=datevec(ndate1);

fprintf('Time range: %i/%i/%i - %i/%i/%i\n',dv1(1:3),dv2(1:3));

pthdat = '/Net/data/ccmp/v02.0/';
pthdt  = '/Net/tholia/ddmitry/ovwst/data_mat/';
pth8   = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig = '/Net/tholia/ddmitry/ovwst/fig_cycl/';
pthtopo= '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

fmat = sprintf('%swrose_Greenland_ccmp020rss_%i_%i.mat',pthdt,ys,ye);
%fmat = sprintf('%swrose_Greenland_ccmp020rss_%i_%iv1.mat',pthdt,ys,ye);

btx = 'wind_ccmp_GrShelf_anls.m';

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

% Regions of interest: SW Greenland shelf
% winds are in math convention U + east, V + north
ir=7;
TM=WRS(1).TM;

id1=datenum(2000,1,1);
id2=datenum(2016,12,31);
I1=find(TM==id1);
I2=find(TM==id2);

U=WRS(ir).U(I1:I2);
V=WRS(ir).V(I1:I2);
TM=TM(I1:I2);
DV=datevec(TM);
yr1=DV(1,1);
yr2=DV(end,1);


% Monthly means:
cc=0;
for iyr=yr1:yr2
  for im=1:12
    I=find(DV(:,1)==iyr & DV(:,2)==im);
    um=mean(U(I));
    vm=mean(V(I));
    cc=cc+1;
    UM(cc)=um;
    VM(cc)=vm;
  end
end

WND.Region=ir;
WND.Indx=WRS(ir).Indices;
WND.Info=sprintf('Monthly U,V for %i-%i',yr1,yr2);
WND.Info2='Source: hycom_NOPP_rivers/anls_atms/wind_ccmp_GrShelf_anls.m';
WND.U=UM;
WND.V=VM;
foutp=sprintf('%swind_ccmp_rgn%2.2i_GrSh.mat',pthdt,ir);
fprintf('Saving %s\n',foutp);
save(foutp,'WND');

bottom_text(btx,'pwd',1);
