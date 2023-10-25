% Plot T/S profiles
% Prepared in TSprofile_cotnour_greenl008
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

%pfld = 'salin'; % field to extract
%pfld = 'temp';
YR = 2005; % year to plot
im = 1;      % month to plot 



%expt = 110; % experiment without runoff
expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
regn = 'ARCc0.08';


rg=9806;  % convert pressure to depth, m
hgg=1e20; 
btx = 'TSprofile_contour_greenl008.m';


fprintf('Plotting T&S, expt %3.3i, mean %i\n',expt,YR);

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';
btx='plot_TSprofile_contour008.m';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section

fmat = sprintf('%s%3.3i_Greenl_TScontour_%i.mat',...
		   pthmat,expt,YR);

fprintf('Loading %s\n',fmat);
load(fmat);

dx = TSGR(1).DistCntr*1e-3; % m->km
Hs = TSGR(1).Hbottom;
%plot(dx,Hs); % plot Bottom profile along contour

%ar = ZZ*0;
T  = TSGR(im).T;
S  = TSGR(im).S;
TM = TSGR(im).TM;
dv = datevec(TM);
ZZ = TSGR(im).ZZ;

[nlr,npts]=size(S);

[Dst,dmb] = meshgrid(dx,[1:nlr]);

tstr=sprintf('arc008_%3.3i, T %4.4i/%2.2i',expt,dv(1:2));
% For plotting add surface layer      
ZZp = [ZZ(1,:);ZZ];
ZZp(1,:)=0;
Tp  = [T(1,:); T];
Dstp = [Dst(1,:);Dst];
c1=-2;
c2=2;
fnmb = 1;
pfld='temp';
Tp=sub_fill_bottom_nans(Tp);
ZZp(isnan(ZZp))=1.001*min(min(ZZp));
sub_plot_TS_Zcntr(Tp,ZZp,Dstp,fnmb,tstr,Hs,c1,c2,pfld);
bottom_text(btx,'pwd',1);

Sp  = [S(1,:); S];
c1=31;
c2=35;
fnmb = 2;
tstr=sprintf('arc008_%3.3i, S %4.4i/%2.2i',expt,dv(1:2));
pfld='salin';
Sp=sub_fill_bottom_nans(Sp);
sub_plot_TS_Zcntr(Sp,ZZp,Dstp,fnmb,tstr,Hs,c1,c2,pfld);
bottom_text(btx,'pwd',1);


