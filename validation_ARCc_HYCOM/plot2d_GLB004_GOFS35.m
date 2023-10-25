% Vertical section in the Arctic Ocean
% GOFS3.5 reanalysis
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

dnmb = datenum(2017,07,01);
DV = datevec(dnmb);


s_fig = 0;
%pfld  = 'temp';
pfld  = 'salin';
fld0  = pfld;
lrplt = 1;  % layer to plot

hg    = 2^100;
rg    = 9806;
% Change Tv and nlev for different ARCc fields
Tv    = 11; % bathym v11
nlev  = 41; % 41

pthmat2 = '/Net/ocean/ddmitry/HYCOM/ARCc0.08/112/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);
btx='plot2d_GLB004_GOFS35.m';

%  cs1=30;
%  cs2=35;  
cs1=3.36;  % log scale S=28.79
cs2=3.56;  % log scale S=35.16
ct1=-2;
ct2=5;


% Subsample:
i1=5000;
i2=8000;
j1=5500;
j2=6900;

yr   = DV(1);
imo  = DV(2);
iday = DV(3);
jday = dnmb-datenum(yr,1,1)+1;

expt = 216;
pthycom = '/Net/kronos/ddmitry/hycom/GOFS3.5/';
fina  = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthycom,expt,yr,jday);
finb  = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthycom,expt,yr,jday);


% Get ARCc grid indices
flgrd = sprintf('%sregional.grid',pthycom);
fltopo = sprintf('%sdepth_GLBc0.04_27.a',pthycom);
GRD = read_grid_bath(flgrd,fltopo);
pln=GRD.PLON;
plt=GRD.PLAT;
hh=GRD.Topo;
I=find(hh>1e20);
hh=-1*hh;
hh(I)=100;

ln11=LN04(1,1);
ln21=LN04(end,1);
ln12=LN04(1,end);
ln22=LN04(end,end);

lt11=LT04(1,1);
lt21=LT04(end,1);
lt12=LT04(1,end);
lt22=LT04(end,end);

% ---------------------------------
% Read GLBc 0.04 GOFS3.5
% Read by layers
% subsample ARCc region
% and rotate the grid
% ---------------------------------
fprintf('Reading %s\n',datestr(dnmb));

ll=41;
SS=[];
%pfld='salin';
if strncmp(pfld,'salin',4);
		fprintf('Reading %s layer %i\n',pfld,lrplt);
		[Fr,n,m,l1] = read_hycom(fina,finb,pfld,'r_layer',lrplt);
		Fr=squeeze(Fr);
		Fr(Fr>hgg)=nan;
  S=Fr(j1:j2,i1:i2);
%		dmm = sub_Glb2Arc(Fr,IJ);
%		SS= dmm;

%  S = squeeze(SS(:,INDs));
  [a1,a2]=size(S);
  % Prepare for plotting - add extra bogus layer
  % at the bottom, Matlab won't show it
  S(end+1,1:a2)=S(end,:);
  AA=S;
  A0=S;
end

figure(1); clf;
pcolor(S); shading flat;
axis('equal');
set(gca,'xlim',[1 1500],...
        'ylim',[1 1400]);
caxis([25 34]);


TT=[];
%pfld='temp';
if strncmp(pfld,'temp',4)
		fprintf('Reading %s layer %i\n',pfld,lrplt);
		[Fr,n,m,l1] = read_hycom(fina,finb,pfld,'r_layer',ltplt);
		Fr=squeeze(Fr);
		Fr(Fr>hgg)=nan;
		dmm = sub_Glb2Arc(Fr,IJ);
		TT= dmm;
  T = squeeze(TT(:,INDs));
  [a1,a2]=size(T);
  % Prepare for plotting - add extra bogus layer
  % at the bottom, Matlab won't show it
  T(end+1,1:a2)=T(end,:);
  AA=T;
end



   
