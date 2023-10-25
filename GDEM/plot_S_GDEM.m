% Plot S field
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

iz0 = 2;  % layer to plot
mplt=[6,7,8]; % Average over these months

pthmat = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/023/atl_water/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthin = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/data_mat/';  % GDEM fields


ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
[X,Y]=meshgrid([1:nn],[1:mm]);

fprintf('Averaging S lr=%i over months: %i - %i\n',iz0,mplt(1),mplt(end));

icc = 0;
for imo=mplt(1):mplt(end)
  ftin = sprintf('%sTave_GDEM_%2.2i.mat',pthin,imo);
  fsin = sprintf('%sSave_GDEM_%2.2i.mat',pthin,imo);

  if ~exist(ftin,'file')
    fprintf('Input fields for %i missing, skipping ...\n',imo);
    continue;
  end

  fprintf('Processing month = %i\n',imo);

  tic;
%  fprintf('Loading %s\n',ftin);
%  TAV=load(ftin);
  fprintf('Loading %s\n',fsin);
  SAV=load(fsin);

%  T = TAV.Tav;
  S = SAV.Sav;
%
% Bug in interpolation code - smaller T,S domain
  [nlr,mt,nt] = size(S);
  if mt~=mm | nt~=nn
    S(:,mt+1:mm,:) = nan;
    S(:,:,nt+1:nn) = nan;
  end
  
  icc=icc+1;
  if (icc==1)
    AA = squeeze(S(iz0,:,:));
  else
    AA = AA+squeeze(S(iz0,:,:));
  end
end 

S = AA/icc;
S(HH>=0)=nan;
fprintf('Total records = %i\n',icc);
fprintf('Min/max Value = %6.2f / %6.2f\n',nanmin(nanmin(S)),nanmax(nanmax(S)));

%
% Do Gaussian smoothing for plotting
Hmsk = HH;
Hmsk(HH<0) = 1;
Hmsk(HH>=0) = 0;
Hmsk(1150:end,:)=0;
SS = S;
pgrd=33;
S = sub_fltr(SS,pgrd,Hmsk);


% subpolar NA
xlim1 = 370;
xlim2 = 1250;
ylim1 = 150;
ylim2 = 1100;
% Greenland
%xlim1 = 450;
%xlim2 = 1100;
%ylim1 = 300;
%ylim2 = 1100;
% SE Greenland:
%xlim1 = 600;
%xlim2 = 860;
%ylim1 = 350;
%ylim2 = 640;

fprintf('Plotting ...\n');
nf = 1;
sb=33;
im1=mplt(1);
im2=mplt(end);
stl = sprintf('GDEMv4 Mean S (cntr0=%3.1f) mo:%i-%i',sb,im1,im2);
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
pfld='salin';
hps = [0.93 0.1 0.02 0.8];
Fpos = [1350  255  1034  1085]; % Figure gcf position
sub_plot_scalar(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
    'c1',30,'c2',35,'cmp',2,'clbpos',hps,'figpos',Fpos);
contour(S,[30:0.5:36],'k');
contour(S,[sb sb],'k','Linewidth',1.8);

txtb = 'plot_S_GDEM.m';
bottom_text(txtb,'pwd',1);




