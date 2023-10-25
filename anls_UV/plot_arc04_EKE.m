% Plot mean EKE, HYCOM ARCc0.04
%
% calculated in meanEKE_arc04.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;

YR1=2005;
YR2=2007;
zz1=0;
zz2=-50;

av=0; % =1 - average over specofoed years, =0 - plot individual years
plr =1;  % U from plr layer
rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.04';
expt = 011;
pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.04/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/data_mat/',expt);
pthmat2 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/data_mat2/',expt);

%fall=1; % overall mean EKE from mean <u'^2>,<v'2>, not same as mean[1993:2016]!

fprintf('Plotting mean EKE %s, %i-%i, %i-%i m\n',regn,abs(zz1),abs(zz2),YR1,YR2);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);

ncc=0;
for iyr=YR1:YR2
  fout = sprintf('%sarc04_%i_annualEKE_%4.4i-%4.4i_%i.mat',...
	       pthmat2,expt,abs(zz1),abs(zz2),iyr);
  fprintf('Loading %s\n',fout);
  load(fout);

  eke=EKE.EKE_m2_s2;
  if ncc==0
    A=eke*0;
  end
  ncc=ncc+1;
  A=A+eke;
end
ekem=A./ncc;

% Greenland
xlim1 = 400;
xlim2 = 2600;
ylim1 = 100;
ylim2 = 2300;

stl = sprintf('0.04-%i, Annual Mean EKE, m2/s2,  %i-%i m,  %i-%i',...
	      expt,abs(zz1), abs(zz2),YR1,YR2);
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
pfld='eke';
hps = [0.92 0.2 0.035 0.7];
nf=1;
ekem(ylim2:end,:)=nan;
Lekem=log(ekem);
sub_plot_scalar(Lekem,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
		'c1',-7,'c2',-2,'clbpos',hps,'cmp',5);
%contour(ekem,[10:10:100],'k');

btx='plot_EKE.m';
bottom_text(btx,'pwd',1,'position',[0.05 0.08 0.4 0.04]);





