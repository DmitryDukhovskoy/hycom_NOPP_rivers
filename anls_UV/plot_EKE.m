% Plot mean EKE
% calculated in meanEKE.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;

YR1=2005;
YR2=2009;
zz1=0;
%zz2=-50;
zz2=-15;

av=0; % =1 - average over specofoed years, =0 - plot individual years
plr =1;  % U from plr layer
rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.08';
expt = 110;
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/%3.3i/fig_meanUV/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat2 =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);


%fall=1; % overall mean EKE from mean <u'^2>,<v'2>, not same as mean[1993:2016]!

fprintf('Plotting mean EKE, %i-%i, %i-%i m\n',abs(zz1),abs(zz2),YR1,YR2);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);

ncc=0;
for iyr=YR1:YR2
  fout = sprintf('%sarc08_%i_annualEKE_%4.4i-%4.4i_%i.mat',...
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
xlim1 = 200;
xlim2 = 1300;
ylim1 = 50;
ylim2 = 1150;

%S. Greenland
xlim1 = 220;
xlim2 = 980;
ylim1 = 175;
ylim2 = 935;


stl = sprintf('0.08-%i, Mean EKE, cm2/s2,  %i-%i m,  %i-%i',...
	      expt,abs(zz1), abs(zz2),YR1,YR2);
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
% set(gcf,'Position',[1315         193        1058        1137]);
pfld='eke';
hps = [0.92 0.2 0.02 0.7];
nf=1;
ekem(ylim2:end,:)=nan;
Lekem=log(ekem);
ekem_cm2=ekem*1e4;  % cm2/s2
sub_plot_scalar(ekem_cm2,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
		'figpos',[1635, 535, 790, 765],'c1',0,'c2',400,'clbpos',hps,'cmp',6);

% Log scale
%sub_plot_scalar(Lekem,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
%    'c1',-7,'c2',-3,'clbpos',hps,'cmp',5);
%contour(ekem,[10:10:100],'k');

btx='plot_EKE.m';

bottom_text(btx,'pwd',1,'position',[0.05 0.08 0.4 0.04]);

if s_fig==1
  fgnm=sprintf('%sarc08_%3.3i_meanEKE_%i-%i',pthfig,expt,YR1,YR2);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r200',fgnm);
end  



