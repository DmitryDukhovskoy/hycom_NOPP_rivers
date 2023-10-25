% Plot mean EKE
% calculated in meanEKE.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=1993;
YR2=2016;
zz1=0;
zz2=-50;

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
pthout='/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';


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
mnEKE=A./ncc;

fmat=sprintf('%shycom008_meanEKE',pthout);
fprintf('Saving %s\n',fmat);
save(fmat,'mnEKE');



