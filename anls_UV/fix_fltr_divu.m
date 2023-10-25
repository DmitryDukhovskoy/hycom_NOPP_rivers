% All fixed

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

YR=[1997:5:2016];
YR=[YR,2016];
nyr = length(YR);
zz0 = 50;
regn = 'ARCc0.08';
expt = 110;

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat2/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell = DX.*DY;

Ideep = find(HH<-500);

% Spatial filtering:
pgrd = 51;
Hmsk = HH;
Hmsk(HH<0)=1;
Hmsk(HH>=0)=0;
Hmsk(1920:end,:)=0;
Hmsk(100:500,1250:end)=0;
Hmsk(450:650,1280:1380)=0;
Hmsk(1:pgrd+1,:) = 0;
Hmsk(mm-pgrd:mm,:) = 0;
Hmsk(:,1:pgrd+1) = 0;
Hmsk(:,nn-pgrd:nn) = 0;
Imsk = find(Hmsk==1);


for ny = 1:nyr
  iyr = YR(ny);
  fmat = sprintf('%s%3.3i_divU%3.3im_v2_%i.mat',...
		   pthmat,expt,abs(zz0),iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  aa = DIVU(12).divU_m3_sec;
  divU = sub_fltr(aa,pgrd,Hmsk);
  DIVU(12).divU_m3_sec = divU;
  
  fprintf('Saving %s\n',fmat);
%  save(fmat,'DIVU');
 keyboard 
  
end;