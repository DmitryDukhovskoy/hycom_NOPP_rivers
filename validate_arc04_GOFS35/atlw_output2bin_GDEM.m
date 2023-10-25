% Average and save atl water characteristics
% extracted in /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/GDEM/
% extr_atlw_GDEM.m
% From mat to binary for python processing
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps/
startup

format long g
clear all
close


pthmat = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/023/atl_water/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

% Get topo:
ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[X,Y]=meshgrid([1:nn],[1:mm]);


% Overall mean
% for all months
icc = 0;
fmat=sprintf('%sGDEM_AtlLr.mat',pthmat);
fprintf('Loading %s\n',fmat);
load(fmat);
imE=length(ATL);
fprintf(' Found %i months\n',imE);

if icc==0
	Iocn = ATL(1).Iocn;
	Tmax = Iocn*0; % Tmax
	Smax = Iocn*0; % S at Tmax
	Zmax = Iocn*0; % depth of T max
	dHL  = Iocn*0; % thickness Atl. layer
	HAtl = Iocn*0;  % heat content in Atl L.
end

for im=1:imE
	icc=icc+1;
	Tmax = Tmax+ATL(im).Tmax;
	Smax = Smax+ATL(im).Smax;
	Zmax = Zmax+ATL(im).Zmax;
	dHL  = dHL+ATL(im).LrThck;
	HAtl = HAtl+ATL(im).HtCnt;
end


Tmax=Tmax/icc;
Smax=Smax/icc;
Zmax=Zmax/icc;
dHL=dHL/icc;
HAtl=HAtl/icc;

Rmax=sw_dens0(Smax,Tmax);

nocn = length(Iocn);
fout = sprintf('%sGDEM_AtlLr_overall.dat',pthmat); 
fid = fopen(fout,'wb');
fwrite(fid,nocn,'int32');
fwrite(fid,Iocn,'int32');
fwrite(fid,Zmax,'float');
fwrite(fid,Tmax,'float');
fwrite(fid,Rmax,'float');
fwrite(fid,dHL,'float');
fwrite(fid,HAtl,'float');
fclose(fid);

fprintf('Output is saved in %s\n',fout);


f_chck=0;
if f_chck==1
  Da = HH*0;
  Da(Iocn) = Zmax;
  Da(HH>=0)=nan;

  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;

  pcolor(Da); shading flat;
  
end







