% Average and save atl water characteristics
% extracted in extr_atlw_hycom004.m
% From mat to binary for python processing
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps/
startup

format long g
clear all
close


YR1 = 2020;
YR2 = 2020;

ixx    = 9; % experiment name and dir - check with EXPT - expt 023
%ixx    = 6;  % expt 022 original 
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;


rg = 9806;
hgg=1e20;

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/atl_water/',expt);

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

% Get topo:
ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[X,Y]=meshgrid([1:nn],[1:mm]);


% Overall mean
icc = 0;
for yr=YR1:YR2

  fmat=sprintf('%shycom004_%3.3i_%s_AtlLr_%4.4i.mat',...
                    pthmat,expt,texpt,yr);

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

end

Tmax=Tmax/icc;
Smax=Smax/icc;
Zmax=Zmax/icc;
dHL=dHL/icc;
HAtl=HAtl/icc;

Rmax=sw_dens0(Smax,Tmax);

nocn = length(Iocn);
fout = sprintf('%shycom004_%3.3i_AtlLr_%4.4i-%4.4i.dat',...
                    pthmat,expt,YR1,YR2); 
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







