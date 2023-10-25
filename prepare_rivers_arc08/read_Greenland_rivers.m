% Read Greenland RUnoff data
%The grid is at 5 km postings and is on a polar stereographic
%projection with standard parallel at 71 degrees. Other details of the projection
%are provided below:
% The grid is 301 x 561 in size. 
%Points are grid corner not grid centre. The file
% grn_5km_lat_lons.txt contains the 
% lon, lat values for each grid cell ordered in rows.
% The values at each grid point are in km^3 water/month. The data run from
% Jan 1958 to Dec 2010, which equals 53*12 (636) monthly grids. The FWF is located
% at the nearest ocean grid point to where the water is routed from the ice sheet.
% This may be within a fjord or it may be at the edge of the island mass.
%
% Bamber's data cover period until 2010
% afr that, I'm adding a trend estimated
% from Bamber's total Greenland Mass Change
%  ~-290 GT/yr for 2010-2014


addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

clear all
close

f_save=1;  % save data in mat files

PTH.river='/Net/ocean/ddmitry/arctic_AOregimes/data/Greenland_rivers/';
fgrd=sprintf('%sgrn_5km_lat_lons.txt',PTH.river);
friv=sprintf('%sGreenland_runoff_calving_5km_monthly.bin',PTH.river);

nn=561;
mm=301;
YS=1958;
YE=2010;
nrec=(YE-YS+1)*12;
% 
dmm=load(fgrd);
LN=reshape(dmm(:,1),[mm,nn]);
LN=LN';
LT=reshape(dmm(:,2),[mm,nn]);
LT=LT';

if f_save==1
  GRgrd=struct;
  GRgrd.LN=LN;
  GRgrd.LT=LT;
  fout=sprintf('%sGreenland_grid.mat',PTH.river);
  fprintf('Saving grid: %s\n',fout);
  save(fout,'GRgrd');
end

fid=fopen(friv,'r','ieee-le');
frewind(fid);

mday=[31;28;31;30;31;30;31;31;30;31;30;31];

for irec=1:nrec
  A = fread(fid,[mm,nn],'float');
  A=A';
  A(A==0)=nan;
  YC=YS+floor((irec-1)/12);
  MC=mod(irec,12);
  nm=MC;
  if MC==0, nm=12; end;
  if MC==1,
    clear GR
    GR = struct;
  end
  
  fprintf('Year=%i, Mo=%i\n',YC,nm);
  
  GR(nm).year=YC;
  GR(nm).units='km3/mo';
  GR(nm).runoff=A;
  if MC==0 & f_save==1;
    fout=sprintf('%sGreenland_runoff_monthly-%4.4i.mat',PTH.river,YC);
    fprintf('Saving runoff: %s\n',fout);
    save(fout,'GR');
  end
 
%  keyboard
  
  
  fplt=logical(0);
% Convert km3/mo to m3/sec:
  km3mo=1; % 1km3/mo
  ndays=30.4167;
  mSv=km3mo*1e9/(3600*24*ndays)*1e-3;  % mSv

  if fplt
    A(A==0)=nan;
    clf
    pcolor(A); shading flat;
    axis('equal');
    hold on;
    contour(LN,[-100:5:100],'Color',[0.4 0.4 0.4]);
    contour(LT,[40:5:90],'Color',[0.4 0.4 0.4]);
    caxis([0 0.1]);
    set(gca,'color',[0 0 0]);
    set(gcf,'Color',[1 1 1]);
    stt=sprintf('Runoff, km3/mo, %4.4i %2.2i',YC,nm);
    title(stt,'Fontsize',14);
    colorbar;
  end

end;
fclose(fid);





