% GOFS3.1 reanalysis
% Expreiments: 53.X
% gridded netcdf 3hrly 
% Calculate mean T/S in Greenland region
% For Greenland - use natie Glb grid
% for Arctic - use extracted and rotated to ARCc grid
% files
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 1;
pfld = 'salin';
%pfld = 'temp';
%mavg = 'DJF'; % winter/summer mean T/S
mavg = 'JJA'; % winter/summer mean T/S
zplt = 0; % depth to plot

rg=9806;  % convert pressure to depth, m
TV = '07';

yr1=2005;
yr2=2005;

pth53X='/nexsan/hycom/GLBv0.08_53X/data/';
%pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % 1993-1994
%pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_fw/';
pthtopo = '/nexsan/hycom/GLBv0.08_53X/topo/';
fmat = sprintf('%sGOFS31_%s_natl_%4.4im_%i-%i_%s.mat',pthmat,pfld,abs(zplt),yr1,yr2,mavg);

% Greenland Region:
lon1 = -100;
lon2 = 40;
lat1 = 45;
lat2 = 88;


% Calculate time-mean:
Sav =[];
ncc =0;
for iyr = yr1:yr2
%  yr = iyr;
%  nday = datenum(yr,12,31)-datenum(yr,1,1)+1;
  
  switch(mavg)
   case('DJF')
    nday=datenum(iyr+1,2,28)-datenum(iyr,12,1)+1;
    id1=datenum(iyr,12,1)-datenum(iyr,1,1)+1;
   case('JJA')
    nday=datenum(iyr,8,31)-datenum(iyr,6,1)+1;
    id1=datenum(iyr,6,1)-datenum(iyr,1,1)+1;
  end  
  dnmb0 = datenum(iyr,1,1)+id1-1;
  
  for idd = 1:5:nday
    tic; 
    dnmb = dnmb0+(idd-1);
    dv   = datevec(dnmb);
    yr   = dv(1);
    mo   = dv(2);
    mday = dv(3);
    iday = dnmb-datenum(yr,1,1)+1;
    
    expt = sub_exptGLBb(dnmb);
    fprintf('Getting data for %s\n',datestr(dnmb));
    
    switch(pfld)
     case('salin');
      flde = 's';
      fldnc = 'salinity';
     case('temp');
      flde = 't';
      fldnc= 'water_temp';
    end

    expt=533;
    if dnmb>datenum(2005,6,30)
      expt=534;
    end
    
    hr=0;
    fin=sprintf('%s%4.4i/hycom_GLBv0.08_%3.3i_%4.4i%2.2i%2.2i12_t%3.3i.nc',...
		pth53X,yr,expt,yr,mo,mday,hr);
    lla=exist(fin,'file');
    while ~lla
      hr=hr+3;
      fin=sprintf('%s%4.4i/hycom_GLBv0.08_%3.3i_%4.4i%2.2i%2.2i12_t%3.3i.nc',...
		pth53X,yr,expt,yr,mo,mday,hr);
      lla=exist(fin,'file');
    end
    
    if ~lla
      fprintf('File %s not found\n');
      continue
    end
    
    
    
    fprintf('Opening %s\n',fin);
    
    
    if ~exist('ZZ','var')
      ZZ = -(nc_varget(fin,'depth'));
      DZ = abs(diff(ZZ));
      ll = length(ZZ);
      iz0= max(find(ZZ>=zplt));
    end


    if ~exist('idx1','var');
      LON = nc_varget(fin,'lon');
      LAT = nc_varget(fin,'lat');
      imx=find(LON>180);
      while ~isempty(imx),
	LON(imx)=LON(imx)-360;
        imx=find(LON>180);
      end	
      
% Find Greenl region: 
      d=abs(LAT-lat1);
      jdx1=find(d==min(d),1);
      d=abs(LON-lon1);
      idx1=find(d==min(d),1);
      
      d=abs(LAT-lat2);
      jdx2=find(d==min(d),1);
      d=abs(LON-lon1);
      idx2=find(d==min(d),1);
      
      d=abs(LAT-lat2);
      jdx3=find(d==min(d),1);
      d=abs(LON-lon2);
      idx3=find(d==min(d),1);
      
      d=abs(LAT-lat1);
      jdx4=find(d==min(d),1);
      d=abs(LON-lon2);
      idx4=find(d==min(d),1);

      nip=idx3-idx1+1;
      njp=jdx2-jdx1+1;

    end

    S = squeeze(nc_varget(fin,fldnc,[0 iz0-1 jdx1-1 idx1-1],[1 1 njp nip]));
    if isempty(Sav);
      Sav=S*0;
    end
    ncc=ncc+1;
    Sav=Sav+S;
    
    smm = Sav./ncc;
    fprintf('Averaged: Smax=%6.4f, Smin=%6.4f\n',max(max(smm)),min(min(smm)));
    fprintf('1 rec: %6.4f min\n\n',toc/60);
  end
  
end
Sav=Sav./ncc;
SAV.regn    = 'Greenland';
SAV.YR      = [yr1, yr2];
SAV.Avrg_Mo = mavg;
SAV.S_mean  = Sav;

fprintf('Saving mean %s: %s\n',pfld,fmat);
save(fmat,'SAV');




