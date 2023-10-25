% Calculate mean T/S in Greenland region
% For Greenland - use natie Glb grid
% for Arctic - use extracted and rotated to ARCc grid
% files
% /nexsan/GLBa0.08/expt_90.9/data
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
zplt = -5; % depth to plot

rg=9806;  % convert pressure to depth, m
TV = '07';

yr1=2015;
yr2=2015;

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_fw/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
%fmat = sprintf('%sFWC_natl_1993-2015.mat',pthmat);
%fmat = sprintf('%sFWC_natl_1993-2016.mat',pthmat);
fmat = sprintf('%s%s_natl_%4.4im_%i-%i_%s.mat',pthmat,pfld,abs(zplt),yr1,yr2,mavg);

% Greenland Region:
lon1 = -80;
lon2 = 10;
lat1 = 50;
lat2 = 85;


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
      flddr = 'salt';
     case('temp');
      flde = 't';
      flddr= 'temp';
      fldnc= 'temperature';
    end
    
    pthbin = sprintf('/nexsan/GLBa0.08/expt_%2.1f/data/%i/%s/',...
		    expt/10,yr,flddr);
 
    if expt<910 & expt>906,
      pthbin=sprintf('/nexsan/GLBa0.08/expt_%2.1f/data/%s/',...
		     expt/10,flddr);
    end
    if expt==906,
      pthbin=sprintf('/nexsan/GLBa0.08/GLBa0.08_906/data/%s/',flddr);
    end
    
    fll = sprintf('archv.%i_%3.3i_00_3z%s.nc4',...
		  yr,iday,flde);
    flnm = sprintf('%s%s',pthbin,fll);
    if ~exist(flnm,'file'),
      fll = sprintf('archv.%i_%3.3i_00_3z%s.nc',...
		  yr,iday,flde);
      flnm = sprintf('%s%s',pthbin,fll);
    end
    
    if ~exist(flnm,'file'),
      fprintf('FIle does not exist, skipping %s\n',flnm);
      continue
    end
    
    fprintf('Opening %s\n',flnm);
    
    if ~exist('ZZ','var')
      ZZ = -(nc_varget(flnm,'Depth'));
      DZ = abs(diff(ZZ));
      ll = length(ZZ);
      iz0= max(find(ZZ>=zplt));
    end


    if ~exist('idx1','var');
      LON = nc_varget(flnm,'Longitude');
      LAT = nc_varget(flnm,'Latitude');
      imx=find(LON>180);
      while ~isempty(imx),
	LON(imx)=LON(imx)-360;
        imx=find(LON>180);
      end	
      
% Find Greenl region: 
      d = sqrt((LAT-lat1).^2+(LON-lon1).^2);
      [jdx1,idx1] = find(d==min(min(d)));
      d = sqrt((LAT-lat2).^2+(LON-lon1).^2);
      [jdx2,idx2] = find(d==min(min(d)));
      d = sqrt((LAT-lat2).^2+(LON-lon2).^2);
      [jdx3,idx3] = find(d==min(min(d)));
      d = sqrt((LAT-lat1).^2+(LON-lon2).^2);
      [jdx4,idx4] = find(d==min(min(d)));
% midpoints for better region boundaries:
      d = sqrt((LAT-lat1).^2+(LON-0.5*(lon1+lon2)).^2);
      [jdx5,idx5] = find(d==min(min(d)));
      d = sqrt((LAT-0.5*(lat1+lat2)).^2+(LON-lon1).^2);
      [jdx6,idx6] = find(d==min(min(d)));
      d = sqrt((LAT-lat2).^2+(LON-0.5*(lon1+lon2)).^2);
      [jdx7,idx7] = find(d==min(min(d)));
      d = sqrt((LAT-0.5*(lat1+lat2)).^2+(LON-lon2).^2);
      [jdx8,idx8] = find(d==min(min(d)));

      vx=[idx1,idx6,idx2,idx7,idx3,idx8,idx4,idx5,idx1];
      vy=[jdx1,jdx6,jdx2,jdx7,jdx3,jdx8,jdx4,jdx5,jdx1];
      
      [mm,nn] = size(LON);
      [XX,YY] = meshgrid((1:nn),(1:mm));
      INP  = inpolygon(XX,YY,vx,vy);
      iBG  = find(INP==1);
      INAN = find(INP==0);
      [J,I] = ind2sub(size(LON),iBG);
      J1 = min(J);
      J2 = max(J);
      I1 = min(I);
      I2 = max(I);
      nj = J2-J1+1;
      ni = I2-I1+1;
      
      f_chck=0;
      if f_chck==1
	figure(10); clf;
	contour(LAT,[40:5:89],'b');
	hold on;
	contour(LON,[-180:5:179],'b');
	plot(vx,vy,'m*');
	contour(LAT,[lat1 lat1],'r');
	contour(LAT,[lat2 lat2],'r');
	contour(LON,[lon1 lon1],'r');
	contour(LON,[lon2 lon2],'r');
        plot(I,J,'g.');
      end
      
    end

    S = squeeze(nc_varget(flnm,fldnc,[0 iz0-1 J1-1 I1-1],[1 1 nj ni]));
    if isempty(Sav);
      Sav=S*0;
    end
    ncc=ncc+1;
    Sav=Sav+S;
    
    smm = Sav./ncc;
    fprintf('Averaged: Smax=%6.4f, Smin=%6.4f\n',max(max(smm)),min(min(smm)));
    fprintf('1 rec: %6.4 min\n\n',toc/60);
  end
  
end
Sav=Sav./ncc;
SAV.regn    = 'Greenland';
SAV.YR      = [yr1, yr2];
SAV.Avrg_Mo = mavg;
SAV.S_mean  = Sav;

fprintf('Saving mean %s: %s\n',pfld,fmat);
save(fmat,'SAV');




