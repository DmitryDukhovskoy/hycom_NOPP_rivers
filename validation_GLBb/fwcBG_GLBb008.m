% Analyze FWC of the BG 
% from the Global HYCOM analysis 0.08
% expt 9[01].X 
% Use interpolated onto Z, Global grid
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 1; % 
s_fig  = 0;
fld    = 'salt';

fprintf('saving mat = %i\n',s_mat);

rg=9806;  % convert pressure to depth, m
Sref=34.8; % N.Atl. is too saline
TV = '07';

pthARC = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mean/';
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
%ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';

% BG coordinates:
lat1 = 73;
lat2 = 83;
lon1 = -160;
lon2 = -130;

for iyr = 2017:2017
  yr = iyr;
  nday = datenum(yr,12,31)-datenum(yr,1,1)+1;
  cc = 0;
  
  for iday = 1:5:nday
    tic; 
    dnmb = datenum(yr,1,1)+iday-1;
    expt = sub_exptGLBb(dnmb);
    fprintf('Getting data for %s\n',datestr(dnmb));
    
    switch(fld)
     case('salt');
      flde = 's';
      fldnc = 'salinity';
     case('temp');
      flde = 't';
    end
    
    pthbin = sprintf('/nexsan/GLBa0.08/expt_%2.1f/data/%i/%s/',...
		    expt/10,yr,fld);
 
    if expt<910 & expt>906,
      pthbin=sprintf('/nexsan/GLBa0.08/expt_%2.1f/data/%s/',...
		     expt/10,fld);
    end
    if expt==906,
      pthbin=sprintf('/nexsan/GLBa0.08/GLBa0.08_906/data/%s/',fld);
    end
    
    fll = sprintf('archv.%i_%3.3i_00_3z%s.nc',...
		  yr,iday,flde);
    flnm = sprintf('%s%s',pthbin,fll);

    if ~exist(flnm,'file'),
      fprintf('FIle does not exist, skipping %s\n',flnm);
      continue
    end
    
    fprintf('Opening %s\n',flnm);
    
    if ~exist('ZZ','var')
      ZZ = -(nc_varget(flnm,'Depth'));
      DZ = abs(diff(ZZ));
      ll = length(ZZ);
    end
    
    if ~exist('idx1','var');
      LON = nc_varget(flnm,'Longitude');
      LAT = nc_varget(flnm,'Latitude');
      imx=find(LON>180);
      while ~isempty(imx),
	LON(imx)=LON(imx)-360;
        imx=find(LON>180);
      end	
      
      I1 = 1000;
      I2 = 3500;
      J1 = 2500;
      J2 = size(LAT,1);
      
      ni = I2-I1+1;
      nj = J2-J1+1;
      
      LON = LON(J1:J2,I1:I2);
      LAT = LAT(J1:J2,I1:I2);
      
% Find BG region: 
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
      
%      fprintf('Cell area\n');
      [DX,DY] = sub_dx_dy(LON,LAT);
      Acell = DX.*DY*1e-6; % km2
    end

    S = squeeze(nc_varget(flnm,fldnc,[0 0 J1-1 I1-1],[1 ll nj ni]));    
% Follow ~Haine et al., 2015 definition of FWC
% The better way would be to interpolate
% between S values to find the exact depth
    Fwc = 0; 
    for k = 1:ll-1
      dz = DZ(k);
      ss1=squeeze(S(k,iBG));
      ss2=squeeze(S(k+1,iBG));
      ac = Acell(iBG);
      Ib1=find(ss1>Sref);
      if k==1, Isrf=Ib1; end; % exclude region with surf S>Sref
      ss1(Ib1)=nan;
      ss2(Ib1)=nan;
      Ib2=find(ss2>Sref); % CHeck next layer see if S > Sref there
 % Interpolate to find depth of Sref
      fwc=dz.*(Sref-ss1)/Sref; % m
% Take care about layers where S becomes S>Sref
% add FWC below layer interface from depth of Sref
% up to the bottom lyaer k, S=0.5(S(k)+Sref)
      if ~isempty(Ib2)
	zz1 = ZZ(k);
	zz2 = ZZ(k+1);
	dZz = abs(zz2-zz1);
	dS = ss2-ss1;
	zzSref = zz1+dZz./dS.*(Sref-ss1);
	dltZs  = abs(zz2-zzSref);
	Smn=0.5*(ss1+Sref); % for integration from 
	dltFwc = dltZs.*(Sref-Smn)/Sref;
	fwc(Ib2) = fwc(Ib2)+dltFwc(Ib2);
      end
      fwc(Ib1)=0;    
      fwc=fwc(:);
      Fwc=Fwc+nansum(fwc.*ac)*1e-3; % BG volume-integrated, km3
    end
    Fwc(Isrf)=0; 
%  Fwc(Ilnd)=nan; % exclude not needed regions
    fprintf('FWC calculation %6.2f min\n',toc/60);
    fprintf('====> %s: Beaufort Gyre, FWC=%8.2f km3\n',datestr(dnmb),Fwc); 
%keyboard
    cc=cc+1;
    FWC.Title         = 'FWC in BG from GLBb0.08 Analysis';
    FWC.Source        = pthbin;
    FWC.TM(cc,1)      = dnmb;
    FWC.Fwc_km3(cc,1) = Fwc;
    FWC.Sref          = Sref;
    FWC.BG_Area_km2   = sum(ac);

  end % day loop
  
% saving data
  if s_mat>0
    fmat = sprintf('%sFWC_BGvol_GLBb008_%i.mat',pthmat,yr);
    fprintf('Saving %s\n\n',fmat);
    save(fmat,'FWC');
  end;
      
end










