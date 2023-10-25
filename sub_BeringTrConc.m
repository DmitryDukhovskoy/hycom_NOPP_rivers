function cBr = sub_BeringTrConc(f_drv);
% This approach is not accurate 
%
%
% estimate what tracer conc
% should be in the Bering Strait
% to make it comparable to other tracers
% Tracer = m3/sec of river runoff
%
% All other FW sources, Tr. conc = FW Flux in this grid point (m3/s)
% Tracers are distributed in the upper 2 layers dz = 1 m each
% It is assumed that after dt = relaxation time scale (1 day)
% Tr. concentration in the grid cell = FWFlux, i.e.
%
% In the simulation Bering Tracer =1
% Sum (all grid points) = N
% N => 0.088e6 m3/sec of FW flux through Bering (model estimate)
% coeff= 0.088e6/N
%
% Need to know the overall vol of the region where
% tracers are seeded (Nseed grid points)
% Assume FWFlx (Sref=34.8) ~=3000 km3/yr or in the model it is ~2770 km3/yr)
% Woodgate, 2012: 2000-2500 in 2001, 3000-3500 km3/yr in 2011
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
rg = 9806; 

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

if isempty(f_drv)
  f_drv = 1; % calculate coefficient
           % =0 - use saved 
end

	   
%cBr = 11.78; % previously estimated value
cBr = 3.49; % previously estimated value

if f_drv==0; return; end;

% Bering FWFlux:
% 1e5 m3/s = 0.1 Sv ~=3153 km3/yr
% 0.88e5 m3/s = 0.088 Sv = 2779 km3/yr - Bering FWFlux in the simulation
%

regn = 'ARCc0.08';
expt = 110;

% Bering Str. seeding points:
i1=639;
i2=656;
j1=1919;
j2=1925;


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

yr = 1999;
pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)
iday = 120;

fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

fld='thknss';
[F,n,m,l] = read_hycom(fina,finb,fld);
F(F>1e10)=nan;
F(F<0.1)=nan;
F=F./rg;
dH = F;

%[ZZ,ZM] = sub_thck2dpth(dH);

% Calculate volume:
Vlm=0;
cnc=0;
Vz1m=0;
AreaVrt=0;
for ii=i1:i2
  for jj=j1:j2
    ar=Acell(jj,ii);
    Asum=Asum+ar;
    dz=squeeze(dH(:,jj,ii));
    dz(dz<1e-6)=nan;
    Inn=find(~isnan(dz));
    d1z=dz;
    d1z(Inn)=1;
    vlm=nansum(dz*ar);
    Vlm=Vlm+vlm;
    Vz1m=Vz1m+nansum(d1z*ar);
    if jj==j1
      AreaVrt=AreaVrt+nansum(DX(jj,ii)*dz);
    end
    
    cnc=cnc+length(Inn); % all active points within the region
  end
end

% Estimate average FWF per 1 grid cell
% 1e5 m3/s = 0.1 Sv ~=3153 km3/yr
% 0.88e5 m3/s = 0.088 Sv = 2779 km3/yr - Bering FWFlux in the simulation
% At instant (1sec) - there is Vfw*1sec m3 of freshwater in 
% N m (width) x 2D section across the strait
% Given mean Velocity Nm = 0.27 m
%
% This volume of freshwater is redistributed over
% Vlm volume of tracers = 1 kg/m3 concentration
%
% Thus need to normalize by the area of the Bering str section
%
%rho_tr = 148; % kg/m3 - estimated Tr. conc for other FW sources
%fwb = 1e5; % m3/s
fwb = 0.88e5; % m3/s
Mfw = fwb*1000; % FW Mass in 1 sec: slub Nm width X Area across Bering
%Mtr = 1*Vlm;
%Mtr = 1*Vz1m; % mass tracer for 1 m dz in all cells
Mtr = AreaVrt;  % mass tracer for 1 m width slab across the strait

cBr = Mfw/Mtr*1/(j2-j1+1); % distributed over nj rows

return