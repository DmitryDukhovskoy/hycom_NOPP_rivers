function [YRS,dP,mnP] = sub_prcp_NCEPR;
% NCEP/DOE AMIP-II Reanalysis (Reanalysis-2) Daily Averages Data
% For the N. Atl region
pthdat = '/Net/mars/ddmitry/hycom/NCEPR2_precip/';
fprc   = sprintf('%sNAtl_precip_rate_NCEPR2.nc',pthdat);

T  = nc_varget(fprc,'time')/24;
TM = datenum(1800,1,1)+T;
DV = datevec(TM);
X  = nc_varget(fprc,'lon');
Y  = nc_varget(fprc,'lat');
Pr = nc_varget(fprc,'prate'); % kg/m^2/s

I=find(X>180);
X(I)=X(I)-360;
dX = diff(X);
dY = abs(diff(Y));
[XX,YY] = meshgrid(X,Y);
[m,n] = size(XX);
for i=1:n
  for j=1:m
    x0=X(i);
    y0=Y(j);
    dx=dX(1);
    dy=dY(1);
    x1=x0+dx;
    y1=y0+dy;
    llx=distance_spheric_coord(y0,x0,y0,x1);
    lly=distance_spheric_coord(y0,x0,y1,x0);
    DX(j,i)=llx;
    DY(j,i)=lly;
  end
end
Acell = DX.*DY; % m2

% area-integrated
rho=1000; 
nt=length(TM);
for it=1:nt
  dmm=squeeze(Pr(it,:,:));
  mass = dmm.*Acell*24*3600;
  vol  = mass/rho; % m3/day
  Vfw(it,1)  = nansum(nansum(vol)); %m3/day
end

% Annual precip flux, m3/year
cc=0;
yr1=1990;
yr2=2016;
clear VFWy
for iyr=yr1:yr2
  I=find(DV(:,1)==iyr);
  dmm=sum(Vfw(I));
  cc=cc+1;
  VFWy(cc,1)=dmm*1e-9;  % km3/yr
end

YRS = [yr1:yr2];
%mnP = mean(VFWy);
mnP=mean(VFWy(1:19)); %mean 1990-2009
dP = VFWy-mnP; % FW flux anomaly

return



