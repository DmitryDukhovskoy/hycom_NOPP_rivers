% Plot depth and Tamx in Atl. layer
% Arctic Ocean
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps/
startup

clear all
close

yr   = 2015;
iday = 15; % check winter day to avoid warm surface layer

pfld  = 'temp';
f_extr = 1;  % =0 - load in extracted depth of Atl. Water, =1 -extract

s_fig=0;

rg = 9806;

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);
pthdat = sprintf('/Net/kronos/ddmitry/EN4_v4.2.1/%4.4i_anls/',yr);
pthmat  = '/nexsan/people/ddmitry/data_mat/';
pthmat2 = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat/';
%pthmat2= '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/fig_EN4/';
en4v   = 'EN.4.2.1.f.analysis.g10'; % EN4 version
txtb   = 'plot_atlw_EN4.m';

% Extract data for:
phi0 = 50;  % southmost latitude
zmin = -200; % depth isobath - within which perform the analysis


% Get topo:
ftopo = sprintf('%s/depth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn] = size(LON);

dnmb = datenum(yr,1,1)+iday-1;
dv   = datevec(dnmb);
im   = dv(2);
fnm  = sprintf('%s%s.%i%2.2i.nc',pthdat,en4v,yr,im);

fprintf('Plotting EN4 %s\n',datestr(dnmb));

LNe  = nc_varget(fnm,'lon');
LT0  = nc_varget(fnm,'lat');
iy   = max(find(LT0<=phi0))-1;
LTe  = LT0(iy+1:end);
n    = length(LNe);
m    = length(LTe);
TT   = squeeze(double(nc_varget(fnm,'temperature')))-273.15;
TT   = TT(:,iy+1:end,:);

ZM  = -1*nc_varget(fnm,'depth');
lz = length(ZM);
dzm = abs(diff(ZM));
dzm(lz) = dzm(lz-1);

clear ZZ
ZZ(1,1) = 0;
for kk=1:lz
  ZZ(kk+1) = -(abs(ZM(kk))+abs(0.5*dzm(kk)));
end
ZZ=ZZ(:);
DZ=abs(diff(ZZ)); % layer thicknesses
k    = length(ZZ);


% Create bottom map for EN4 from HYCOM ARCc0.08 topo:
fen4=sprintf('%sEN4_topo.mat',pthmat2);
f_extr = 0;
if f_extr
  Hen4 = sub_intrp_Hen4(LAT,LON,HH,LNe,LTe);
  fprintf('Saving EN4 topo %s\n',fen4);
  save(fen4,'Hen4');
else
  fprintf('Loading GDEM topo %s\n',fen4);
  load(fen4); % interpoalted for 
end


[LONe,LATe] = meshgrid(LNe,LTe);
Iarc = find(Hen4<=-100 & LATe>70);
nI = length(Iarc);
%[Y,X] = ind2sub(size(Hen4),Iarc);
%plot(X,Y,'r.');

[ny,nx] = size(Hen4);

%I = find(LONe>180);
%LONe(I)=LONe(I)-360;

fprintf('Seaching for atl. water depth, EN4, month%i ...\n',im);
Zmax = zeros(ny,nx);
Zt0  = Zmax;
Tmax = Zmax-999;

for i=1:nI
  if mod(i,1000)==0;
    fprintf('Processing, done %6.2f%%\n',i/nI*100);
  end;
  i0=Iarc(i);

  if (Hen4(i0)>zmin); continue; end;
  
  tt = squeeze(TT(:,i0));
  z0 = ZZ;
  iz0 = max(find(z0>-100));
  if ~isempty(iz0)
    tt(1:iz0) = -10;
  end
  izb = max(find(z0>=Hen4(i0)));
  tt(izb+1:end) = nan;
  
  tm = max(tt);
%    fprintf('tm=%6.4f\n',tm);
  if tm<0, continue; end;

  iz = find(tt==max(tt),1);
  Zmax(i0) = ZZ(iz);
  iz = min(find(tt>=0));

% Interpolate to find exact z
  iz1=iz-1;
  iz2=iz;
  t1 = tt(iz1);
  t2 = tt(iz2);
  z2 = ZZ(iz2);
  z1 = ZZ(iz1);
  t0 = 0;
  dtdz = (t2-t1)/(z2-z1);
  dz0 = (t0-t1)/(dtdz);
  zt0 = z1+dz0;
  Zt0(i0) = zt0;
  Tmax(i0)= tm;

end


xlim1 = 50;
xlim2 = nn-10;
ylim1 = 400;
ylim2 = mm-450;

% Interpolate GDEM into ARCc grid:
I=find(LNe>180);
if ~isempty(I); LNe(I)=LNe(I)-360; end;
ipiv=find(LNe==-179);
dmm=LNe;
LNe=dmm*0;
LNe(1:n-ipiv+1,1) = dmm(ipiv:end);
LNe(n-ipiv+2:n)   = dmm(1:ipiv-1);
dl = abs(dmm(2)-dmm(1));
LNe = [LNe(1)-dl; LNe; LNe(n)+dl];


% Plot 0-dgr depth
Zt1 = Zt0; 
A2=Zt0;
A2(:,1:n-ipiv+1,1) = Zt0(:,ipiv:end);
A2(:,n-ipiv+2:n)   = Zt0(:,1:ipiv-1);
A2 = [Zt0(:,ipiv-1), A2, Zt0(:,ipiv)];
A2(A2==0) = nan;
A=interp2(LNe,LTe,A2,LON,LAT);
A(Hen4>zmin)=nan;
A(A==0) = nan;
Zt0 = A;

stl = sprintf('EN4, Depth upper 0C m, B.Thick = 250m+/-50,%i/%i',im,yr);
txtb = 'plot_atlw_EN4.m';
pfld = 'AtlZ';
nf = 1;
sub_plot_scalar(Zt0,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld);
contour(Zt0,[-500:50:-150],'k');
contour(Zt0,[-250 -250],'k','linewidth',1.6)
bottom_text(txtb,'pwd',1,'Fontsize',8);

if s_fig>0
  fgnm = sprintf('%sGDEM_0dgr_depth_mo%2.2i',pthfig,im);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end


% Plot Max T in Atl. layer:
Tmx1 = Tmax;
A2=Tmax;
A2(:,1:n-ipiv+1,1) = Tmax(:,ipiv:end);
A2(:,n-ipiv+2:n)   = Tmax(:,1:ipiv-1);
A2 = [Tmax(:,ipiv-1), A2, Tmax(:,ipiv)];
A2(A2<-2) = nan;
A=interp2(LNe,LTe,A2,LON,LAT);
A(Hen4>zmin)=nan;
%A(A==0) = nan;
Tmax = A;

stl = sprintf('EN44, Tmax Atl. L., B.Thick = 1C+/-0.5C, %i/%i',im,yr);
pfld = 'AtlT';
nf = 2;
sub_plot_scalar(Tmax,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld);
contour(Tmax,[0:0.5:2],'k');
contour(Tmax,[1 1],'k','linewidth',1.6)
bottom_text(txtb,'pwd',1,'Fontsize',8);

if s_fig>0
  fgnm = sprintf('%sEN4_AtlL_Tmax_mo%2.2i',pthfig,im);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end

