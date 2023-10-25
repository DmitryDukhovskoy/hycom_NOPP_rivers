% Plot relative Vorticity/f of the bottom layer and SSH contours
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = '110';
TV = 11;
rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 
dnmb = datenum(2016,9,10); % date to plot

sfig = 0;
%zz0   = -5;  %  depth of calculation
zz0 = -10;
%txtb = 'plot_deepU_ssh.m';
btxt = 'plot_vort008.m';

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/fig_2D/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Fcor = 2*omg*sind(LAT);
Lmsk = HH*0;
Lmsk(HH<zz0) = 1;

DV =datevec(dnmb);
yr = DV(1);
iday = dnmb-datenum(yr,1,1)+1;

pthbin = sprintf('/nexsan/archive/ARCc0.08_%s/data/%4.4i/',expt,yr);

fina = sprintf('%s%s_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%s_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

if ~exist(fina,'file') | ~exist(finb,'file');
  fprintf('Not found *a or *b: %s\n',fina);
  fprintf('                     %s\n',finb);
%  continue;
end

%cnc=cnc+1;
%TM(cnc,1)=dnmb;

fprintf('Getting data expt %s: %4.4i_%2.2i_%2.2i: %s\n',expt,DV(1:3),fina);

[ZM,ZZ] = sub_zz_zm(fina,finb,HH);
ZZ(isnan(ZZ))=100;
ZM(isnan(ZM))=100;

nlr = size(ZM,1);
Iz=find(HH<zz0);
for k=1:nlr
  zav(k) = nanmean(ZM(k,Iz));
end
dzav = abs(zav-zz0);
Iz0 = find(dzav==min(dzav));

[F,n,m,nlr] = read_hycom(fina,finb,'u-vel.','r_layer',Iz0);
%  F=F(:,jnc1:jnc2,inc1:inc2);
F(F>huge)=0;
UN=squeeze(F);


[F,n,m,nlr] = read_hycom(fina,finb,'v-vel.','r_layer',Iz0);
F(F>huge)=0;
VN=squeeze(F);

% Calculate dv/dx
% central differencing to have du and dv collocated
fprintf('Calculating dv/dx & du/dy\n');
dVdX = HH*nan;
dUdY = HH*nan;
for j=1:m
  dvp1 = squeeze(VN(j,3:end));
  dvm1 = squeeze(VN(j,1:end-2));
  dvp1(isnan(dvp1))=0;
  dvm1(isnan(dvm1))=0;
  dx  = DX(j,2:end-1);
  dvm1=dvm1(:)';
  dvp1=dvp1(:)';
  dx=dx(:)';
  dVdX(j,2:nn-1) = (dvp1-dvm1)./(2*dx);
end
for i=1:n
  dup1 = squeeze(UN(3:end,i));
  dum1 = squeeze(UN(1:end-2,i));
  dup1(isnan(dup1))=0;
  dum1(isnan(dum1))=0;
  dy  = DY(2:end-1,i);
  dup1=dup1(:);
  dum1=dum1(:);
  dUdY(2:mm-1,i) = (dup1-dum1)./(2*dy);
end



cl1 = colormap_blue(100);
cl2 = colormap_orange(100);
cl2(1,:) = [1 1 1];
cl2(2,:) = [1 1 1];
cl2(3,:) = [1 1 1];
cl2(4,:) = [1 1 1];
cl2(5,:) = [1 1 1];
cl1(1,:) = [1 1 1];
cl1(2,:) = [1 1 1];
cl1(3,:) = [1 1 1];
cl1(4,:) = [1 1 1];
cl1(5,:) = [1 1 1];
cl1 = flipud(cl1);
cmp=[cl1;cl2];
cmp =smooth_colormap(cmp,5);
nint=length(cmp);

ZtF = (dVdX-dUdY)./Fcor;
ZtF(HH>zz0)=nan;

xl1=[];
xl2=[];
yl1=[];
yl2=[];


% Greenland:
%xl1=350;
%xl2=1250;
%yl1=50;
%yl2=1150;

% Arctic Ocean:
xl1=200;
xl2=1580;
yl1=600;
yl2=1900;

% Lofoten:
% 66-73 N latitude
% 10W - 18E long.
lft=[];
if ~isempty(lft)
  d1=distance_spheric_coord(66,-10,LAT,LON);
  [j1,i1]=find(d1==min(min(d1)));
  d1=distance_spheric_coord(73,-10,LAT,LON);
  [j2,i2]=find(d1==min(min(d1)));
  d1=distance_spheric_coord(73,18,LAT,LON);
  [j3,i3]=find(d1==min(min(d1)));
  d1=distance_spheric_coord(66,18,LAT,LON);
  [j4,i4]=find(d1==min(min(d1)));
  xl1=950;
  xl2=1250;
  yl1=500;
  yl2=820;
end

% Fram Strait:
% find AWI grid boundary:
%lt1=75;
%lt2=82.5;
%ln1=-20;
%ln2=20;

%Nordic Seas:
lt1=60;
lt2=81.;
ln1=-35;
ln2=15;
%ln2=30;

% SPNA
%lt1=50;
%lt2=80;
%ln1=-80;
%ln2=15;
%xl1=270;
%xl2=1310;
%yl1=25;
%yl2=1100;

if isempty(xl1)
  d1=distance_spheric_coord(lt1,ln1,LAT,LON);
  [j1,i1] = find(d1==min(min(d1)),1);
  d1=distance_spheric_coord(lt2,ln1,LAT,LON);
  [j2,i2] = find(d1==min(min(d1)),1);
  d1=distance_spheric_coord(lt2,ln2,LAT,LON);
  [j3,i3] = find(d1==min(min(d1)),1);
  d1=distance_spheric_coord(lt1,ln2,LAT,LON);
  [j4,i4] = find(d1==min(min(d1)),1);
  jL1=min([j1,j2,j3,j4]);
  jL2=max([j1,j2,j3,j4]);
  %jL2=2200;
  iL1=min([i1,i2,i3,i4]);
  iL2=max([i1,i2,i3,i4]);
  %iL2=1145;
  xl1 = iL1;
  xl2 = iL2;
  yl1 = jL1;
  yl2 = jL2;
end

%iV = [iL1,iL1,iL2,iL2];
%jV = [jL1,jL2,jL2,jL1];
iV = [xl1,xl1,xl2,xl2];
jV = [yl1,yl2,yl2,yl1];
[II,JJ]=meshgrid([1:nn],[1:mm]);
%dmm = inpolygon(II,JJ,iV,jV);
%IN  = find(dmm==1);
%IOUT= find(dmm==0);

%HH(IOUT)=nan;
%ZtF(IOUT)=nan;

ZtF(HH>zz0)=nan;


figure(1); clf;
set(gcf,'Position',[1262 431 931 900]);
axes('Position',[0.05 0.1 0.8 0.8]);
pcolor(ZtF); shading flat;
hold on;
contour(HH,[0 0],'w');
%contour(HH,[-5000:500:-200],'Color',[0.7 0.7 0.7]);
caxis([-0.25 0.25]);
colormap(cmp);
%plot([i1 i2],[j1 j2],'m');
%plot([i2 i3],[j2 j3],'m');
%plot([i3 i4],[j3 j4],'m');
%plot([i1 i4],[j1 j4],'m');

axis('equal');
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[],...
	'ytick',[],...
	'Color',[0 0 0]);
ch = colorbar;
set(ch,'Position',[0.85 0.11 0.02 0.8],...
       'TickLength',0.02,...
       'Fontsize',16);
stl=sprintf('ARCc0.08-%s, Zt/f, depth=%5.1fm, %i/%2.2i/%2.2i',expt,zz0,DV(1:3));
title(stl);

% set(gca,'xlim',[900 1600],'ylim',[750 1400]);

bottom_text(btxt,'pwd',1,'fontsize',10);

if sfig>0
  fgnm=sprintf('%sarc008_%s_vort_%i%2.2i%2.2i',pthfig,expt,DV(1:3));
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end
