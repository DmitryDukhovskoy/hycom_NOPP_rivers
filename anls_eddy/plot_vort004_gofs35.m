% Plot relative Vorticity/f of the bottom layer and SSH contours
% Simulation from the HYCOM-CICE5 GOFS3.5 
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = '022';
TV = '17DD';
rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 
dnmb = datenum(2017,5,30); % date to plot

sfig = 0;
zz0   = -10;  %  depth of calculation
%txtb = 'plot_deepU_ssh.m';
btxt = 'plot_vort004_gofs35.m';

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.04/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%s/fig_2D/',...
		  expt);
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);

ftopo = sprintf('%sdepth_ARCc0.04_%s.nc',pthtopo,TV); % 
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


pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i_mean/',yr);

% 
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

%fld='thknss';
%[F,n,m,nlr] = read_hycom(fina,finb,fld);
%F(F>1e10)=nan;
%F(F<0.01)=0;
%dH=F/rg;

[ZM,ZZ] = sub_zz_zm(fina,finb,HH);
ZZ(isnan(ZZ))=100;
ZM(isnan(ZM))=100;

[nlr,a1,a2]=size(ZM);

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
fprintf('Calculating dv/dx & du/dy\n');
for j=1:m
  dvp = squeeze(VN(j,2:end));
  dv0 = squeeze(VN(j,1:end-1));
  dx  = DX(j,1:end-1);
  dv0=dv0(:)';
  dvp=dvp(:)';
  dx=dx(:)';
  dVdX(j,:) = (dvp-dv0)./dx;
end
dVdX(:,n)=nan;
for i=1:n
  dup = squeeze(UN(2:end,i));
  du0 = squeeze(UN(1:end-1,i));
  dy  = DY(1:end-1,i);
  dup=dup(:);
  du0=du0(:);
  dUdY(:,i) = (dup-du0)./dy;
end
dUdY(m,:)=nan;

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
xl1=700;
xl2=2500;
yl1=100;
yl2=2300;
% Arctic Ocean:
xl1=400;
xl2=3160;
yl1=1200;
yl2=3800;
% Fram Strait:
% find AWI grid boundary:
%lt1=75;
%lt2=82.5;
%ln1=-20;
%ln2=20;

%Nordic Seas:
%lt1=60;
%lt2=81.;
%ln1=-35;
%ln2=15;

% SPNA
%lt1=50;
%lt2=80;
%ln1=-80;
%ln2=15;
%xl1=540;
%xl2=2620;
%yl1=50;
%yl2=2200;

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
  %iL2=2290;
  xl1 = iL1;
  xl2 = iL2;
  yl1 = jL1;
  yl2 = jL2;
end

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
caxis([-0.25 0.25]);
colormap(cmp);
axis('equal');
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[],...
	'ytick',[],...
	'Color',[0 0 0]);
ch = colorbar;
set(ch,'Position',[0.87 0.11 0.02 0.8],...
       'TickLength',0.02,...
       'Fontsize',16);

stl=sprintf('ARCc0.04-%s, Zt/f, depth=%5.1fm, %i/%2.2i/%2.2i',expt,zz0,DV(1:3));
title(stl);

bottom_text(btxt,'pwd',1,'fontsize',10);

if sfig>0
  fgnm=sprintf('%sarc004_%s_vort_%i%2.2i%2.2i',pthfig,expt,DV(1:3));
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end
