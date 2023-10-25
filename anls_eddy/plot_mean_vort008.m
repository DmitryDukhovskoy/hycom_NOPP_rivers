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

sfig = 0;
zz0   = -100;  %  depth of calculation
%txtb = 'plot_deepU_ssh.m';
btxt = 'plot_mean_vort008.m';

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/fig_2D/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

yr = 2006;
fmat = sprintf('%smean_vort008_%i.mat',pthmat,yr);

fprintf('Loading %s\n',fmat);
load(fmat);

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


% Lofoten:
% 66-73 N latitude
% 10W - 18E long.
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

ZtF = MVRT.Vrt_Fcor;
zz0 = MVRT.Depth;

figure(1); clf;
pcolor(ZtF); shading flat;
hold on;
contour(HH,[0 0],'k');
contour(HH,[-5000:1000:-200],'Color',[0.7 0.7 0.7]);
caxis([-0.25 0.25]);
colormap(cmp);
plot([i1 i2],[j1 j2],'m');
plot([i2 i3],[j2 j3],'m');
plot([i3 i4],[j3 j4],'m');
plot([i1 i4],[j1 j4],'m');

axis('equal');
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[],...
	'ytick',[]);
ch = colorbar;
stl=sprintf('ARCc0.08-%s, mean Zt/f, depth=%5.1fm, %i',expt,zz0,yr);
title(stl);

bottom_text(btxt,'pwd',1,'fontsize',10);

% Plot mean U,V
U = MVRT.U_mean;
V = MVRT.V_mean;
U(HH>0) = nan;
V(HH>0) = nan;
S = sqrt(U.^2+V.^2);

CMP = create_colormap_WBYR(200,-0.3, 0.3);
cmpS=CMP.colormap;
%for k=1:10
%  cmpS(k,:)=[1 1 1];
%  cl2(k,:)=[1,1,1];
%end
%cmpS = smooth_colormap(cmpS,15);


%figure('Visible','off');
figure(2); clf;
pcolor(S); shading flat;
colormap(cmpS);
hold on
caxis([0 0.3]);
contour(HH,[-5000:1000:-100],'Color',[1 1 0.3]);
axis('equal');
%set(gca,'xlim',[1000 1200],...
%	'ylim',[600 750],...
%	'color',[0 0 0]);
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'color',[0 0 0]);
set(gca,'xtick',[],...
	'ytick',[]);
% Plot vectors
fprintf('Plotting vectors ...\n');
[mm,nn]=size(U);

% Plot domain from Volkov et al., 2015
plot([i1 i2],[j1 j2],'m');
plot([i2 i3],[j2 j3],'m');
plot([i3 i4],[j3 j4],'m');
plot([i1 i4],[j1 j4],'m');

scl=3;
cf=0.3;
beta=18;
lwd=1.2;
v_col=[0 0 0];
dii=2;
for ii=xl1:dii:xl2
  for jj=yl1:dii:yl2
    clear u v
    u = U(jj,ii);
    v = V(jj,ii);
    s = S(jj,ii);
    if isnan(u), continue; end;
  %    if res>0,
    u=u/s;
    v=v/s;
%    scl=25;
  %    end

    x0=ii;
    y0=jj;

    x1=x0+u*scl;
    y1=y0+v*scl;
    draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
  end
end
colorbar;

stl=sprintf('ARCc0.08-%s, mean|U|, depth=%5.1fm, %i',expt,zz0,yr);
title(stl);

bottom_text(btxt,'pwd',1,'fontsize',10);




