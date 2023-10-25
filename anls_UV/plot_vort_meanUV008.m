% Plot relative Vorticity/f of the bottom layer and SSH contours
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

plr=1;
expt = 110;
sexpt = sprintf('%3.3i',expt);
TV = 11;
rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 
dnmb = datenum(2005,8,21); % date to plot

sfig = 0;
%zz0   = -5;  %  depth of calculation
%txtb = 'plot_deepU_ssh.m';
btxt = 'plot_vort_meanUV008.m';

%pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthfig  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/fig_2D/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);

YRS = 2005; % more than 1 year - mean over these years
np = length(YRS);
imo   = 12; % >12 - annual mean
if np>1, imo=13; end;

%if imo>13
%  fprintf('Streamlines for annual mean %i, mat saved 0/1: %i\n',YRS,f_mat);
%end
%if imo<13
%  fprintf('Streamlines for month %i/%i, mat saved 0/1: %i\n',imo,YRS,f_mat);
%end




ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Fcor = 2*omg*sind(LAT);
Lmsk = HH*0;
Lmsk(HH<0)=1;
%Lmsk(HH<zz0) = 1;

hmsk=HH;
hmsk(HH<0)=nan;

[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);
xlim1 = 20;
xlim2 = nn-1;
ylim1 = 100;
ylim2 = mm-100;

usm=zeros(mm,nn);
vsm=zeros(mm,nn);
cc = 0;
for ik=1:np
  iyr = YRS(ik);
  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  if imo>12
    for ik=1:12
      cc = cc+1;
      U = meanUV(ik).U;
      V = meanUV(ik).V;
      usm = usm+U;
      vsm = vsm+V;
    end
  else
    cc=1;
    usm = meanUV(ik).U;
    vsm = meanUV(ik).V;
  end
end;

U = usm./cc;
V = vsm./cc;

S = sqrt(U.^2+V.^2);

xlim1 = 400;
xlim2 = 1250;
ylim1 = 130;
ylim2 = 1200;


% Calculate dv/dx
fprintf('Calculating dv/dx & du/dy\n');
for j=1:mm
  dvp = squeeze(V(j,2:end));
  dv0 = squeeze(V(j,1:end-1));
  dx  = DX(j,1:end-1);
  dv0=dv0(:)';
  dvp=dvp(:)';
  dx=dx(:)';
  dVdX(j,:) = (dvp-dv0)./dx;
end
dVdX(:,nn)=nan;
for i=1:nn
  dup = squeeze(U(2:end,i));
  du0 = squeeze(U(1:end-1,i));
  dy  = DY(1:end-1,i);
  dup=dup(:);
  du0=du0(:);
  dUdY(:,i) = (dup-du0)./dy;
end
dUdY(mm,:)=nan;

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
%ZtF(HH>zz0)=nan;
ZtF(xlim2:end,:)=nan;
ZtF(:,ylim2:end)=nan;


figure(1); clf;
pcolor(ZtF); shading flat;
hold on;
contour(HH,[0 0],'k');
contour(HH,[-5000:1000:-200],'Color',[0.7 0.7 0.7]);
caxis([-0.1 0.1]);
colormap(cmp);

axis('equal');
set(gca,'xlim',[xlim1 xlim2],...
	'ylim',[ylim1 ylim2],...
	'xtick',[],...
	'ytick',[]);
ch = colorbar;
if imo>12
  stl=sprintf('ARCc0.08-%s, Zt/f, meanUV, upper 50m, %i',sexpt,YRS);
else
  stl=sprintf('ARCc0.08-%s, Zt/f, meanUV, upper 50m, %2.2i/%i',sexpt,imo,YRS);
end

title(stl);

bottom_text(btxt,'pwd',1,'fontsize',10);

if sfig>0
  fgnm=sprintf('%sarc008_%s_vort_%i%2.2i%2.2i',pthfig,expt,DV(1:3));
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end
