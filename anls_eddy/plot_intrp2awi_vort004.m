% Plot relative Vorticity/f of the bottom layer and SSH contours
% interpolated fields 2 awi grid
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

dnmb =  datenum(2006,1,10);
DV   = datevec(dnmb);
yr   = DV(1);
mo   = DV(2);
mday = DV(3);
iday = dnmb-datenum(yr,1,1)+1;

zz0  = -100;
expt = '011';
TV = '17DD';
rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 
dnmb = datenum(2005,8,21); % date to plot

sfig = 0;
%txtb = 'plot_deepU_ssh.m';
btxt = 'plot_interp2awi_vort004.m';

%pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.04/data_mat/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/data_eddy_AWI/';
pthout  = '/Net/tholia/ddmitry/hycom/ARCc0.04/data_awi_intrp/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthfig  = sprintf('/Net/mars/ddmitry/hycom/ARCc0.04/%s/fig_vort/',...
		  expt);

fintrp = sprintf('%sarc004_%s_archm2awi.%4.4i_%2.2i.mat',...
		       pthout,expt,yr,mo);
fprintf('Loading %s\n',fintrp);
load(fintrp);
%keyboard

% Get AWI grid:
fmat = sprintf('%sgrid_eddy_tracking.mat',pthmat);
fprintf('Loading %s\n',fmat);
AWI = load(fmat);
alat=AWI.lat;
alon=AWI.lon;
lt1=min(alat);
lt2=max(alat);
ln1=min(alon);
ln2=max(alon);

n=length(alon);
m=length(alat);

Fcor = zeros(m,n);
Fcor(:,1)=2*omg*sind(alat);
Fcor(1,:)=2*omg*sind(alat(1));
for ii=2:n
  x0=alon(ii-1);
  x1=alon(ii);
  y1=alat(1);
  Fcor(1,ii) = 2*omg*sind(y1);
  DX(1,ii) = distance_spheric_coord(y1,x0,y1,x1);
  for jj=2:m
    y0=alat(jj-1);
    y1=alat(jj);
    dst=distance_spheric_coord(y0,x1,y1,x1);
    DY(jj,ii) = dst;
    if ii==2
      DY(jj,1)=dst;
    end
    dst=distance_spheric_coord(y1,x0,y1,x1);
    DX(jj,ii)=dst;
    Fcor(jj,ii) = 2*omg*sind(y1);
  end
end
DX(:,1)=DX(:,2);
DY(1,:)=DY(2,:);


UN = squeeze(HYCOM.U(mday,:,:));
VN = squeeze(HYCOM.V(mday,:,:));
S  = sqrt(UN.^2+VN.^2);
%keyboard


% Calculate dv/dx
fprintf('Calculating dv/dx & du/dy\n');
dVdX=zeros(m,n)*nan;
dUdY=zeros(m,n)*nan;
%for j=1:m
%  dvp = squeeze(VN(j,2:end));
%  dv0 = squeeze(VN(j,1:end-1));
%  dx  = DX(j,1:end-1);
%  dv0=dv0(:)';
%  dvp=dvp(:)';
%  dx=dx(:)';
%  dVdX(j,:) = (dvp-dv0)./dx;
%end
%dVdX(:,n)=nan;
%for i=1:n
%  dup = squeeze(UN(2:end,i));
%  du0 = squeeze(UN(1:end-1,i));
%  dy  = DY(1:end-1,i);
%  dup=dup(:);
%  du0=du0(:);
%  dUdY(:,i) = (dup-du0)./dy;
%end
%dUdY(m,:)=nan;

for i=1:n-1
  for j=1:m-1
    dx=DX(j,i);
    dy=DY(j,i);
    dv=VN(j,i+1)-VN(j,i);
    dVdX(j,i)=dv/dx;
    du=UN(j+1,i)-UN(j,i);
    dUdY(j,i)=du/dy;
  end
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
%ZtF(HH>zz0)=nan;

% Arctic Ocean:
xl1=1;
xl2=800;
yl1=1;
yl2=825;


fprintf('Plotting %s\n',datestr(dnmb));
%figure('Visible','off'); clf;
figure(1); clf;
pcolor(ZtF); shading flat;
hold on;
%contour(HH,[0 0],'k');
caxis([-0.25 0.25]);
colormap(cmp);
axis('equal');
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[],...
	'ytick',[]);
set(gca,'Color',[0 0 0]);
ch = colorbar;
stl=sprintf('ARCc0.04-%s, AWI grid, Zt/f, depth=%5.1fm, %i/%2.2i/%2.2i',expt,zz0,DV(1:3));
title(stl);

bottom_text(btxt,'pwd',1,'fontsize',10);

if sfig>0
  fgnm=sprintf('%sarc004_%s_intrpAWI_vort_%i%2.2i%2.2i',pthfig,expt,DV(1:3));
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end
