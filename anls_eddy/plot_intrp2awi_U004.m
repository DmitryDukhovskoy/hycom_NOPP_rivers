% Plot U at zz0 depth from HYCOM 0.04
% interpolated into AWI grid
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
btxt = 'plot_interp2awi_U004.m';

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

CMP = create_colormap3(200,-0.3, 0.3);
cmpS=CMP.colormap;
for k=1:10
  cmpS(k,:)=[1 1 1];
%  cl2(k,:)=[1,1,1];
end
cmpS = smooth_colormap(cmpS,15);

%figure('Visible','off');
figure(1); clf;
pcolor(S); shading flat;
colormap(cmpS);
caxis([0 0.5]);
set(gca,'Color',[0 0 0]);
stl=sprintf('|U|, HYCOM 0.04 AWI grid %i/%2.2i/%2.2i',DV(1:3));

% Plot vectors
fprintf('Plotting vectors ...\n');
[mm,nn]=size(UN);

scl=9;
cf=0.3;
beta=20;
lwd=1.2;
v_col=[0 0 0];
dii=15;
for ii=10:dii:nn
  for jj=10:dii:mm
    clear u v
    u = UN(jj,ii);
    v = VN(jj,ii);
    s = S(jj,ii);
    if isnan(u), continue; end;
  %    if res>0,
    u=u/s;
    v=v/s;
    scl=25;
  %    end

    x0=ii;
    y0=jj;

    x1=x0+u*scl;
    y1=y0+v*scl;
    draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
  end
end

%txtb = 'plot_meanUV.m';
bottom_text(btxt,'pwd',1,'fontsize',8);


title(stl);
colorbar;
if sfig>0
  fgnmS=sprintf('%shycom004_intrpAWI_U_%4.4i%2.2i%2.2i',pthfig,DV(1:3));
  fprintf('Saving %s\n',fgnmS);
  print('-dpng','-r250',fgnmS);
end

