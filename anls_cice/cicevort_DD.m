addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
%addpath '/home/bhaynes/matlabscripts';
startup; 

clear
close all

sfig = 1;
expt='110';
pthfig = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/frames_cice/';
%pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthbin  = '/nexsan/people/ddmitry/hycom/ARCc0.08/110/cice/'; 

%ocn=[0.0040,0.7296,0.7701];
%land=[0.8057,0.6772,0.5620];
%cmap=[ocn;bluewhitered;land];
%close(gcf);

%v=VideoWriter('cicevort_1994_03.avi');
%v.FrameRate=3;
%v.Quality=95;
%open(v);

pthtopo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
ftopo=sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH      = nc_varget(ftopo,'Bathymetry');
LON     = nc_varget(ftopo,'Longitude');
LAT     = nc_varget(ftopo,'Latitude');
[mm,nn] = size(HH);
om      = 0.000072921;
Fcor    = 2*om.*sind(LAT);

[DX,DY]      = sub_dx_dy(LON,LAT);
Lmsk         = HH*0;
Lmsk(HH<0) = 1;


YRPLT=[];
cc=0;
for iyr=2005:2005
  for idd=1:366
%    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

%cl1 = colormap_blue(100);
%cl2 = colormap_orange(100);
cl1 = colormap_cold(100);
%cl2 = colormap_red(100);
cl2 = colormap_hot(100);
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

lmp = [0 0 0; 0.8 0.8 0.8];

fcnt = 0;
for ip=1:np
  tic;
  yr   = YRPLT(ip,1);
  iday = YRPLT(ip,2);
  dnmb = datenum(yr,1,1)+iday-1;
  dv   = datevec(dnmb);
  mo   = dv(2);
  mday = dv(3);
  
  pthbin = sprintf('/nexsan/archive/ARCc0.08_%s/data/%i_cice/',expt,yr);
  fin    = sprintf('%s%s_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',...
      pthbin,expt,yr,mo,mday);
  
  if ~exist(fin,'file');
    fprintf('Not found %s\n',fin);
    continue;
  end

  ui = squeeze(nc_varget(fin,'uvel'));
  vi = squeeze(nc_varget(fin,'vvel'));
  hi = squeeze(nc_varget(fin,'hi'));
  ci = squeeze(nc_varget(fin,'aice'));

  fprintf('Calculating dv/dx & du/dy\n');
  dVdX = HH*nan;
  dUdY = HH*nan;
  for j=1:mm
    dvp = vi(j,2:end);
    dv0 = vi(j,1:end-1);
    dx  = DX(j,1:end-1);
    dvp = dvp(:)';
    dv0 = dv0(:)';
    dx  = dx(:)';
    dVdX(j,1:end-1) = (dvp-dv0)./dx;
  end

  for i=1:nn
    dup = ui(2:end,i);
    du0 = ui(1:end-1,i);
    dy  = DY(1:end-1,i);
    dup = dup(:);
    du0 = du0(:);
    dUdY(1:end-1,i) = (dup-du0)./dy;
  end
  
  ZtF = (dVdX-dUdY)./Fcor; %k-component of vorticity = zetta
  ZtF(hi<1e-2 | ci<0.15)=nan;
  ZtF(HH>=0)=nan;
  

  figure(1); clf;
  axes('position',[0.06 0.08 0.84 0.84]);
  pcolor(Lmsk); shading flat;
  colormap(lmp);
  freezeColors;
  hold on;
  
  pcolor(ZtF); shading flat;
  colormap(cmp);
  caxis([-0.1 0.1]);
  
  axis('equal');
  set(gca,'xlim',[600 1200],...
	  'ylim',[500 1100]);
  set(gca,'xtick',[],'ytick',[]);
  stl = sprintf('ARCc0.08_%s, sea ice vort/f, %4.4i/%2.2i/%2.2i',...
		expt,dv(1:3));
%  clr=[0.9 0.9 0.9];
%  plot_gridlines(45,10,1,clr,LON,LAT);
  title(stl,'Fontsize',12,'Interpreter','none');
  
  hb=colorbar;
  set(hb,'position',[0.89 0.1 0.02 0.8],...
	 'TickLength',0.025,...
	 'FontSize',12);
  
  fcnt = fcnt+1;
  if sfig
    fnm=sprintf('arc08_%s_icevort_%4.4i',expt,fcnt);
    fg=[pthfig,fnm];
    fprintf('Saving figure %s\n',fg);
    print('-dpng','-r300',fg);
  end
  
  fprintf('1 record processed: %8.6f min\n\n',toc/60);

end
%close(v);

