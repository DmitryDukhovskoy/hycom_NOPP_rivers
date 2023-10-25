% Plot monthly mean UV fields
% derived in mnthly_arc04_UV.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;
f_mnth = 9; % plot <12 - month, >12 - annual mean

plr =10;  % U from plr layer ~50 m
rg = 9806;


regn = 'ARCc0.04';
expt = 022;
pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.04/011/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/data_mat/',expt);
pthmat  = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/u_mnth/';


%YRPLT = [2011,2012,2013,2014,2015];
YRPLT = [2016];
np = length(YRPLT);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);
% North Atlantic
%xlim1 = 100;
%xlim2 = 2800;
%ylim1 = 5;
%ylim2 = 2200;

xlim1 = 100;
xlim2 = 3200;
ylim1 = 10;
ylim2 = 5040;

for ik=1:np
  iyr = YRPLT(ik);
%  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
    
  if f_mnth>12
    usm=zeros(mm,nn);
    vsm=zeros(mm,nn);
    cc = 0;
    for imo=1:12
      cc = cc+1;
      U = meanUV(imo).U;
      V = meanUV(imo).V;
      usm = usm+U;
      vsm = vsm+V;
    end
    U = usm./cc;
    V = vsm./cc;
  else
    imo=f_mnth;
    U = meanUV(imo).U;
    V = meanUV(imo).V;
  end

  S = sqrt(U.^2+V.^2);
  nf = 1;
  if f_mnth>12
    stl = sprintf('ARCc0.04-%i, Mean U Lr %i, %i',expt,plr,iyr);
  else
    stl = sprintf('ARCc0.04-%i, Mean U Lr %i, %i/%2.2i',expt,plr,iyr,f_mnth);
  end
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
  pfld='speed';
  hps = [0.91 0.1 0.015 0.8];
%keyboard
%  sub_plot_scalar(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
%		  'c1',0,'c2',0.5,'clbpos',hps);
%  contour(HH,[-6000:1000:0],'Color',[0.6 0.6 0.6]);
  Lmsk = HH*0;
  Lmsk(HH<0)=1;
  lcmp=[0.5 0.5 0.5; 0.9 0.9 0.9];

  c1=0;
  c2=0.5;
%  CMP = create_colormapBGY(200,c1,c2);
%  cmp = CMP.colormap;
  CMP = create_colormap_WBYR(200,c1,c2);
  cmp = CMP.colormap;


		figure(nf); clf;
		axes('Position',[0.05 0.1 0.8 0.8]);
		hold on;
		pcolor(Lmsk); shading flat;
		colormap(lcmp);
		freezeColors;

		pcolor(S); shading flat;
		caxis([c1 c2]);
		colormap(cmp);

		axis('equal');

		set(gca,'xlim',[xlim1 xlim2],...
										'ylim',[ylim1 ylim2],...
										'xtick',[],...
										'ytick',[]);
		%        'Color',[0 0 0]);

		title(stl,'Fontsize',12,'Interpreter','none');
				
		hb=colorbar;
		set(hb,'Position',[0.87 0.11 0.02 0.8],...
									'Ticklength',0.03,...
									'Fontsize',14);




  i0=10;
  j0=10;
%  Psi = stream_fn(U,V,DX,DY,i0,j0,'simps');
%keyboard

% Greenland, scaling
%  scl=16;
%  cf=0.3;
%  beta=20;
%  lwd=1.;
%  v_col=[0.6 0.6 0.6];
%  dii=15;
% SE Greenland
  scl=150;
  cf=0.3;
  beta=15;
  lwd=1.;
  v_col=[0. 0. 0.];
  dii=40;
  for ii=xlim1:dii:xlim2
    for jj=ylim1:dii:ylim2
      clear u v
      u = U(jj,ii);
      v = V(jj,ii);
      s = S(jj,ii);
      if isnan(u) | s<0.02, continue; end;
  %    if res>0,
%	u=u/s;
%	v=v/s;
%	scl=25;
  %    end

      x0=ii;
      y0=jj;

      x1=x0+u*scl;
      y1=y0+v*scl;
      draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
    end
  end

  txtb = 'plot_arc04_meanUV.m';
  bottom_text(txtb,'pwd',1);
  
end



