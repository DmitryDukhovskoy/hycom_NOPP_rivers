% Plot monthly mean UV fields
% derived in mnthly_arc08_UV.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;
f_mnth = 13; % plot <12 - month, >12 - annual mean

plr =1;  % U from plr layer
rg = 9806;


regn = 'ARCc0.04';
expt = 011;
pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.04/011/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/data_mat/',expt);

%YRPLT = [2011,2012,2013,2014,2015];
YRPLT = [2005];
np = length(YRPLT);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);
% Greenland
%xlim1 = 900;
%xlim2 = 2200;
%ylim1 = 600;
%ylim2 = 2200;
% SE Greenland:
xlim1 = 1200;
xlim2 = 1720;
ylim1 = 700;
ylim2 = 1280;

for ik=1:np
  iyr = YRPLT(ik);
  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  
  if f_mnth>12
    usm=zeros(mm,nn);
    vsm=zeros(mm,nn);
    cc = 0;
    for ik=1:12
      cc = cc+1;
      U = meanUV(ik).U;
      V = meanUV(ik).V;
      usm = usm+U;
      vsm = vsm+V;
    end
    U = usm./cc;
    V = vsm./cc;
  end

  S = sqrt(U.^2+V.^2);
  nf = 1;
  stl = sprintf('ARCc0.04-%i, Mean U %i',expt,iyr);
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
  pfld='speed';
  hps = [0.95 0.1 0.035 0.8];
  sub_plot_scalar(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
		  'c1',0,'c2',0.5,'clbpos',hps);
%  contour(HH,[-500:100:-5],'Color',[0.9 0.9 0.9]);

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
  scl=12;
  cf=0.3;
  beta=15;
  lwd=1.;
  v_col=[0.2 0.2 0.2];
  dii=10;
  for ii=xlim1:dii:xlim2
    for jj=ylim1:dii:ylim2
      clear u v
      u = U(jj,ii);
      v = V(jj,ii);
      s = S(jj,ii);
      if isnan(u) | s<0.02, continue; end;
  %    if res>0,
	u=u/s;
	v=v/s;
%	scl=25;
  %    end

      x0=ii;
      y0=jj;

      x1=x0+u*scl;
      y1=y0+v*scl;
      draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
    end
  end

  txtb = 'plot_arc04_Greenl_meanUV.m';
  bottom_text(txtb,'pwd',1,'fontsize',8);
  
  if s_fig>0
    fgnm = sprintf('%sarc08_110_meanUV_Lr%2.2i_%i',...
		   pthfig,plr,iyr);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r350',fgnm);
  end
  
  
end



