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
s_mat = 0; % save overall mean; =0 - load saved and plot
f_mnth = 13; % plot <12 - month, >12 - annual mean

plr =1;  % U from plr layer
rg = 9806;


regn = 'ARCc0.08';
expt = 110;
pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

YRPLT = [1993:2016];
%YRPLT = [1998];
np = length(YRPLT);

fout = sprintf('%sarc08_%i_allmeanUV_lr%2.2i_%i-%i.mat',...
	       pthmat,expt,plr,YRPLT(1),YRPLT(end));

ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);

if s_mat>0
  cc = 0;
  usm=zeros(mm,nn);
  vsm=zeros(mm,nn);
  for ik=1:np
    iyr = YRPLT(ik);
    fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
    fprintf('Loading %s\n',fmat);
    load(fmat);

    if f_mnth>12
      for ik=1:12
	cc = cc+1;
	U = meanUV(ik).U;
	V = meanUV(ik).V;
	usm = usm+U;
	vsm = vsm+V;
      end
    end
  end
  U = usm./cc;
  V = vsm./cc;

  fprintf('Saving mean %s\n',fout);
  save(fout,'U','V','HH','LON','LAT');

else
  fprintf('Loading mean %s\n',fout);
  load(fout);
end

xlim1 = 100;
xlim2 = 1400;
ylim1 = 10;
ylim2 = 1200;

S = sqrt(U.^2+V.^2);
nf = 1;
stl = sprintf('ARCc0.08-%i, Mean U %i-%i',expt,YRPLT(1),YRPLT(end));
c1=0;
c2=0.2;
f_cmp=3;
sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,...
		   ylim1,ylim2,LON,LAT,...
		   stl,'c1',c1,'c2',c2,'cmp',f_cmp);

scl=30;
cf=0.35;
beta=25;
lwd=1.2;
v_col=[0 0 0];
dii=40;
for ii=10:dii:nn
  for jj=10:dii:mm
    clear u v
    u = U(jj,ii);
    v = V(jj,ii);
    s = S(jj,ii);
    if isnan(u), continue; end;
%    if res>0,
      u=u/s;
      v=v/s;
%    end

    x0=ii;
    y0=jj;

    x1=x0+u*scl;
    y1=y0+v*scl;
    draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
  end
end

txtb = 'overall_meanUV.m';
bottom_text(txtb,'pwd',1,'fontsize',8);

if s_fig>0
  fgnm = sprintf('%sarc08_110_allmeanUV_Lr%2.2i_%i',...
		 pthfig,plr,iyr);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r350',fgnm);
end




