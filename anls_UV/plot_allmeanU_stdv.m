% Plot overall-mean UV fields
% with standard deviation ellipses at selected locations
% derived in mnthly_arc08_UV.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;

fprintf('s_fig= %i, long-term mean UV\n',s_fig);

plr =1;  % U from plr layer
rg = 9806;


regn = 'ARCc0.08';
expt = 110;
%pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/110/fig_meanUV/';
pthfig  = '/Net/mars/ddmitry/hycom/ARCc0.08/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

%YRPLT = [1994:2016];
%YRPLT = [2008];
YRPLT = [1993:2016];
np = length(YRPLT);

ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);
xlim1 = 20;
xlim2 = nn-1;
ylim1 = 100;
ylim2 = mm-100;

% Select locations where to plot StDev Ell.
IJp =[894        1159
        1025        1461
        1019        1621
        1221        1589
        1013        1746
         908        1430
         865        1729
         458        1589
         629        1660
         740        1515
         586        1410
         822        1632
         757        1296
        1013        1241];
npp = length(IJp);

cc = 0;
clear up vp
for ik=1:np
  iyr = YRPLT(ik);
  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  if cc==0
    usm=zeros(mm,nn);
    vsm=zeros(mm,nn);
  end
  
  for im=1:12
    cc = cc+1;
    U = meanUV(im).U;
    V = meanUV(im).V;
    usm = usm+U;
    vsm = vsm+V;
    for ipp = 1:npp
      i0=IJp(ipp,1);
      j0=IJp(ipp,2);
      up(cc,ipp)=meanUV(im).U(j0,i0);
      vp(cc,ipp)=meanUV(im).V(j0,i0);
    end
  end
  
  fprintf('cc=%i\n',cc);
end

U = usm./cc;
V = vsm./cc;

S = sqrt(U.^2+V.^2);
nf = 1;
stl = sprintf('ARCc0.08-%i, Mean U %i-%i',expt,YRPLT(1), YRPLT(end));
sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);

i0=10;
j0=10;
%  Psi = stream_fn(U,V,DX,DY,i0,j0,'simps');
%keyboard

scl=10;
cf=0.3;
beta=20;
lwd=1.2;
v_col=[0.3 0.3 0.3];
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
      scl=25;
%    end

    x0=ii;
    y0=jj;

    x1=x0+u*scl;
    y1=y0+v*scl;
    draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
  end
end

% Calculate StDev ellipses
% using monthly data
fprintf('Plotting StDev ellipses ...\n');
for ipp = 1:npp
  i0=IJp(ipp,1);
  j0=IJp(ipp,2);

  clear uu vv
  uu = up(:,ipp);
  vv = vp(:,ipp);

  sub_std_ellipse(uu,vv,i0,j0);

end

set(gca,'xlim',[280 1600],...
	'ylim',[580 1950]);

txtb = 'plot_allmeanUV_stdv.m';
bottom_text(txtb,'pwd',1,'fontsize',8);

if s_fig>0
  fgnm = sprintf('%sarc08_110_meanUVstd_Lr%2.2i_%i-%i',...
		 pthfig,plr,YRPLT(1),YRPLT(end));
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r350',fgnm);
end






