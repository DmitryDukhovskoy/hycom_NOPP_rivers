% Approximate tracer/FW pathways
% following the maximum speed
% and using spline interpolation
% to connect these points
%
% The connection points = maximum speed
% around the subpolar gyre
% similar to how the GSA propagated
% using annual mean UV
%
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
s_mat  = 0;
f_mnth = 13; % plot <12 - month, >12 - annual mean

plr =1;  % U from plr layer
rg = 9806;
npth = 1; % # of pathway that is tracked


regn = 'ARCc0.08';
expt = 110;
pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

%YRPLT = [2011,2012,2013,2014,2015];
YRPLT = [2005];
np = length(YRPLT);

ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);
% Greenland
xlim1 = 300;
xlim2 = 1300;
ylim1 = 100;
ylim2 = 1100;
% SE Greenland:
%xlim1 = 600;
%xlim2 = 860;
%ylim1 = 350;
%ylim2 = 640;

% Specify points of max U
% based on the plot of the speed U
PTS = sub_define_GSApts(npth);

ik=1;
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

% Do parameteric spline interpolation
X=PTS.IJ(:,1);
Y=PTS.IJ(:,2);
[Xi,Yi]=sub_parametric_spline(X,Y);
Xi=Xi(:);
Yi=Yi(:);
JJs=round(Yi);
IIs=round(Xi);

% Eliminate repeated points
nS=length(Xi);
clear Dpp
Dpp(1)=0.1;
for ii=1:nS-1
  i1=IIs(ii);
  j1=JJs(ii);
  i2=IIs(ii+1);
  j2=JJs(ii+1);
  Dpp(ii)=sqrt((i2-i1).^2+(j2-j1).^2);
end

I=find(Dpp>1e-20);
IIs=IIs(I);
JJs=JJs(I);
IJs=[IIs,JJs];

nS=length(IIs);
clear Xl Yl
for ii=1:nS
  i0=IJs(ii,1);
  j0=IJs(ii,2);
  Xl(ii,1)=LON(j0,i0);
  Yl(ii,1)=LAT(j0,i0);
end;
INDs=sub2ind(size(HH),JJs,IIs);

nn=length(Xl);
x1=Xl(1);
y1=Yl(1);
clear Dst
clear dXX
for ii=1:nn
  x2=Xl(ii);
  y2=Yl(ii);
  dXX(ii)=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
  x1=x2;
  y1=y2;
  Dst(ii)=sum(dXX);
end

PTS.IJpathway=IJs;
PTS.Ipathway=INDs;
PTS.Dist=Dst;

if s_mat==1
  fmout=sprintf('%strcr_pathway_%2.2i.mat',pthmat,npth);
  fprintf('Saving %s\n',fmout);
  save(fmout,'PTS');
end





% ===========================
% Plotting
% ===========================


S = sqrt(U.^2+V.^2);
nf = 1;
stl = sprintf('ARCc0.08-%i, Mean U %i',expt,iyr);
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
pfld='speed';
hps = [0.92 0.1 0.03 0.8];
sub_plot_scalar(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
		'cmp',4,'c1',0,'c2',0.5,'clbpos',hps);
contour(HH,[-800 -800],'Color',[0. 0. 0.],'Linewidth',1.6);

i0=10;
j0=10;
%  Psi = stream_fn(U,V,DX,DY,i0,j0,'simps');
%keyboard

scl=13;
%scl=26;
cf=0.3;
beta=20;
lwd=1.;
v_col=[0.3 0.3 0.3];
dii=12;
% SE Greenland
%  scl=7;
%  cf=0.3;
%  beta=15;
%  lwd=1.;
%  v_col=[0.2 0.2 0.2];
%  dii=8;
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

plot(IIs,JJs,'c.','Markersize',12);
% Plot original nodes:
plot(X,Y,'k^','Markersize',5,...
          'MarkerFaceColor','k');
for ik=1:11
  xx=(ik-1)*1000;
  D=abs(Dst-xx);
  i0=find(D==min(D),1);
  plot(IIs(i0),JJs(i0),'b.','Markersize',26);
  stl=sprintf('%i',xx);
%  text(IIs(i0),JJs(i0),stl,'Fontsize',12);
end


txtb = 'tracer_pathways_meanUV.m';
bottom_text(txtb,'pwd',1,'fontsize',8);

if s_fig>0
  fgnm = sprintf('%sarc08_110_pthws_Lr%2.2i_%i',...
		 pthfig,plr,iyr);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r350',fgnm);
end

  




