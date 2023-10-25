% Plot monthly climatology mean UV fields
% Can do 1 month at a time:
% im1=im2 for all years (monthly climatology)
% or average over all months
% im1:im2 if im1/=im2
%
% with standard deviation ellipses at selected locations
% derived in mnthly_arc08_UV.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat  = 6; % =0 - load saved,  =1 - save
im1 = 12; % im1=im2 plot month=1,...,12 - 1 month at a time, 
im2 = 12; % otherwise - average over N years and all months


plr =1;  % U from plr layer
rg = 9806;

regn = 'ARCc0.08';
expt = 110;
%pthfig  = '/Net/ocean/ddmitry/hycom/ARCc0.08/110/fig_meanUV/';
pthfig  = '/Net/mars/ddmitry/hycom/ARCc0.08/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fout = sprintf('%s0.08_%3.3i_meanUV50m_Fram_%2.2i_%2.2i-1993-2016.mat',pthmat,expt,im1,im2);

YRPLT = [1993:2016];
%YRPLT = [2013];
np = length(YRPLT);

fprintf('Mean U: %i-%i/ %i-%i\n',im1,im2,YRPLT(1),YRPLT(end));

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);

% Select locations where to plot StDev Ell.
% Arctic Ocean and N. Atlantic
% may need to adjust scales
% as speed is very different in 
% N. Atl and Arctic O. 
IJp =[         958         977
         977         972
        1001         970
        1019         968
        1038         968
        1060         969
         944         959
         970         934
         990         922
        1006         921
        1031         922
        1055         927
        1072         932
         934         930
         954         920
         927         889
         955         881
         978         869
         998         872
        1029         880
        1069         886
        1091         897
         924         841
         949         830
         977         821
         997         819
        1046         820
        1076         827
        1108         839
        1132         856
        1155         873];
    
npp = length(IJp);

fprintf('Selected points for  mean U and std ellipse\n');
fprintf('  Lon     Lat\n');
for ipp=1:npp
  i1=IJp(ipp,1);
  j1=IJp(ipp,2);
  
  xx=LON(j1,i1);
  yy=LAT(j1,i1);
  
  fprintf('%7.3f   %7.3f\n',xx, yy);
end


if s_mat==1
  usm=zeros(mm,nn);
  vsm=zeros(mm,nn);
  ikk=0;
  for ik=1:np
    iyr = YRPLT(ik);
    fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
    fprintf('Loading %s\n',fmat);
    load(fmat);

    for imo=im1:im2
      ikk=ikk+1;
      U = meanUV(imo).U;
      V = meanUV(imo).V;
      usm = usm+U;
      vsm = vsm+V;

    % Save data for Calculating StDev
      for ipp = 1:npp
	i0=IJp(ipp,1);
	j0=IJp(ipp,2);

	Up(ikk,ipp)=meanUV(imo).U(j0,i0);
	Vp(ikk,ipp)=meanUV(imo).V(j0,i0);

      end

    end

  end
  U=usm/ikk;
  V=vsm/ikk;

  UMN.U=U;
  UMN.V=V;
  UMN.yr1=YRPLT(1);
  UMN.yr2=YRPLT(end);
  UMN.tser_IJ=IJp;
  UMN.Utser=Up;
  UMN.Vtser=Vp;
  fprintf('Saving %s\n',fout);
  save(fout,'UMN');

else
  fprintf('Loading %s\n',fout);
  load(fout);
  
end

U=UMN.U;
V=UMN.V;
Up=UMN.Utser;
Vp=UMN.Vtser;

xlim1 = 890;
xlim2 = 1150;
ylim1 = 780;
ylim2 = 990;

S = sqrt(U.^2+V.^2);
nf = 1;
stl = sprintf('ARCc0.08-%i, U, mo=%i/%i, %i-%i',...
	      expt,im1,im2,min(YRPLT),max(YRPLT));
c1=0;
c2=0.2;
f_cmp=3;

%sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,...
%		   ylim1,ylim2,LON,LAT,...
%		   stl,'c1',c1,'c2',c2,'cmp',f_cmp);

figure(1); clf;
set(gcf,'Position',[1321 219 1010 1098])
axes('Position',[0.05 0.2 0.9 0.7]);
hmsk=HH;
hmsk(HH<0)=nan;
pcolor(hmsk); shading flat;
hold on;
%colormap([0.1 0.1 0.1]);
colormap([0.8 0.8 0.8]);
%freezeColors;

contour(HH,[-8000:500:-10],...
	'Color',[0.6 0.6 0.6],...
	'Linewidth',1.2);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],...
	'ylim',[ylim1 ylim2]);
set(gca,'xtick',[],'ytick',[]);
%clr=[0.9 0.9 0.9];
%plot_gridlines(45,10,1,clr,LON,LAT);
title(stl,'Fontsize',12,'Interpreter','none');

% Plot moorings
SBE = Fram_moorAWI;
nf=length(SBE);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  XY=[x,y];
  IJ = sub_XY2indx(XY,LON,LAT);
  plot(IJ(1),IJ(2),'r.','Markersize',20);
  nm=SBE(ik).Name;
%  text(IJ(1),IJ(2),nm,'Fontsize',12);
end




scl=1.1;
cf=0.35;
beta=13;
lwd=1.6;
v_col=[0.3 0.3 0.3];
dii=10;
for ii=xlim1:dii:xlim2
  for jj=ylim1:dii:ylim2
    clear u v
    u = U(jj,ii)*100; % m/s -> cm/s
    v = V(jj,ii)*100;
    s = sqrt(u*u+v*v);
    if isnan(u); continue; end;
    
    u=u/s;  % scale to speed, i.e. unit vectors
    v=v/s;
    if s>=1 & s<5
      scl=4;
      v_col=[0.5 0.5 0.5];
      cf=0.8;
%      beta=30;
    elseif s>=5 & s<10
      scl=10;
      v_col=[0. 0.6 1];
      cf=0.4;
    elseif s>=10 & s<15
      scl=21;
      v_col=[0.9 0.5 0.];
      cf=0.3;
    elseif s>=15
      scl=25;
      v_col=[0.5 0. 0.];
      cf=0.3;
    end

    x0=ii;
    y0=jj;

    x1=x0+u*scl;
    y1=y0+v*scl;
    draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
    
%    keyboard
    
  end
end

% Calculate StDev ellipses
% using monthly data
fprintf('Plotting StDev ellipses ...\n');
%um2cm=100;
%sclv=20;
sclv=1.;
lwd=1.8;
v_col=[0 0.2 0.5];
for ipp = 1:npp
  clear u v
  i0=IJp(ipp,1);
  j0=IJp(ipp,2);
  if i0>=xlim1 & i0<=xlim2 & ...
     j0>=ylim1 & j0<=ylim2
    u=Up(:,ipp)*100;
    v=Vp(:,ipp)*100;
    sub_std_ellipse_v2(u,v,i0,j0,sclv,v_col,cf,lwd,beta);
  end
end


%set(gca,'xlim',[280 1600],...
%	'ylim',[10 1950]);

% plot scale vectors:
axes('Position',[0.05 0.05 0.9 0.1]);
ylim22=ylim1+(ylim2-ylim1)/7;
x0=xlim1+10;
y0=ylim22-0.5;
dxx=5;
dyy=7;

u=1; v=0;
% 1-5 cm
scl=4;
v_col=[0.5 0.5 0.5];
cf=0.8;
%beta=30;

x1=x0+u*scl;
y1=y0+v*scl;
draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
text(x1+dxx,y1,'1-5 cm/s','Fontsize',14,'Color',v_col);

%beta=20;
% 5-10 cm
scl=10;
v_col=[0. 0.6 1];
cf=0.4;
y0=y0-dyy;
x1=x0+u*scl;
y1=y0+v*scl;
draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
text(x1+dxx,y1,'5-10 cm/s','Fontsize',14,'Color',v_col);

%10-15
scl=21;
v_col=[0.9 0.5 0.];
%v_col=[0. 0. 0.6];
cf=0.3;
y0=y0-dyy;
x1=x0+u*scl;
y1=y0+v*scl;
draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
text(x1+dxx,y1,'10-15 cm/s','Fontsize',14,'Color',v_col);

%>15
scl=25;
v_col=[0.5 0. 0.];
cf=0.3;
y0=y0-dyy;
x1=x0+u*scl;
y1=y0+v*scl;
draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
text(x1+10,y1,'>15 cm/s','Fontsize',14,'Color',v_col);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],...
        'ylim',[ylim1 ylim22]);


f_legend=1;
if f_legend>0
  lwd=1.8;
  v_col=[0 0.2 0.5];

  y0=ylim22-0.5;
  x0=1000;
  text(x0,y0,'<U> and StDev Ellipses:','Fontsize',14);
  ulg=20;
  vlg=0;
  ii=x0;
  jj=y0-dyy;
  col=v_col;
  [Xh,Yh]=get_arrowF(ii,ii+ulg*sclv,jj,jj+vlg*sclv,cf,beta,col,lwd);
  ptx=text(Xh(1)+5,Yh(1),sprintf('<U> %2.2i cm/s',ulg));
  set(ptx,'Color',col,'Fontsize',14);

  ulg=10;
  vlg=0;
  ii=x0;
  jj=jj-dyy;
  col=v_col;
  [Xh,Yh]=get_arrowF(ii,ii+ulg*sclv,jj,jj+vlg*sclv,cf,beta,col,lwd);
  ptx=text(Xh(1)+5,Yh(1),sprintf('<U> %2.2i cm/s',ulg));
  set(ptx,'Color',col,'Fontsize',14);
end
set(gca,'visible','off');



txtb = 'plot_meanUVmonth_stdv_Fram.m';
bottom_text(txtb,'pwd',1,'fontsize',7);





