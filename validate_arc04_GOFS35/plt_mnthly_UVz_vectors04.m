% Plot monthly vectors
% no underlying colormap
% Length/color of the vectors are proportional to 
% the vector magnitude
% 
% Select regions
%
% from 0.04 HYCOM-CICE
% averaged for depth layer zz1-zz2
% calculated in mnthly_UVz_arc04.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2017;
YR2 = 2020;
imo1 = 1;
imo2 = 12;
basin = 'AO'; % AO - artctic ocean, SPNA, Baff, GIN, Bering


fmt_nm = 'mnthUVz'; % fmat name 


% Averaged over the depth layers:
%zz1=0;
%zz2=-50;
zz1=-200;
zz2=-400;


% Choose experiment:
ixx    = 9; % experiment name and dir - check with EXPT - expt 023
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;



regn = 'ARCc0.04';
rg = 9806;
hgg=1e20;

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/meanUV/',expt);

btx='plt_mnthly_UVz_vectors04.m';


%ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

xlim1 = 100;
xlim2 = nn-1;
ylim1 = 200;
ylim2 = mm-400;
dii = 20;
fct = 5;

fprintf(' Mean U/V vectors %s, %i-%i %i-%i\n',basin,YR1,YR2,imo1,imo2);

switch(basin)
 case('AO')
  XY = [ 597        3900
         597        1892
        2930        1892
        2930        3900];
  dii = 35; 
  fct = 15;

 case('SPNA')
  XY = [ 550        1450
         550          80
        2350          80
        2350        1450];
  dii=25;
  fct=15;

 case('Baff')
  XY = [600         2400
        600         915
        1300        915
        1300        2400];
  dii=20;
  fct=15;

 case('GIN')
  XY = [1500        2032
        1500         920
        2500         920
        2500        2032];
  dii=20;
  fct=15;

 case('Bering')
  XY = [760        4571
        760        3160
        1900       3160
        1900       4571];
 dii=20;
 fct=15;

end
xlim1 = min(XY(:,1));
xlim2 = max(XY(:,1));
ylim1 = min(XY(:,2));
ylim2 = max(XY(:,2));

xl1=xlim1;
xl2=xlim2;
yl1=ylim1;
yl2=ylim2;



usm=zeros(mm,nn);
vsm=zeros(mm,nn);
cc = 0;
for yr=YR1:YR2
	fmat = sprintf('%s%s_%4.4i-%4.4i_%i.mat',pthmat,fmt_nm,abs(zz1),abs(zz2),yr);
	fprintf('Loading %s\n',fmat);
	load(fmat);

	for imo=imo1:imo2
    fprintf('Calculating mean: %i/%i\n',yr,imo);
		for ik=1:12
			cc = cc+1;
			U = meanUV(ik).U;
			V = meanUV(ik).V;
			usm = usm+U;
			vsm = vsm+V;
		end
		U = usm./cc;
		V = vsm./cc;
	%  U = meanUV(imo).U;
	%  V = meanUV(imo).V;
	end
end

%S = sqrt(U.^2+V.^2);
nf = 1;
stl = sprintf('0.04 %3.3i-HYCOM2.3-CICE5, Mean U %i-%im, %i-%i, mo=%i-%i',...
              expt,abs(zz1),abs(zz2),YR1,YR2,imo1,imo2);

figure(nf); clf;
%ps=[1400         303        1002        1029];
ps=[1010  27 1259 1310];
set(gcf,'Position',ps);

xpos=[0.1,0.2,0.8,0.7];
axes('Position',xpos);
hmsk=HH;
hmsk(HH<0)=nan;
pcolor(hmsk); shading flat;
hold on;
%colormap([0.1 0.1 0.1]);
colormap([0.8 0.8 0.8]);
%freezeColors;

contour(HH,[-8000:500:-10],...
  'Color',[0.85 0.85 0.85],...
  'Linewidth',1.2);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],...
  'ylim',[ylim1 ylim2]);
set(gca,'xtick',[],'ytick',[]);
%clr=[0.9 0.9 0.9];
%plot_gridlines(45,10,1,clr,LON,LAT);
title(stl,'Fontsize',12,'Interpreter','none');


ubin=[0;5;10;15;20;400];
SCL=[1;2;2.5;3;3.5;4];
VCLR=[0.3,0.3,0.3;
      0., 0.2, 0.9;
      0, 0.7, 1;
      0.8, 0.5, 0;
      0.7, 0, 0];
CF=[0.9;0.4;0.3;0.3;0.3];
BT=[20;20;20;20;20];

if zz1<-100
	ubin=[0;2.5;5;7.5;10;400];
	SCL=[1;2;2.5;3;3.5;4];
	VCLR=[0.3,0.3,0.3;
				0., 0.2, 0.9;
				0, 0.7, 1;
				0.8, 0.5, 0;
				0.7, 0, 0];
	CF=[0.9;0.4;0.3;0.3;0.3];
	BT=[20;20;20;20;20];
end

%scl=10;
cf=0.3;
%beta=20;
lwd=1.6;
v_col=[0 0 0];
%dii=10;
fprintf('Plotting U vectors every %i grid points\n\n',dii);
for ii=xlim1:dii:xlim2
	for jj=ylim1:dii:ylim2
		clear u v
		u = U(jj,ii)*100;  % cm/s
		v = V(jj,ii)*100;  % cm/s
		if isnan(u) | isnan(v); continue; end;
    s = sqrt(u*u+v*v);

    beta=20;
    u=u/s;
    v=v/s;

    ibn=max(find(ubin<=s));
    
		a0    = SCL(ibn);
		v_col = VCLR(ibn,:);
		scl   = fct*a0;
		cf    = CF(ibn);
		beta  = BT(ibn); 

		x0=ii;
		y0=jj;

		x1=x0+u*scl;
		y1=y0+v*scl;
		draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);

	end
end


% Plot scale vector
xp1=xpos(1);
dxp=xpos(3);
yp1=xpos(2);
dyp=xpos(4);
axes('Position',[xp1,0.08,dxp,dyp]);
hold on;
y0=ylim1;
x0=xlim1+100;
dyy=40;
for ivv=1:5
  y0 = y0+dyy;
  a0 = SCL(ivv);
  cf = CF(ivv);
  v_col=VCLR(ivv,:);
  beta = BT(ivv);

  x1=x0+fct*a0;
  y1=y0;
  draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
  u1=ubin(ivv);
  u2=ubin(ivv+1);
  slgd = sprintf('%4.1f-%4.1f cm/s',u1,u2);
  text(x1+4*fct,y1,slgd,'Fontsize',14,'Color',v_col);

end

axis('equal');
set(gca,'xlim',[xlim1 xlim2],...
        'ylim',[ylim1 ylim2],...
        'visible','off'); 



bottom_text(btx,'pwd',1,'fontsize',8);





