% Quick plot monthly UV fields
% from 0.04 HYCOM-CICE
% averaged for depth layer zz1-zz2
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2018;
YR2 = 2020;
imo1 = 1;
imo2 = 6;


fmt_nm = 'mnthUVz'; % fmat name 


% Averaged over the depth layers:
zz1=0;
zz2=-50;


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

btx='plt_mnthly_UVz_arc04.m';


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


usm=zeros(mm,nn);
vsm=zeros(mm,nn);
cc = 0;
for yr=YR1:YR2
	fmat = sprintf('%s%s_%4.4i-%4.4i_%i.mat',pthmat,fmt_nm,abs(zz1),abs(zz2),yr);
	fprintf('Loading %s\n',fmat);
	load(fmat);

	for imo=imo1:imo2
    fprintf('Calculating mean: %i\%i\n',yr,imo);
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

S = sqrt(U.^2+V.^2);
nf = 1;
stl = sprintf('0.04 %3.3i-HYCOM2.3-CICE5, Mean U %i-%im, %i-%i, mo=%i-%i',...
              expt,abs(zz1),abs(zz2),YR1,YR2,imo1,imo2);

figure(nf); clf;
%ps=[1400         303        1002        1029];
ps=[1010  27 1259 1310];
set(gcf,'Position',ps);

sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,'c1',0,'c2',0.1,'cmp',8);


scl=10;
cf=0.3;
beta=20;
lwd=1.2;
v_col=[0 0 0];
dii=40;
fprintf('Plotting U vectors every %i grid points',dii);
for ii=xlim1:dii:xlim2
	for jj=ylim1:dii:ylim2
		clear u v
		u = U(jj,ii);
		v = V(jj,ii);
		s = S(jj,ii);
		if isnan(u), continue; end;
%    if res>0,
    u=u/s;
    v=v/s;
    scl=50;
%    end

		x0=ii;
		y0=jj;

		x1=x0+u*scl;
		y1=y0+v*scl;
		draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
	end
end

bottom_text(btx,'pwd',1,'fontsize',8);





