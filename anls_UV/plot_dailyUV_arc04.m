% Plot daily UV fields
% fmor instant output archive files archv
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear


f_ice = 1; % contour sea ice conc ontop
plr =5;  % U from plr layer
rg = 9806;


regn = 'ARCc0.04';
expt = 022;
%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

day_plot = datenum(2017,6,1);
dnmb = day_plot;
dv = datevec(dnmb);


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
xlim1 = 400;
xlim2 = 3160;
ylim1 = 1200;
ylim2 = 3800;

xl1=xlim1;
xl2=xlim2;
yl1=ylim1;
yl2=ylim2;


yr=dv(1);
iday=dnmb-datenum(yr,1,1)+1;
%pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output/'; % test simulation with coupling
pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i/',yr);
%pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_0901gofs35/'; % sea ice restart on 09/01 from GOFS3.5 GLBc


fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);



if ~exist(fina,'file');
		fprintf('Not found: %s\n\n',fina);
end


if f_ice==1
  YR=yr;
  DV = datevec(dnmb);
  mo = DV(2);
  mday = DV(3);
  pthbini = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%i_cice/',expt,YR);
  flnm = sprintf('%s022_cice.%i-%2.2i-%2.2i.nc',pthbini,YR,mo,mday);
  Ci = squeeze(nc_varget(flnm,'aice'));
  Ci(HH>=0)=nan;
  Ci(Ci<1e-10)=0;
end


DV=datevec(dnmb);


fprintf('U/V V.Lr=%i, %4.4i_%2.2i_%2.2i: %s\n',plr,DV(1:3),fina);
% in archm u_vel=utot; total velocity
% in archv - need to add barotropic u
%  keyboard
%    ilr=1;
[A,n,m,l] = read_hycom(fina,finb,'u_btrop');
[B,n,m,l] = read_hycom(fina,finb,'v_btrop');
Ub=squeeze(A);
Ub(Ub>1e6)=nan;
Vb=squeeze(B);
Vb(Vb>1e6)=nan;

[Fu,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',plr);
[Fv,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',plr);
U=squeeze(Fu);
U(U>1e6)=nan;
Ut=U+Ub;
V=squeeze(Fv);
V(V>1e6)=nan;
Vt=V+Vb;

SS=sqrt(Ut.^2+Vt.^2);

%keyboard

cc1=0;
cc2=0.5;
CMP = create_colormap_WBYR(200,cc1,cc2);
cmp = CMP.colormap;


isl = max(strfind(fina,'/'));
nf = 1;
%stl = sprintf('CPL TEST ARCc0.04-%3.3i, U Lr %i, %s',expt,plr,fina(isl+1:end));
stl = sprintf('ARCc0.04-%3.3i, U Lr %i, %s',expt,plr,fina(isl+1:end));
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
pfld='speed';
hps = [0.93 0.1 0.035 0.8];
%keyboard


Lmsk = HH*0;
Lmsk(HH<0)=1;
lcmp=[0.5 0.5 0.5; 0.9 0.9 0.9];


figure(nf); clf;
set(gcf,'Position',[1262 431 931 900]);
axes('Position',[0.05 0.1 0.8 0.8]);
hold on;
pcolor(Lmsk); shading flat;
colormap(lcmp);
freezeColors;

pcolor(SS); shading flat;
caxis([cc1 cc2]);
colormap(cmp);

if f_ice==1
  contour(Ci,[0.01 0.01],'Color',[0.8 0.3 0],'Linewidth',2);
end

axis('equal');

set(gca,'xlim',[xl1 xl2],...
        'ylim',[yl1 yl2],...
        'xtick',[],...
        'ytick',[]);
%        'Color',[0 0 0]);

title(stl,'Fontsize',12,'Interpreter','none');

hb=colorbar;
set(hb,'Position',[0.87 0.11 0.02 0.8],...
       'Ticklength',0.03,...
       'Fontsize',14);



f_vectors=0;
if f_vectors==1
		i0=10;
		j0=10;

		scl   = 150;
		cf=0.3;
		beta=15;
		lwd=1.;
		v_col=[0.1 0.1 0.1];
		dii=40;
		for ii=xlim1:dii:xlim2
				for jj=ylim1:dii:ylim2
						clear u v
						u = U(jj,ii);
						v = V(jj,ii);
						s = SS(jj,ii);
						if isnan(u) | s<0.02, continue; end;
						x0=ii;
						y0=jj;

						x1=x0+u*scl;
						y1=y0+v*scl;
						draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
				end
		end

end

btx = 'plot_dailyUV_arc04';
bottom_text(btx,'pwd',1);
  


