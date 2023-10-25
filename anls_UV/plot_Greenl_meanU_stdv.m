% Plot overall-mean UV fields for Greenland
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
f_getu = 0; 

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
%xlim1 = 250;
%xlim2 = 1210;
%ylim1 = 150;
%ylim2 = 1150;
xlim1 = 330;
xlim2 = 1280;
ylim1 = 75;
ylim2 = 1090;

% Select locations where to plot StDev Ell.
IJp = sub_Green_stdPts;    
npp = length(IJp);

cc = 0;
clear up vp
fmatU = sprintf('%smeanU_StdDev_Greenl.mat',pthmat);
if f_getu==1
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

  save(fmatU,'up','vp','U','V');

else
  load(fmatU);
end

  
S = sqrt(U.^2+V.^2);

%
% ==================
% Plot map
% ==================
hmsk=HH*0+1;
hmsk(HH<0)=nan;

% Colormap:
c1 = 0;
c2 = 0.5;
nint = 200;
%CMP = colormap_sclr1(nint,c1,c2);
%CMP = create_colormap_WBYR(nint,c1,c2);
CMP = create_colormap_paleWBYR(nint,c1,c2);
cmp = CMP.colormap;
for ik=1:10;
  cmp(ik,:)=[1 1 1];
end
cmp = smooth_colormap(cmp,7);
cnt = CMP.intervals;
nint=length(cmp);

figure(1); clf;
axes('Position',[0.05 0.08 0.85 0.85]);

pcolor(hmsk); shading flat;
colormap([0.3 0.3 0.3]);
freezeColors;
hold on


pcolor(S); shading flat;
%hold on;
%contour(HH,[0 0],'k','Linewidth',1);
%contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',1);
%caxis([0 2]);
caxis([c1 c2]);
colormap(cmp);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],...
	'ylim',[ylim1 ylim2]);
set(gca,'xtick',[],'ytick',[]);

% =================
%
% Plot streamlines
%
% =================
yrF = 2005;
arrow_dst = 5; % min dst btw arrows on strmlines
smx_L = 350; % max length of strml
smn_L = 40; % min length of strml
dlim = 8; % min dist between strmlines
dip  = 2;   % step btw prtcles to plot strmlines
strm_clr = [0.6 0.6 0.6];
v_col = [0.5 0.5  0.5];
lhead = 8; % size of arrow head 
sub_plot_strml_meanU(expt,yrF,arrow_dst,smx_L,smn_L,...
			      dlim,dip,strm_clr,v_col,lhead);



% ==================================
%
% Plot mean vectors + ellipses
%
% ==================================
% Calculate StDev ellipses
% using monthly data
fprintf('Plotting StDev ellipses ...\n');
uscl=10;
for ipp = 1:npp
  i0=IJp(ipp,1);
  j0=IJp(ipp,2);

  clear uu vv
  uu = up(:,ipp);
  vv = vp(:,ipp);

  sub_std_ellipse(uu,vv,i0,j0,uscl);

end






%clr=[0.9 0.9 0.9];
%plot_gridlines(45,10,1,clr,LON,LAT);
stl=sprintf('ARCc0.08-%3.3i, mean U, StDev Ellipses',expt);
title(stl,'Fontsize',12,'Interpreter','none');

hght=[];
lngth=[];
mint=20;
mbx=mint;
fsz=12;
bxc='k';
posc=[0.9 0.1 0.8 0.05];
aend=0;
[az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);


txtb = 'plot_Greenl_meanUV_stdv.m';
bottom_text(txtb,'pwd',1);







