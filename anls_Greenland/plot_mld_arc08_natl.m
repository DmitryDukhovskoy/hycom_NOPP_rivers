% Plot MLD
% see mld_arc08_natl.m where monthly means are calculated:
% Calculate MLD from the model
% Using Kara and my methods
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

yr1=2010;
yr2=2016;
mthd='Kara';
%mthd='DD';
%mthd='dRho';

regn = 'ARCc0.08';
expt = 110; % no Greenland runoff  
%expt = 112;  % Greenland runoff

s_fig = 0;
rg = 9806;

iyr=2012;
im=2;

fprintf('%3.3i MLD-%s, %4.4i/%2.2i\n',expt,mthd,iyr,im);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

%[DX,DY]=sub_dx_dy(LON,LAT);
%Acell=DX.*DY; % Grid cell area, m2

cc=0;
for iyr=yr1:yr2
  cc=cc+1;
  fmat = sprintf('%sarc08_%3.3i_mnthly_MLD_%i%2.2i.mat',...
		     pthmat,expt,iyr,im);
  fprintf('Loading %s\n',fmat);
  load(fmat); % struct array MILD
  switch(mthd);
   case('DD')
    mld=MILD.MLD_DD; % DD
   case('Kara')
    mld=MILD.MLD_Kara; % Kara
   case('dRho');
    mld=MILD.MLD_dRho; % dlt Rho
  end
  if cc==1, msum=mld*0; end;
  I=find(isnan(mld));
  mld(I)=-50;
  msum=msum+mld;
end
if cc>1
  mld=msum/cc;
end

Iocn=MILD.Iocn;
MLD=HH*nan;
MLD(Iocn)=mld;


hmsk=HH;
hmsk(HH<0)=nan;

c1=-500;
c2=-50;
dcc=round(abs(c2-c1)/10);
CMP = colormap_WBP(200,c1,c2);
cmp = flipud(CMP.colormap);

figure(1); clf;
pcolor(MLD); shading flat;
hold on;
contour(HH,[-800 -800],'Color',[0.7 0.7 0.7],'Linewidth',1);
%caxis([0 2]);
colormap(cmp);
caxis([c1 c2]);

xlim1=390;
xlim2=1230;
ylim1=10;
ylim2=1150;
axis('equal');
set(gca,'xlim',...
	[xlim1 xlim2],...
	'ylim',[ylim1 ylim2],...
	'xtick',[],'ytick',[],...
	'Color',[0.5 0.5 0.5]);
%clr=[0.9 0.9 0.9];
%plot_gridlines(45,10,1,clr,LON,LAT);
%stl=sprintf('arc08 %3.3i MLD_DD %i/%2.2i',expt,iyr,im);
stl=sprintf('ARCc0.08-%3.3i MLD-%s %i-%i/%2.2i',expt,mthd,yr1,yr2,im);
title(stl,'Fontsize',12,'Interpreter','none');

hbb=colorbar;
set(hbb,'Position',[0.86 0.12 0.02 0.8],...
	'TickLength',0.02,...
	'Ticks',[c1:dcc:c2],...
	'Fontsize',14);

%  hb = colorbar;
%set(hb,'position',hps,...
%	 'TickLength',0.04,...
%	 'Fontsize',14);


freezeColors;
pcolor(hmsk); shading flat;
colormap([0 0 0]);
freezeColors;

colormap(cmp);

set(gcf,'Position',[610 18 1186 1321]);
btx='plot_mld_arc08_natl.m';
bottom_text(btx,'pwd',1);
