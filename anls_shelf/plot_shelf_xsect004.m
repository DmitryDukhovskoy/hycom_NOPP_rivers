% HYCOM0.04
%
% analyze fluxes:
% Plot seasonal(or monthly) sections of U,T,S
% and  Vol, heat, FW fluxes
%
% averaged over specified years
% across sections on the Arctic shelf
% extracted in shelf_xsect004.m by months
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear


YR1=2006;
YR2=2006;
mo1=3;
mo2=3;
%nsct=1; % section to plot

Sref=34.9; % Similar to Le Bras and Straneo et al 2018
Sc0=33.8;

% =============
% Plot/anls fields:
% =============
f_plts=1;  % S
f_pltt=1;  % T


Tref=0;
rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.04';
%expt = 011;  % no Gr FWF
expt = 012;  % Gr FWF
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%3.3i/data_shelf/',expt);
%pthmat2=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%3.3i/data_GrSect/',expt);
btx = 'plot_shelf_xsect004.m';

fprintf('Plotting %s-%3.3i Monthly Greenl Shelf Sections %i/%i-%i/%i\n\n',...
	regn,expt,YR1,mo1,YR2,mo2);

ftopo = sprintf('%sdepth_%s_17DD.nc',pthtopo,regn); % 
fprintf('Getting bathymetry %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
%[DX,DY]=sub_dx_dy(LON,LAT);

%GC = sub_greenl_isobath(HH,LON,LAT);
%Ig=GC.cntr_Iindx;
%Jg=GC.cntr_Jindx;



TM=[];
cc=0;

for nsct=1:4

cc=0;
Tm=[];
Sm=[];
for YR=YR1:YR2
  yr=YR;
  
  for mo=mo1:mo2
    cc=cc+1;
    fmatu=sprintf('%s%3.3i_ArctShelf_xsct_dayTSZ_%4.4i%2.2i.mat',pthmat,expt,YR,mo);
    fprintf('Loading %s\n',fmatu);
    load(fmatu);

    if cc==1
      Tm=squeeze(nanmean(TSZ(nsct).Temp,1));
      Sm=squeeze(nanmean(TSZ(nsct).Saln,1));
    else
      Tm=Tm+squeeze(nanmean(TSZ(nsct).Temp,1));
      Sm=Sm+squeeze(nanmean(TSZ(nsct).Saln,1));
    end    

  end
  
end

Tw=Tm./cc;
Sw=Sm./cc;



f_map=0;
if f_map==1
  fn=20;
  figure(fn); clf;
  hold on;
  contour(HH,[0 0],'k');
  contour(HH,[-5000:1000:-10],'Color',[0.7 0.7 0.7]);

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%	 'Linewidth',2.5,'Color',[1. 0.6 0]);
    IIs=TSZ(ip).I;
    JJs=TSZ(ip).J;
    plot(IIs,JJs,'-',...
	 'Linewidth',2.5,'Color',[1. 0.4 0]);
  end
  axis('equal');
  set(gca,'xlim',[750 3200],...
	  'ylim',[2000 3800]);
  bottom_text(btx,'pwd',1);
end




LL=TSZ(nsct).Dist_origin*1e-3; % m
dx=1.0;
dL=diff(LL);
dL=[dL;dL(end)];
Ip=find(dL>dx);
%LLu=LL(Ip);
%Tw=Tw(:,Ip);
%Sw=Sw(:,Ip);

xl1=0;
xl2=round(max(LL));
yl1=-400;
yl2=0;
c1=26;
c2=35;
tc1=-2;
tc2=2;
Scntr=[20:1:33,33.5:0.5:35];
Tcntr=[-2:0.5:10];
dxx=20;
switch (nsct)
 case(1)
  yl1=-300;
  xl1=10;
 case(2)
  dxx=50;
 case(3)
  yl1=-350;
  dxx=40;
 case(4);
  yl1=-500;
  dxx=50;
end

% ============================
% Salinity
% Summer and winter months
% ============================
Sc0=30.0;
if f_plts==1
  nf=nsct;
  stlS=sprintf('004-%3.3i, S, %i-%i, mo:%i-%i, Sfront=%4.2f',expt,YR1,YR2,mo1,mo2,Sc0);
  sub_plot_S_Shelf004(nf,TSZ,Sw,btx,stlS,Sc0,nsct,xl1,xl2,yl1,yl2,Scntr,c1,c2,dxx);
%  contour(LLu,ZM,Us,[-1:0.1:-0.001],'--','Color',[0.8 0.8 0.8]);
%  contour(LLu,ZM,Us,[-0.5 -0.5],'--','Color',[0.8 0.8 0.8],'Linewidth',1.6);
  %set(gcf,'Position',[912 529 1544 804]);
  bottom_text(btx,'pwd',1,'Position',[0.08 0.15 0.4 0.05]);

  
end

% ============================
% Temperature - not finished 
% need sub-function 
% ============================
Tc0=0;
if f_pltt==1
  nf=10+nsct;
  stlS=sprintf('004-%3.3i, T, %i-%i, mo:%i-%i, Tfront=%4.2f',expt,YR1,YR2,mo1,mo2,Tc0);
  sub_plot_T_Shelf004(nf,TSZ,Tw,btx,stlS,Tc0,nsct,xl1,xl2,yl1,yl2,Tcntr,tc1,tc2,dxx);
%  sub_plot_T_Shelf004(nf,UTSZ,Ss,btx,stlS);
  bottom_text(btx,'pwd',1,'Position',[0.08 0.15 0.4 0.05]);
end

end



