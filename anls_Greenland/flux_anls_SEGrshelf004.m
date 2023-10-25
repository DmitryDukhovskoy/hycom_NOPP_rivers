% HYCOM0.04
%
% analyze fluxes:
% Plot seasonal(or monthly) sections of U,T,S
% and  Vol, heat, FW fluxes
%
% averaged over specified years
% across sections on the Gr Shelf
% see also anls_EGC004.m
% extracted in flux_SEGrShelf_xsct004.m
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
YRS = [YR1:YR2];
  
Sref=34.9; % Similar to Le Bras and Straneo et al 2018
Sc0=33.8;

% =============
% Plot/anls fields:
% =============
f_plts=1;  % S
f_pltt=0;  % T


Tref=0;
rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.04';
%expt = 011;  % no Gr FWF
expt = 012;  % Gr FWF
pthfig  = '/Net/ocean/ddmitry/hycom/ARCc0.04/110/fig_GrSct/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%3.3i/data_GrSect/',expt);
%pthmat2=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%3.3i/data_GrSect/',expt);
btx = 'flux_anls_SEGrshelf004.m';

fprintf('Plotting %s-%3.3i Monthly Greenl Shelf Sections %i-%i\n\n',...
	regn,expt,YR1,YR2);

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


fmatu=sprintf('%s%3.3i_SEGreenlSh_xsct_dayUTSZ_%i.mat',pthmat,expt,YR1);
fprintf('Loading %s\n',fmatu);
load(fmatu);
ZZ=UTSZ(1).ZZlevels;
ZM=UTSZ(1).ZM;
nzm=length(ZM);
nzz=length(ZZ);
dZ=abs(diff(ZZ));

% Summer and winter months
js1=7;
js2=7;
jw1=1;
jw2=1;

% Winter mean
TM=[];
nsct=1;
cc=0;
for YR=YR1:YR2
  yr=YR;
  
  cc=cc+1;

  TM=UTSZ.Time;
  DV=datevec(TM);
  Iw=find(DV(:,2)>=jw1 | DV(:,2)<=jw2);
  U=UTSZ.Unrm(Iw,:,:);
  T=UTSZ.Temp(Iw,:,:);
  S=UTSZ.Saln(Iw,:,:);
  
  if cc==1
    Um=squeeze(nanmean(U,1));
    Tm=squeeze(nanmean(T,1));
    Sm=squeeze(nanmean(S,1));
  else
    Um=Um+squeeze(nanmean(U,1));
    Tm=Tm+squeeze(nanmean(T,1));
    Sm=Sm+squeeze(nanmean(S,1));
  end
  
end
Uw=Um./cc;
Tw=Tm./cc;
Sw=Sm./cc;

% Summer mean:
cc=0;
clear Um Tm Sm
for YR=YR1:YR2
  yr=YR;
  
  cc=cc+1;

  TM=UTSZ.Time;
  DV=datevec(TM);
  Iw=find(DV(:,2)>=js1 & DV(:,2)<=js2);
  U=UTSZ.Unrm(Iw,:,:);
  T=UTSZ.Temp(Iw,:,:);
  S=UTSZ.Saln(Iw,:,:);
  
  if cc==1
    Um=squeeze(nanmean(U,1));
    Tm=squeeze(nanmean(T,1));
    Sm=squeeze(nanmean(S,1));
  else
    Um=Um+squeeze(nanmean(U,1));
    Tm=Tm+squeeze(nanmean(T,1));
    Sm=Sm+squeeze(nanmean(S,1));
  end
  
end
Us=Um./cc;
Ts=Tm./cc;
Ss=Sm./cc;

% Project U on normal plain to get rid of
% "steps" in +/- flux across zigzagin transect
LL=UTSZ.Dist_origin*1e-3; % m
dx=2.0;
dL=diff(LL);
dL=[dL;dL(end)];
Ip=find(dL>dx);
LLu=LL(Ip);
Us=Us(:,Ip);
Uw=Uw(:,Ip);

LL=UTSZ.Dist_origin*1e-3; % km
% ============================
% Salinity
% Summer and winter months
% ============================
if f_plts==1
  nf=2;
  stlS=sprintf('004-%3.3i, S, %i-%i, mo:%i-%i, Sfront=%4.2f',expt,YR1,YR2,js1,js2,Sc0);
  sub_plot_S_GrShelf004(nf,UTSZ,Ss,btx,stlS,Sc0);
  contour(LLu,ZM,Us,[-1:0.1:-0.001],'--','Color',[0.8 0.8 0.8]);
  contour(LLu,ZM,Us,[-0.5 -0.5],'--','Color',[0.8 0.8 0.8],'Linewidth',1.6);

  %set(gcf,'Position',[912 529 1544 804]);
  bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.4 0.05]);

  nf=3;
  stlS=sprintf('004-%3.3i, S, %i-%i, mo:%i-%i, Sfront=%4.2f',expt,YR1,YR2,jw1,jw2,Sc0);
  sub_plot_S_GrShelf004(nf,UTSZ,Sw,btx,stlS,Sc0);
  contour(LLu,ZM,Uw,[-1:0.1:-0.001],'--','Color',[0.8 0.8 0.8]);
  contour(LLu,ZM,Uw,[-0.5 -0.5],'--','Color',[0.8 0.8 0.8],'Linewidth',1.6);

  %set(gcf,'Position',[912 529 1544 804]);
  bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.4 0.05]);

  
  
end

% ============================
% Temperature - not finished 
% need sub-function 
% ============================
if f_pltt==1
% Summer and winter months

  nf=7;
  stlS=sprintf('004-%3.3i, T, %i-%i, mo:%i-%i',expt,YR1,YR2,js1,js2);
  sub_plot_T_GrShelf004(nf,UTSZ,Ss,btx,stlS);
  
  stlS=sprintf('004-%3.3i, T, %i-%i, mo:%i-%i',expt,YR1,YR2,jw1,jw2);
  sub_plot_T_GrShelf004(nf,UTSZ,Ss,btx,stlS);
end




