% Plot difference of 
% Extract S fields and calc mean over specified days for given year 
% Try to compare anomalies from the 2 experiments
% by removing seasonality based on 1993
% instead of direct S fields
% 
% Monthly mean S is averaged within the layers
% and exactly within the specified depth layers
% for 2 experiments (with & without Greenland runoff)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1  = 2012;
YR2  = 2016;
imo  = 13; % 13 - yearly average

LRS = load('LRS.dat');
nlrs= length(LRS)-1; % skip whole depth  
       
regn = 'ARCc0.08';
expt = 110; % no Greenland runoff  
%expt = 112;  % Greenland runoff

s_fig = 0;

fprintf('Years: %i-%i\n',YR1,YR2);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

%fmat = sprintf('%sarc08_%3.3i_mnthly_%s_%i.mat',pthmat,expt,pfld,iyr);

%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

%[DX,DY]=sub_dx_dy(LON,LAT);
%Acell=DX.*DY; % Grid cell area, m2

mo=imo;
if imo==13
  mo1=1;
  mo2=12;
else
  mo1=imo;
  mo2=imo;
end

% Get reference years;
YR0=1993;
cc=0;
for iyr=YR1:YR2
  YR=iyr;
  for mo=mo1:mo2
    cc=cc+1;
    expt1=110;
    pthmat1  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt1);
    pthm1 = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt1);
% Reference field:
%    fmat1 = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthm1,expt1,YR0,mo);
%    fprintf('Loading %s\n',fmat1);
%    load(fmat1);
%    S1r=meanS;

    fmat1 = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthm1,expt1,YR,mo);
    fprintf('Loading %s\n',fmat1);
    load(fmat1);
    if cc==1
      S1=meanS;
    else
      for ilv=1:nlrs
	dmm=meanS(ilv).Savrg;
	S1(ilv).Savrg=S1(ilv).Savrg+dmm;
      end
    end


    expt2=112;
    pthmat2  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt2);
    pthm2 = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt2);
% Reference field:
%    fmat2 = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthm2,expt2,YR0,mo);
%    fprintf(':: Loading %s\n',fmat2);
%    load(fmat2);
%    S2r=meanS;
    
    fmat2 = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthm2,expt2,YR,mo);
    fprintf(':: Loading %s\n',fmat2);
    load(fmat2);
    if cc==1
      S2=meanS;
    else
      for ilv=1:nlrs
	dmm=meanS(ilv).Savrg;
	S2(ilv).Savrg=S2(ilv).Savrg+dmm;
      end
    end
%keyboard
  end % month
end % year

if cc>1
  for ilv=1:nlrs
    dmm=S1(ilv).Savrg;
    S1(ilv).Savrg=dmm/cc;
    
    dmm=S2(ilv).Savrg;
    S2(ilv).Savrg=dmm/cc;
  end
end


%cl1=colormap_blue(200);
cl1 = colormap_cold(200);
%cl2 = colormap_red(200);
%cl3 = colormap_purple(100);
for ik=1:10;
  cl1(ik,:)=[1 1 1];
%  cl3(ik,:)=[1 1 1];
end
cl1=flipud(cl1);
%cl3=flipud(cl3);
%cmp = [cl1;cl3];
cmp=cl1;
cmp = smooth_colormap(cmp,10);
cmp = smooth_colormap(cmp,10);

c1 = -0.1;
c2 = 0.;


xl1=300;
xl2=1300;
yl1=200;
yl2=1100;

Lmsk = HH*0;
Lmsk(HH<0)=1;
%lmp=[0 0 0; 0.8 .8 .8];
lmp=[0.4 0.4 0.4; 0.9 .9 .9];

for ilr=1:nlrs
  fprintf('Plotting Lr %i\n',ilr);
  zz1 = LRS(ilr,1);
  zz2 = LRS(ilr,2);
  F1 = S1(ilr).Savrg;
  F2 = S2(ilr).Savrg;
  F1(F1==0)=nan;
  F2(F2==0)=nan;
  dF = F2-F1;
  
  sgmx = 4;
  sgmy = sgmx;
  npnts = 11;
  dFf = sub_gauss_filter(dF,sgmx,sgmy,npnts,xl1,xl2,yl1,yl2);
  dF = dFf;
  
  zb2 = LRS(ilr,2);
%  dF(HH>zb2)=nan;
  dF(HH>zz1)=nan;
  
  figure(ilr); clf;
  pcolor(Lmsk); shading flat;
  colormap(lmp);
  freezeColors;
  hold on;
  
  pcolor(dF); shading flat;
  caxis([c1 c2]);
%  hold on;
%  contour(HH,[0 0],'w','linewidth',1.2);
  axis('equal');
  set(gca,'xlim',[xl1 xl2],...
	  'ylim',[yl1 yl2],...
	  'Color',[0.8 0.8 0.8],...
	  'xtick',[],...
	  'ytick',[]);
  colormap(cmp);
  chb=colorbar('location','eastoutside');
  set(chb,'Fontsize',14,'ticklength',0.02);
  
  if imo==13, imo=0; end;
  ctl=sprintf('dltS=Gr-NoGr, %4.1f-%4.1f m, %4.4i-%4.4i',abs(zz1),abs(zz2),YR1,YR2);
  title(ctl);
  
  btx='dltS_LrAvrg_GreenldExp_v2.m';
  bottom_text(btx,'pwd',1);
end


