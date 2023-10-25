% Plot extracted T/S fields
% Save and Plot trends across the specified years
% by layers annual averaged and 
% Annual mean T/S - saved by years
% S is integrated exactly over specified layers 0-50m, 50-200m, etc 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 1993;
YR2 = 2009;

regn = 'ARCc0.08';
expt = 110; % no Greenland runoff  
%expt = 112;  % Greenland runoff

pfld  = 'temp';
%pfld  = 'salin';

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt);

% New experiments GOFS3.X use Topo 11
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
%[mm,nn]=size(LON);

% Cut the region
i1=460;
j1=320;
i2=1010;
j2=1140;
H=HH(j1:j2,i1:i2);
LN=LON(j1:j2,i1:i2);
LT=LAT(j1:j2,i1:i2);
[mm,nn]=size(LN);

LRS=[0, -50; ...
     -50, -200; ...
     -200, -800];

%LRS=[0, -50; ...
%     -50, -200];

icc=0;
for YR=YR1:YR2
  yr = YR;
  flnm_out = sprintf('arc08_%3.3i_yr-%s_lrs_%4.4i',expt,pfld,YR);
  fmat = sprintf('%s%s.mat',pthmat,flnm_out);

  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  nlv=length(meanF);
  icc=icc+1;
  for ilv=1:nlv
    A=meanF(ilv).Savrg(j1:j2,i1:i2);
    F(ilv).A(icc,:,:)=A;
  end
end

for ilv=1:nlv
  A=F(ilv).A;
  Sm=squeeze(nanmean(A,1));
  I=find(~isnan(Sm));
  
  fprintf('Finding trends ...\n');
  Tr=Sm*nan;
  for ik=1:length(I);
    if mod(ik,10000)==0,
      fprintf(' ... Done %3.1f%%\n',ik/length(I)*100);
    end
    
    i0=I(ik);
    [jj,ii]=ind2sub(size(Sm),i0);
    Y=A(:,jj,ii);
    X=[Y*0+1,[1:length(Y)]'];
    B=regress(Y,X);
    Tr(i0)=B(2); % trend deg/yr
  end
  
  FMEAN(ilv).Depth=H;
  FMEAN(ilv).LON=LN;
  FMEAN(ilv).LAT=LT;
  FMEAN(ilv).Layer=LRS(ilv,:);
  FMEAN(ilv).Years(1:2)=[YR1,YR2];
  FMEAN(ilv).mean_Field=Sm;
  FMEAN(ilv).trend=Tr;
  
end

f_plt=1;
btx='plot_annual_layers_TS_arc008.m';
if f_plt==1
  for ilv=1:nlv
    T=FMEAN(ilv).mean_Field;
    zz=FMEAN(ilv).Layer;
    figure(ilv); clf;
    contour(H,[0 0],'k'); hold on;
    pcolor(T); shading flat;
    colorbar;
    axis('equal');
    
    stl=sprintf('008-%3.3i, mean %s, %i-%i, %3.1f-%3.1f',...
		expt,pfld,YR1,YR2,zz);
    title(stl);
  bottom_text(btx,'pwd',1);
  end
  
  c1=-0.3;
  c2=0.3;
  D=create_colormap3(200,c1,c2);
  cmp=D.colormap;
  for ilv=1:nlv
    T=FMEAN(ilv).trend;
    zz=FMEAN(ilv).Layer;
    figure(ilv+10); clf;
    contour(H,[0 0],'k'); hold on;
    pcolor(T); shading flat;
    caxis([c1 c2]);
    colormap(cmp);
    colorbar;
    axis('equal');
    
    stl=sprintf('008-%3.3i, d%s/dT (per yr), %i-%i, %3.1f-%3.1f',...
		expt,pfld,YR1,YR2,zz);
    title(stl);
  bottom_text(btx,'pwd',1);
  end
  
end;  

f_fmat=1;
fmatout=sprintf('%sarc008_%3.3i_mean_trend_%s_%i-%i.mat',pthmat,expt,pfld,YR1,YR2);
if f_fmat==1
  fprintf('Saving %s\n',fmatout);
  save(fmatout,'FMEAN');
end


    
    
    