% Plot depth of the maximum T>0
% 
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps/
startup

format long g
clear all
close

yr   = 2016;
iday = 60; % check winter day to avoid warm surface layer
regn = 'ARCc0.08';
expt = 110;

pfld  = 'temp';
%f_extr = 1;  % =0 - load in extracted depth of Atl. Water, =1 -extract
s_fig  = 0;

sfig=0;

rg = 9806;

pthfig  = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/%3.3i/fig_AtlLayer/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);

% Get topo:
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn] = size(LON);

Adpth=zeros(mm,nn)*nan;
nxx=0;
%mm=2520;
%nn=1600;
%ll=33;

Iarc = find(HH<=-200 & LAT>70);
nI = length(Iarc);

%if f_extr == 1
  disp('Searching for atl. water depth ...');
  pthbin = sprintf('/nexsan/archive/ARCc0.08_110/data/%4i/',yr);
 
  dnmb=datenum(yr,1,1)+iday-1;
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
  
 % Layer thickness:
%  [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
%  F=squeeze(F(plr,:,:));
%  F=F./rg;
%  F(F>1e10)=0;
%  F(F<1e-2)=0;
%  dH=squeeze(F); 

  tic;
  [F,n,m,l] = read_hycom(fina,finb,pfld);
  toc;
  F(F>1e6)=nan;
  T = F;

  Zmax = HH*0;
  Zt0  = HH*0;
  Tmax = HH*0-999;
  for i=1:nI
    i0 = Iarc(i);

    tt = squeeze(T(:,i0));
    z0 = ZZ(:,i0);
    iz0 = max(find(z0>-100));
    if ~isempty(iz0)
      tt(1:iz0) = -10;
    end
    
%    tt(1:10) = -10;
    tm = max(tt);
%    fprintf('tm=%6.4f\n',tm);
    if tm<0, continue; end;
    
    iz = find(tt==max(tt),1);
    Zmax(i0) = ZM(iz,i0);
    iz = min(find(tt>=0));

% Interpolate to find exact z
    iz1=iz-1;
    iz2=iz;
    t1 = tt(iz1);
    t2 = tt(iz2);
    z2 = ZZ(iz2,i0);
    z1 = ZZ(iz1,i0);
    t0 = 0;
    dtdz = (t2-t1)/(z2-z1);
    dz0 = (t0-t1)/(dtdz);
    zt0 = z1+dz0;
    Zt0(i0) = zt0;
    Tmax(i0) = tm;
    
    
  end
  
%end

xlim1 = 50;
xlim2 = nn-10;
ylim1 = 400;
ylim2 = mm-450;

dv = datevec(dnmb);
stl = sprintf('ARCc0.08_110, Dpth 0C m, B.Thick = 250m+/-50, %i/%2.2i/%2.2i',dv(1:3));
txtb = 'plot_atlw2.m';
pfld = 'AtlZ';
nf = 1;
Zt0(Zt0==0) = nan;
sub_plot_scalar(Zt0,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,'cmp',2);
contour(Zt0,[-500:50:-150],'k');
contour(Zt0,[-250 -250],'k','linewidth',1.6)
bottom_text(txtb,'pwd',1,'Fontsize',8);
if s_fig>0
  fgnm = sprintf('%sARCc008_%3.3i_0dgr_dpth_%i%2.2i%2.2i',pthfig,expt,dv(1:3));
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end

nf = 2; 
stl = sprintf('ARCc0.08_110, Tmax, B.Thick = 1C+/-0.5C, %i/%2.2i/%2.2i',dv(1:3));
pfld = 'AtlT';
Tmax(Tmax<-10)=nan;
sub_plot_scalar(Tmax,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
		'c1',-1,'c2',4.3,'cmp',3);
%Zmax(Zmax<-0.1)=nan;
%contour(Zmax,[-600:50:-10],'k');
contour(Tmax,[0:0.5:2],'k');
contour(Tmax,[1 1],'k','linewidth',1.6)

bottom_text(txtb,'pwd',1,'Fontsize',8);
if s_fig>0
  fgnm = sprintf('%sARCc008_%3.3i_AtlL_Tmax_%i%2.2i%2.2i',pthfig,expt,dv(1:3));
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end


%Zt0(Zt0==0)=nan;
%figure(1);
%pcolor(Zt0); shading flat;
%hold on;
%contour(HH,[0 0],'k');
%title(stt);



    








