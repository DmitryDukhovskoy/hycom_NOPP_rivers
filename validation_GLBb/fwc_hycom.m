% Read instant output archive files:
% Note that original global HYCOM GLBb0.08 outputs 
% have been "remapped" onto ARCc0.08 grid
% see ../REMAP_ARCc/remap_gridGLB2ARC_archv
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 1; % overwritten by s_extr
s_fig  = 0;
s_extr = 1; %=0 - plot saved FWC, =1 - calculate FWC

rg=9806;  % convert pressure to depth, m
Sref=34.8;

%pthbin  = '/nexsan/GLBb0.08/';
pthbin  = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mean/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/GLBb0.08/data_mat/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/GLBb0.08/fig_climatology/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

fmat    = sprintf('%sfwc_hycom.mat',pthmat);

ftopo = sprintf('%sdepth_ARCc0.08_07.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);

% Mask of the region of interest
% exclude Pacific Ocean and North.Atl.
Lmsk            = HH;
Lmsk(Lmsk>=0)   = 0;
Lmsk(Lmsk<0)    = 1;
Lmsk(1:250,:)   = 0;
Lmsk(1935:end,:)= 0;
Iocn=find(Lmsk==1);
Ilnd=find(Lmsk==0);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% Create masks for the regions
[IH,JH] = meshgrid([1:nn],[1:mm]);
RG = regions;
IN = inpolygon(IH,JH,RG(1).IJ(:,1),RG(1).IJ(:,2));
IpC = find(IN==1); % Canada Basin
%[jj,ii]=ind2sub([mm,nn],IpC);
IN = inpolygon(IH,JH,RG(2).IJ(:,1),RG(2).IJ(:,2));
IpE = find(IN==1); % Eurasian Basin
%[jj,ii]=ind2sub([mm,nn],IpE);


if s_extr==0
%  YR=[1993,2002;2003,2012]; % make YR=[1993] to plot 1993 year
                            % YR=[yr1,yr2; yr3,y4] - will average yr1-yr2 
			    % and yr3-yr4 and plot those averages
  YR=1993;
  
  fprintf('Loading %s\n',fmat);
  load(fmat);
  sub_plot_fwc(YR,FWC,pthfig,s_fig,HH,LON,LAT);
  fprintf('Done\n');
  return
end


% Create long-term means from yearly means:
%for ik=1:2
yr1=1993;
yr2=2012;

FWCs=[];
cc=0;

%  if ik==2
%    yr1=2003;
%    yr2=2012;
%  end

for year=yr1:yr2
  fprintf('%i\n',year);
  if year<1995
    EE=190;
  else
    EE=191;
  end
  E=sprintf('GLBb0.08_%3.3i',EE);

%    pthin=sprintf('%s%s/data/meanstd/',pthbin,E);
%    fina=sprintf('%s190_archMN.%i_01_%i_12.a',pthin,year,year);
%    finb=sprintf('%s190_archMN.%i_01_%i_12.b',pthin,year,year);
  fina = sprintf('%s%3.3i_archMN_GLBb2ARCc.%i_01_%i_12.a',...
		 pthbin,EE,year,year);
  finb = sprintf('%s%3.3i_archMN_GLBb2ARCc.%i_01_%i_12.b',...
		 pthbin,EE,year,year);

  if ~exist(fina,'file');
    fprintf('doesnot exist %s\n',fina);
    continue;
  end

  cc=cc+1;
  fprintf('Reading %s\n',fina);

% Get layer thickness to be able
% to construct depth arrays of model layers
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  DP=F./rg;
  DP(DP<0.1)=nan; % 0-m layers, vanished
% Interf. depths and depths of the middle of the layers (m)
% NOTE: sign convention: depths are negative
% dP - Pa or m
  [ZZ,ZM] = sub_thck2dpth(DP); % 

  fld='salin';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  S=F;

% Get Canada Basin:
% after interpolating S into fixed z-levels
% (same as NGDM)
  for k=1:32
    dmm = squeeze(S(k,:,:));
    Slr = dmm(IpC);      % subset Canada Basin points
    ms  = nanmean(nanmean(Slr));
    m25 = prctile(Slr,0.25);
    m75 = prctile(Slr,0.75);
    S(k,1)=ms;
    S(k,2)=m25;
    S(k,3)=m75;
  end
  
    
  
% Follow ~Haine et al., 2015 definition of FWC
% The better way would be to interpolate
% between S values to find the exact depth
% Here, I simply igonre the whole layer if its S>34.8
% This is much faster and perhaps not such a big of an error
  Fwc=zeros(mm,nn);
%keyboard
  for k=1:l
    dz=abs(squeeze(DP(k,:,:)));
    dz(isnan(dz))=0;
    ss=squeeze(S(k,:,:));
    Ib=find(ss>Sref);
    ss(Ib)=Sref;
    fwc=dz.*(Sref-ss)/Sref; % m
    Fwc=Fwc+fwc;
  end
  Fwc(Ilnd)=nan; % exclude not needed regions

  fprintf('%i: BeaufGyre, FWC=%8.2f m\n',year,Fwc(1500,600)); 
  fprintf('%i: GreenGyre, FWC=%8.2f m\n',year,Fwc(800,1050)); 

  FWC(cc).Title = 'Annual mean FWC  from GLBb0.08';
  FWC(cc).Year  = year;
  FWC(cc).Fwc_m = Fwc;
% saving data
  if s_mat>0
    fprintf('Saving %s\n',fmat);
    save(fmat,'FWC');
  end;

end          % years
%end          % 10-yr cycles


