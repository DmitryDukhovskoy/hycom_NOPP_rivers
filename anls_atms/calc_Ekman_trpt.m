% Calculate Ekman transport - 
% same as upwelling index
% using CCMp winds
% extraceted at
% several locations
% on the Greenland shelf
%
% along Greenland coast
% using Bakun method
% M = 1/(f*rho_w)*(tau x k)
% see also Picket and Padun,
% upwelling Cal. Current, JGR
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

rhoa = 1.2;
Cd   = 0.0013;
rhow = 1027;

%s_mat = 1;

pthdat = '/Net/data/ccmp/v02.0/';
pthdt  = '/Net/tholia/ddmitry/ovwst/data_mat/';
pth8='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthdat1 = '/Net/kronos/ddmitry/ncep_cfsr/';
%pthdat2 = '/Net/kronos/ddmitry/ncep_cfsv2/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';

%pthout  = ('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/');
%fmat    = sprintf('%sarc008_Greenl_upwl_CFSR_month.mat',pthout);
fmat = sprintf('%swrose_Greenland_ccmp020rss_1997_2016.mat',pthdt);
btx='calc_Ekman_trpt.m';

% Get indices for interpolating onto HYCOM ARCc0.72 grid:
%fgrd = sprintf('%sccmp_gridindx_arc008.mat',pth8);
%INDX=load(fgrd);

% Topo for ARCc0.08
ftopo = sprintf('%s/depth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
%[DX,DY]=sub_dx_dy(LN,LT);
%[II,JJ] = meshgrid([1:nn],[1:mm]);
%Fcor = 2*7.292e-5*sind(LAT);

fprintf('Loading %s\n',fmat);
load(fmat); % WCLIM WRS

nbx=length(WRS);

%
% Plot boxes
f_rgnmap=0;
if f_rgnmap==1
  figure(1); clf;
  contour(HH,[0 0],'Color',[0.4 0.4 0.4]);
  hold on;
  contour(HH,[-500 -500],'Color',[0.8 0.8 0.8]);

  nrg=length(WRS);
  for irg=1:nrg
    xv=WRS(irg).Box_VX;
    yv=WRS(irg).Box_VY;
    [IC,JC] = sub_find_indxHYCOM(LON,LAT,xv,yv);
    plot(IC,JC,'r.-');
    text(IC(1),JC(1),sprintf('%i',irg),'Fontsize',12);
  end
  axis('equal');
  set(gca,'xlim',[450 1000],...
	  'ylim',[300 1100]);
  title('Wind Roses Averaged Regions');
  
  contour(LON,[-110:10:0],'Color',[0.8 0.8 0.8]);
  contour(LAT,[50:10:89],'Color',[0.8 0.8 0.8]);
  
  bottom_text(btx,'pwd',1);
end

% Find general directions 
% of the coast for each box
% by specifing endponts of
% the vector 
% c/clckwise direction
% around Greenland such that
% positive wind projection = upwelling
% Indices are for 0.08 grid! 
Dr(1).XY=[924   787
   949   961];
Dr(2).XY=[882   645
   922   796];
Dr(3).XY=[753   559
   878   645];
Dr(4).XY=[656   500
   789   586];
Dr(5).XY=[617   415
   651   497];
Dr(6).XY=[581   424
	 619   415];
Dr(7).XY=[535   537
	  565   440];
Dr(8).XY=[535   660
	  538   518];
Dr(9).XY=[582   780
	  525   638];
Dr(10).XY=[609   882
	   570   756];
Dr(11).XY=[581   957
	   614   914];

% Calculate orientation of the
% local isobath directions
% wrt to 0 degree = Eastward axis
for ik=1:length(Dr)
  i1=Dr(ik).XY(1,1);
  i2=Dr(ik).XY(2,1);
  j1=Dr(ik).XY(1,2);
  j2=Dr(ik).XY(2,2);
  ln1=LON(j1,i1);
  ln2=LON(j2,i2);
  lt1=LAT(j1,i1);
  lt2=LAT(j2,i2);
  L=distance_spheric_coord(lt2,ln2,lt1,ln1);
  Lx=dsphc(lt1,ln2,lt1,ln1);
  Lx=sign(ln2-ln1)*Lx;
  Ly=dsphc(lt2,ln1,lt1,ln1);
  Ly=sign(lt2-lt1)*Ly;
  
 % plot([i1 i2],[j1 j2],'-');
%
% Vector in 0 direction - East:
%  Nx=1;
%  Ny=0;
  thtd=atan2d(Ly,Lx);
  Dr(ik).Theta_dgr=thtd;
  
% Approximate segment length btw the boxes  
  x0=WRS(ik).Box_XY(1);
  y0=WRS(ik).Box_XY(2);
  if ik==1,
    x1=350;
    y1=82;
    d1=dsphc(y1,x1,y0,x0);
  else
    x1=WRS(ik-1).Box_XY(1);
    y1=WRS(ik-1).Box_XY(2);
    d1=0.5*dsphc(y1,x1,y0,x0);
  end
  
  if ik==length(Dr)
    x2=290;
    y2=75;
    d2=dsphc(y2,x2,y0,x0);
  else
    x2=WRS(ik+1).Box_XY(1);
    y2=WRS(ik+1).Box_XY(2);
    d2=0.5*dsphc(y2,x2,y0,x0);
  end

  D=d1+d2;
  Dr(ik).Segment_Length=D;
  
end  
  
  
% Calc Ekman transport
for ik=1:nbx
  U=WRS(ik).U;
  V=WRS(ik).V;
%
  S = sqrt(U.^2+V.^2);
  
  lat0=WRS(ik).Box_XY(2);
  Fcor = 2*7.292e-5*sind(lat0);

% Upwelling index, Ekman transport
% m3/s
% Uwelling/downwelling transport
% is obtained by projecting the Ekman 
% transp vector onto a vector
% aligned with the bathymetric gradient
% positive toward deeper water
% Wind-drag - Large and Pond 1981
  Cd=U*0+1.2e-3; % |U|<11 m/s
  I=find(S>=11 & S<25);
  Cd(I)=(0.49+0.065*S(I))*1e-3;
  I=find(S>=25);
  Cd(I)=2.115*1e-3;
  
  tx = rhoa*Cd.*S.*U;
  ty = rhoa*Cd.*S.*V;
  Mx = 1./(rhow*Fcor).*ty;
  My = -1./(rhow.*Fcor).*tx;

% Upwelling favorable winds are negative
% Negative Ekman transp projected on bath grad
% is Downwelling 
% If c/clckwise direction of isobath orientation
% then positive depth grad is +90 dgr
  thtD=Dr(ik).Theta_dgr-90;
  Pcst=[cosd(thtD), sind(thtD)];
  Wcst=[Mx,My]*Pcst'; % m3/s per 1 m of coast line
  
  ETRP(ik).Transp_m3sec=Wcst;
%keyboard

end

for ik=1:nbx
  U=WRS(ik).U;
  V=WRS(ik).V;
%
  S = sqrt(U.^2+V.^2);
  
  lat0=WRS(ik).Box_XY(2);
  Fcor = 2*7.292e-5*sind(lat0);

% Upwelling index, Ekman transport
% m3/s
% Uwelling/downwelling transport
% is obtained by projecting the Ekman 
% transp vector onto a vector
% aligned with the bathymetric gradient
% positive toward deeper water
% Wind-drag - Large and Pond 1981
  Cd=U*0+1.2e-3; % |U|<11 m/s
  I=find(S>=11 & S<25);
  Cd(I)=(0.49+0.065*S(I))*1e-3;
  I=find(S>=25);
  Cd(I)=2.115*1e-3;
  
  tx = rhoa*Cd.*S.*U;
  ty = rhoa*Cd.*S.*V;
  Mx = 1./(rhow*Fcor).*ty;
  My = -1./(rhow.*Fcor).*tx;

% Upwelling favorable winds are negative
% Negative Ekman transp projected on bath grad
% is Downwelling 
% If c/clckwise direction of isobath orientation
% then positive depth grad is +90 dgr
  thtD=Dr(ik).Theta_dgr-90;
  Pcst=[cosd(thtD), sind(thtD)];
  Wcst=[Mx,My]*Pcst'; % m3/s per 1 m of coast line
  
  ETRP(ik).Transp_m3sec=Wcst;
  
end

TM=WRS(1).TM;
DV=datevec(TM);
yr1=DV(1,1);
yr2=DV(end,1);
%
%
% Plotting EKman Transport
%
for ik=1:nbx;
  Wcst=ETRP(ik).Transp_m3sec;

  YRS=[];
  for iyr=yr1:yr2
    J=find(DV(:,1)==iyr);
    nj=length(J);
    y0=DV(J(1));
    dmm=[y0:1/nj:y0+(nj-1)/nj];
    YRS=[YRS,dmm];
  end
  YRS=YRS(:);

  YM=[];  % time array for plotting means
  cc=0;
  for iyr=yr1:yr2
    for im=1:12
      I=find(DV(:,1)==iyr & DV(:,2)==im);
      wm=mean(Wcst(I));
      cc=cc+1;
      WcstM(cc)=wm;
      dmm=YRS(I);
      YM(cc,1)=dmm(1);
      YM(cc,2)=dmm(end);
    end
  end

  % Convert to km3/mo per 10 km for monthly means
  Lc=10e3;
  Wcst_mo=WcstM*3600*24*30/1e9*Lc;

  nyr=length(Wcst_mo)/12;
  Wcst_mo=reshape(Wcst_mo,[12,nyr]);
  Wcst_clm=mean(Wcst_mo,2);

  di=0.1;
  figure(ik); clf;
  axes('position',[0.1 0.5 0.6 0.4]);
  hb=bar(Wcst_clm);
  hold on;
  for im=1:12
    y1=prctile(Wcst_mo(im,:),10);
    y2=prctile(Wcst_mo(im,:),90);
    plot([im im],[y1 y2],'Color',[0 0 0],'Linewidth',1.8);
    plot([im-di im+di],[y1 y1],'Color',[0 0 0],'Linewidth',1.8);
    plot([im-di im+di],[y2 y2],'Color',[0 0 0],'Linewidth',1.8);
  end

  set(gca,'tickdir','out',...
	  'xlim',[0.5 12.5],...
	  'xtick',[1:12],...
	  'Fontsize',14,...
	  'xgrid','on',...
	  'ygrid','on');

  stl=sprintf('Bx %i, CCMP Ekman Trnspt, km^3/mo / %i km',ik,Lc*1e-3);
  title(stl,'Fontsize',12);

% Cumulative Ekm trpt:
  Ftrp=cumsum(Wcst_mo,1)/Lc*Dr(ik).Segment_Length; % km3
  mF=mean(Ftrp,2);
  for im=1:12
    prc1(im)=prctile(Ftrp(im,:),10);
    prc2(im)=prctile(Ftrp(im,:),90);
  end
  
  axes('position',[0.1 0.1 0.6 0.3]);
  plot(mF,'Color',[0 0.2 0.8],'Linewidth',2);
  hold on;
  plot(prc1,'Color',[0.8 0.9 1],'Linewidth',2);
  plot(prc2,'Color',[0.8 0.9 1],'Linewidth',2);
  yl1=1.05*min([min(prc1),0]);
  yl2=max([max(prc2), 0]);
  dyy=round(round(ceil(yl2-yl1)/100)*100/8);
  
  set(gca,'tickdir','out',...
	  'xlim',[0.5 12.5],...
	  'xtick',[1:12],...
	  'ylim',[yl1 yl2],...
 	  'Fontsize',14,...
	  'xgrid','on',...
	  'ygrid','on');
%	  'ytick',[-10*dyy:dyy:10*dyy],...
  
  axes('position',[0.78 0.2 0.2 0.1]);
  stl=sprintf('Cumulative Ekm trpt over segm \n btw boxes, km^3, 10,90 prctl');
  text(0,0,stl);
  set(gca,'visible','off');
  
  bottom_text(btx,'pwd',1);
  
%  set(gcf,'Position',[933 35 1437 1297]);
  drawnow
end

