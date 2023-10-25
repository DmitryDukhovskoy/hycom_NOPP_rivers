% Plot xsections of T and S 
% annual mean - average extracted profiles
% extracted in extact TS.m

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 1;
expt = '110';
TV   = 11;  % topo version

isgm = 1; % segment to plot
yplt = 1997; % year to plot; =0 - all years

SEGM{1} = 'BeringS';
SEGM{2} = 'FramS';
SEGM{3} = 'BarOp';
SEGM{4} = 'DavisS';
nsgm = length(SEGM);
sname = SEGM{isgm};
fprintf('=====   Section:  %s  ======\n',sname);

%Zmn = 100; % average over the top Zmn m
rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%s/fig_sections/',expt);

btxt = 'plot_TSxsct.m';

fmat = sprintf('%s%s_TS_straits_%4.4i_%4.4i.mat',...
	       pthmat,expt,1993,2015);

fprintf('Loading %s\n',fmat);
load(fmat);

% HYCOM topo
ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mh,nh]=size(LON);


TM = SGMTS(isgm).TM;
DV = datevec(TM);
nrc = length(TM);
I  = SGMTS(isgm).I;
J  = SGMTS(isgm).J;
IJ = sub2ind(size(HH),J,I);
Hs = squeeze(HH(IJ));
if (I(1)==I(2))
  xsct = 0;
  XL = LAT(IJ);
else
  xsct = 1;
  XL = LON(IJ);
end


% Group by years
cc = 0;
yr_old = 0;
Iyr = [];
for j=1:nrc
  yr = DV(j,1);
  if yr~=yr_old
    cc = cc+1;
    yr_old = yr;
    Iyr(cc,1) = j;
  end
end
nyrs=cc;

yr_old = DV(1,1);
cyr = 0;
irc = 0;
Tsm = squeeze(SGMTS(isgm).Temp(1,:,:))*0;
Ssm = Tsm;
dZsm= Tsm;
[ll,mm] = size(Tsm);
YRS = [];
for it=1:nrc
  yr = DV(it,1);
  if yr~=yr_old
    cyr = cyr+1;
    Tav = Tsm./dZsm;
    Sav = Ssm./dZsm;
    dZav= dZsm/irc; % layer thickness
    TSAV(cyr).Year=yr_old;
    TSAV(cyr).Tav=Tav;
    TSAV(cyr).Sav=Sav;
    TSAV(cyr).dZ =dZav;
    YRS(cyr) = yr_old;
%keyboard
    yr_old = yr_old+1;
    Tsm = Tsm*0;
    Ssm = Ssm*0;
    dZsm= dZsm*0;
    irc = 0;
  end
  tt = squeeze(SGMTS(isgm).Temp(it,:,:));
  ss = squeeze(SGMTS(isgm).Saln(it,:,:));
  zz = squeeze(SGMTS(isgm).ZZ(it,:,:));
%  zm = squeeze(SGMTS(isgm).ZM(it,:,:));
  dz = ss*0;
  for k=1:ll
    dz(k,:)=zz(k+1,:)-zz(k,:);
    if (abs(dz)<1e-5), dz=0; end;
  end;
  Ib = find(abs(dz)<1e-5);
% For plotting - extend values below bottom
  dz(Ib)=-0.001;
%  tt(Ib)=nan;
%  ss(Ib)=nan;
  
  irc= irc+1;
  Tsm = Tsm+tt.*dz;
  Ssm = Ssm+ss.*dz;
  dZsm= dZsm+dz;
end

cyr = cyr+1;
Tav = Tsm./dZsm;
Sav = Ssm./dZsm;
dZav= dZsm/irc; % layer thickness
TSAV(cyr).Year=yr_old;
TSAV(cyr).Tav=Tav;
TSAV(cyr).Sav=Sav;
TSAV(cyr).dZ =dZav;

yr_old = yr_old+1;
Tsm = Tsm*0;
Ssm = Ssm*0;
dZsm= dZsm*0;
irc = 0;

%keyboard

% Plot mean T/S sections
%pfld='temp';
if yplt>0;
  ii1 = find(YRS==yplt);
  ii2 = ii1;
  if isempty(ii1),
    fprintf('Year %i is not found in output %i-%i\n',...
	    yplt, YRS(1), YRS(end));
  end
else
  ii1=1;
  ii2=cyr;
end
  
switch(sname)
 case('BeringS')
  dy = 5;
  mlt = mean(LAT(IJ));
  loct = sprintf('%4.2fN',mlt);
 case('FramS')
  dy = 250;
  mlt = mean(LAT(IJ));
  loct = sprintf('%4.2fN',mlt);
 case('BarOp')
  dy = 250;
  mlt = mean(LON(IJ));
  loct = sprintf('%4.2fW',mlt);
 case('DavisS')
  dy = 250;
  mlt = mean(LAT(IJ));
  loct = sprintf('%4.2fN',mlt);
end


for ii=ii1:ii2
  dZ = TSAV(ii).dZ;
  TT = TSAV(ii).Tav;
  SS = TSAV(ii).Sav;
  ZZ=TT*0;
  ZZ(end+1,:)=ZZ(end,:);
  for k=1:ll
    ZZ(k+1,:)=ZZ(k,:)+dZ(k,:);
  end
  TT(ll+1,:)=TT(ll,:);
  SS(ll+1,:)=SS(ll,:);
  [XL2,dmm]=meshgrid(XL,ZZ(:,10));
  clear dmm
  
  yr = TSAV(ii).Year;
%  switch(pfld),
%   case('temp');
%    F = TT;
%%    F2=SS;
%   case('saln');
%    F = SS;
%  end

  nfg = ii;
  pfld='temp';
  sub_plot_xsct(nfg,ZZ,XL,TT,pfld,dy,sname);
  contour(XL2,ZZ,TT,[0 0],'k','linewidth',1.2);
  contour(XL2,ZZ,TT,[0.:0.5:10],'k','linewidth',1.);
  contour(XL2,ZZ,TT,[-2:0.5:-0.1],'w','linewidth',1);
%  if strcmp(sname,'BarOp')
%    set(gca,'xdir','reverse');
%  end
  
  stl=sprintf('ARCc0.08-%s, %s, %s, %s, %i',...
	      expt,pfld,sname,loct,yr);
  title(stl);
  bottom_text(btxt,'pwd',1);

  if s_fig>0
    fgnm=sprintf('%sarc08-%s_%s_%s_%i',pthfig,expt,sname,pfld,yr);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r300',fgnm);
  end
  
  
  pfld='saln';
  nfg = ii+10;
  sub_plot_xsct(nfg,ZZ,XL,SS,pfld,dy,sname);
  switch(sname);
   case('FramS')
    contour(XL2,ZZ,SS,[34.8 34.8],'k','linewidth',1.6);
    contour(XL2,ZZ,SS,[34.9 34.9],'k');
    contour(XL2,ZZ,SS,[34.92 34.92],'k');
    contour(XL2,ZZ,SS,[34.95 34.95],'k');
   case('BarOp');
    contour(XL2,ZZ,SS,[34.8 34.8],'k','linewidth',1.6);
    contour(XL2,ZZ,SS,[34.7 34.7],'k');
    contour(XL2,ZZ,SS,[34.9 34.9],'k');
    contour(XL2,ZZ,SS,[35 35],'k');
    contour(XL2,ZZ,SS,[35.1 35.1],'k');
   case('BeringS');
    contour(XL2,ZZ,SS,[31.:0.1:33],'k','linewidth',1.);
    contour(XL2,ZZ,SS,[31.5 31.5],'k','linewidth',1.6);
   case('DavisS');
    contour(XL2,ZZ,SS,[32:0.2:35],'k','linewidth',1.);
    contour(XL2,ZZ,SS,[33 33],'k','linewidth',1.6);
  end
  
%  if strcmp(sname,'BarOp')
%    set(gca,'xdir','reverse');
%  end

  stl=sprintf('ARCc0.08-%s, %s, %s, %s, %i',...
	      expt,pfld,sname,loct,yr);
  title(stl);
  bottom_text(btxt,'pwd',1);
  
  if s_fig>0
    fgnm=sprintf('%sarc08-%s_%s_%s_%i',pthfig,expt,sname,pfld,yr);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r300',fgnm);
  end
  
  
end


  

