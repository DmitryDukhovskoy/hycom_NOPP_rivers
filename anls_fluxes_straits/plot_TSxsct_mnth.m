% Plot xsections of T and S 
% month mean - average extracted profiles
% extracted in extract_TS.m

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;
expt = '110';
TV   = 11;  % topo version

isgm = 1; % segment to plot
yplt = 2005; % year to plot;
mplt = 8;    % month to plot
%yplt = 2010; % year to plot;
%yplt = 2013; % year to plot;
%mplt = 9;    % month to plot

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

btxt = 'plot_TSxsct_mnth.m';

fmat = sprintf('%s%s_TS_straits_%4.4i_%4.4i.mat',...
	       pthmat,expt,1993,2016);
%fmat = sprintf('%s%s_TS_straits_%4.4i_%4.4iv2.mat',...
%	       pthmat,expt,1993,2016);

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

% Find month & year
% if d2>d1 - average
% otherwise - instant output
%d1 = datenum(yplt,mplt,1);
%d2 = d1+31; % average over 1 month
d1 = datenum(yplt,mplt,31);
d2 = d1;  % same day
if d2>d1
  Im = find(TM>=d1 & TM<=d2);
  if isempty(Im),
    fprintf('Year %i & month %i is not found in output %i-%i\n',...
	    yplt, mplt, YRS(1), YRS(end));
  end
else
  dm1=d1-15;
  dm2=d2+15;
  
  Im = find(TM>=dm1 & TM<=dm2);
  if isempty(Im),
    fprintf('Year %i & month %i is not found in output %i-%i\n',...
	    yplt, mplt, YRS(1), YRS(end));
  end

  dtm=abs(TM-d1);
  Im=find(dtm==min(dtm));
  
end

nrc = length(Im);

yr_old = DV(1,1);
Tsm = squeeze(SGMTS(isgm).Temp(1,:,:))*0;
Ssm = Tsm;
dZsm= Tsm;
[ll,mm] = size(Tsm);
irc = 0;
for ii=1:nrc
  it = Im(ii);
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

Tav = Tsm./dZsm;
Sav = Ssm./dZsm;
dZav= dZsm/irc; % layer thickness
TSAV.Year=yplt;
TSAV.Month=mplt;
TSAV.Tav=Tav;
TSAV.Sav=Sav;
TSAV.dZ =dZav;


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


dZ = TSAV.dZ;
TT = TSAV.Tav;
SS = TSAV.Sav;
ZZ=TT*0;
ZZ(end+1,:)=ZZ(end,:);
for k=1:ll
  ZZ(k+1,:)=ZZ(k,:)+dZ(k,:);
end
TT(ll+1,:)=TT(ll,:);
SS(ll+1,:)=SS(ll,:);
[XL2,dmm]=meshgrid(XL,ZZ(:,10));
clear dmm

yr = yplt;
mo = mplt;

nfg = 1;
pfld='temp';
sub_plot_xsct(nfg,ZZ,XL,TT,pfld,dy,sname);
contour(XL2,ZZ,TT,[0 0],'k','linewidth',1.2);
contour(XL2,ZZ,TT,[0.:0.5:20],'k','linewidth',1.);
contour(XL2,ZZ,TT,[-2:0.5:-0.1],'w','linewidth',1);
%  if strcmp(sname,'BarOp')
%    set(gca,'xdir','reverse');
%  end

stl=sprintf('ARCc0.08-%s, %s, %s, %s, %i',...
	    expt,pfld,sname,loct,yr);
title(stl);
bottom_text(btxt,'pwd',1);

if s_fig>0
  fgnm=sprintf('%sarc08-%s_%s_%s_%2.2i%i',pthfig,expt,sname,pfld,mo,yr);
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
  contour(XL2,ZZ,SS,[28.:0.2:33],'k','linewidth',1.);
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
  fgnm=sprintf('%sarc08-%s_%s_%s_%2.2i%i',pthfig,expt,sname,pfld,mo,yr);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end



  

