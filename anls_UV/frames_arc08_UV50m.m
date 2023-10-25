% Extract and calculate Monthly mean U,V fields
% averaged over exact depth levels zz1:zz2
%
% see anls_Greenland/monthly_mean_flds.m
% exact averaging: extr_MassTrcr_month.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 1993;
YR2 = 1998;
dday= 7;
plr = 11;  % ~50m deep ocean

s_frm = 1; % =0 - do not save


% Averaged over the depth layers:
%zz1=0;
%zz2=-50;

regn = 'ARCc0.08';
expt = 110;
rg = 9806;

pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/frames_UV';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);

fprintf('Monthly mean UV in the layer: %i-%im\n',zz1,zz2);


rg=9806;  % convert pressure to depth, m

%figure(1); clf;
%set(gcf,'Visible','off');

%YPLT=[2004,299;2005,351;2006,295;2007,354;2008,286;2009,236];
% Animation:
YRPLT=[];
cc=0;
dday = 7;
for iyr=YR1:YR2
  for idd=1:dday:365
%    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);


ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);


hmsk=HH;
hmsk(HH<0)=nan;

if ~exist('DX','var')
  [DX,DY]=sub_dx_dy(LON,LAT);
end


for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  fmat = sprintf('%smnthUV_%4.4i-%4.4i_%i.mat',pthmat,abs(zz1),abs(zz2),yr);
  
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);
  imo=DV(2);
  fprintf('Processing Month %i\n',imo);
  
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);

  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  fprintf('%4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

  tic;
  [F,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',plr);
  F(F>1e6)=nan;
  U=squeeze(F);
  
  [F,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',plr);
  F(F>1e6)=nan;
  V=squeeze(F);
  V(HH>-52)=nan; % near coast 
  
%  fprintf('Getting layer depths ...\n');
%  [ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',1);
%  
%  dz=abs(zz2-zz1);
%  zbtm=-9000;
%  Uav = sub_zLayer_average(HH,ZZ,U,zbtm,zz1,zz2); % depth-average field 
%  Vav = sub_zLayer_average(HH,ZZ,V,zbtm,zz1,zz2); % depth-average field 
  
  %
  nf=1;
  stl=sprintf('U %4.4i/%2.2i/%2.2i, at 50m',DV(1:3));
  u=U;
  v=V;
  s=sqrt(u.^2+v.^2);
  c1=0;
  c2=0.5;
  xl1=200;
  xl2=1300;
  yl1=50;
  yl2=1200;
  sub_plot_divu(s,LON,LAT,nf,HH,u,v,stl,xl1,xl2,yl1,yl2,c1,c2);


 
end;  % time loop

if s_mat>0
  fprintf('===End:   Saving meanUV month %i ...\n',mold);

%  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,yr);
  fprintf('Saving %s\n',fmat);
  save(fmat,'meanUV');
end




