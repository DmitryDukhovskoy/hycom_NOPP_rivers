% Extract S at some locations
% 
% Monthly mean T/S averaged within the layers
% for 2 experiments (with & without Greenland runoff)
% expt_110 - no Greenland runoff
% expt 112 - with Greenland runoff
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

yr1=2005;
yr2=2005;

rg = 9806;

regn = 'ARCc0.08';
%expt = 110; % no Greenland runoff  
expt = 112;  % Greenland runoff

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);

s_mat = 0; % =2 - load and start from last saved


ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

fprintf('expt %3.3i\n\n',expt);

fprintf('Cascading %i-%i\n',yr1,yr2);
YRPLT=[];
cc=0;
for iyr=yr1:yr2
  for idd=1:5:366
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

fmat = sprintf('%sarc08_%3.3i_cascading_%i.mat',pthmat,expt,yr1);
cnc=0;
ip1=1;
mold=0;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i/',expt,yr);  
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
  if expt~=110,
    pthbin = sprintf('/nexsan/hycom/%s_%3.3i/data/%i/',regn,expt,yr);
  end
  
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);
  mo=DV(2);
  mday=DV(3);

% Layer thickness:
  [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
  F=squeeze(F);
  F=F./rg;
  F(F>1e10)=0;
  F(F<1e-2)=0;
  dH=squeeze(F); 

  tic;
  [F,n,m,l] = read_hycom(fina,finb,pfld,'r_layer',plr);
  toc;
  F(F>1e6)=nan;

  AA = squeeze(F);

end
  
  
  