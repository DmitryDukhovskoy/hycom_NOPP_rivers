% Show spatial variability of S
% within specifed box 
% - satellite footprint
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;

startup;

close all
clear

regn = 'ARCc0.04';
%expt = 011;  
expt = 012;  

iyr = 2005;
%mplt = [12,1,2]; % average over N months, if not, mplt=1month
mplt = [6,7,8]; % average over N months, if not, mplt=1month
nav = length(mplt);

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


yr=iyr;
iday=180;
%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i/',expt,yr);  
%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data012/%i/',yr);  % Greenland on exp

fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

dnmb=datenum(yr,1,1)+iday-1;
DV=datevec(dnmb);
mo=DV(2);
mday=DV(3);

plr=1;
[F,n,m,l] = read_hycom(fina,finb,'salin','r_layer',plr);
F(F>1e6)=nan;
SS = squeeze(F);

% Define boxes:
dbx=100e3; % approximate footprint dimensions
BX=[849        1171
         854        1171
         989         915
        1140         859
        1663        1148
        1659         694
        2017        1536];

nbx=size(BX,1);
for ibb=1:nbx
  i0=BX(ibb,1);
  j0=BX(ibb,2);
  dx0=DX(j0,i0);
  dy0=DY(j0,i0);
  nx=round(dbx/dx0);
  ny=round(dbx/dy0);
  nx2=round(nx/2);
  ny2=round(ny/2);
  h=HH(j0-ny2:j0+ny2,i0-nx2:i0+nx2);
  inn=find(isnan(h));
  if ~isempty(inn); 
    fprintf('ibb=%i, land points in the box Nland=',ibb,length(inn));
  end
  
  S0=SS(j0-ny2:j0+ny2,i0-nx2:i0+nx2);
  S0=S0(:);
  [ha,hx]=hist(S0,30);
  dhx=hx(2)-hx(1);
%  han=ha/(dhx*sum(ha)); % probability
  han=ha/length(S0);  % frequency
  
% find mean:
  Smn=han*hx';
  Smn2=mean(S0);
  sgm=std(S0);
  sg1=Smn-sgm;
  sg2=Smn+sgm;
  
  figure(ibb); clf;
  axes('Position',[0.1 0.6 0.8 0.32]);
  
  hbb=bar(hx,han);
  hold on;
  plot([Smn Smn],[0 max(han)],'r--');
  plot([Smn2 Smn2],[0 max(han)],'g--');
  plot([sg1 sg1],[0 max(han)],'k--');
  plot([sg2 sg2],[0 max(han)],'k--');
  
  set(hbb,'Facecolor',[0.6 0.6 0.6]);
  set(gca,'tickdir','out',...
	  'xtick',hx);
  
  stl=sprintf('i=%i, j=%i, S, hycom04-012, %s',i0,j0,datestr(dnmb));
  title(stl);
  
   
end

  