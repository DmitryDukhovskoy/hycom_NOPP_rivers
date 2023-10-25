% Area-mean MLD
% over specified regions
% Output from the model
% time series extracted in extr_mld.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

yr1=1993;
yr2=2016;
nbx = 5;

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

btx = 'plot_mld_tser.m';

% Define Regions:
f_pltbox=0;
BX = sub_define_boxes(HH,LON,LAT,f_pltbox);

[XX,YY] = meshgrid((1:nn),(1:mm));
for ib=1:nbx
  iBG = BX(ib).IJ;
  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
  IN = find(INP==1);
  BX(ib).IN = IN;
end

cc=0;
for iyr=yr1:yr2
  tic;
  fmat = sprintf('%smld_%i.mat',pthmat,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  for imo=1:12
    fprintf('%i/%2.2i\n',iyr,imo);
    F = MLD(imo).mld;
    cc=cc+1;
    for ib=1:nbx % regions
      IN=BX(ib).IN;
      Arg = sum(Acell(IN)); % total area, region
      mld = nansum(F(IN).*Acell(IN))./Arg; % avrg Tr conc in region
      BX(ib).TM(cc)   = datenum(iyr,imo,15);
      BX(ib).mld(cc) = mld; % avrg Tr conc in region
    end
  end
  
  fprintf('1 yr processed %6.3fmin\n',toc/60);
  
end

% =================
% Plot time series 
% of heat fluxes
% =================
TM  = BX(1).TM;
nrc = length(TM); 
yrs = [0:nrc-1]/12+yr1;

POS = [0.08 0.71 0.85 0.22; ...
       0.08 0.4 0.85 0.22; ...
       0.08 0.1 0.85 0.22];


for ib=1:nbx
  nr=mod(ib,3);
  if (nr==0), nr=3; end;
  pos=POS(nr,:);
  nk=floor(ib/3)+10;
  
  mld = BX(ib).mld;
  nm  = BX(ib).Name;
  yl1=1.1*(min(mld));
  yl2=0;
  
  if (nr==1), 
    figure(nk); clf; 
  end;
  axes('position',pos);
  plot(yrs,mld);
  set(gca,'tickdir','out',...
	  'ylim',[yl1 yl2],...
	  'xlim',[yr1 max(yrs)],...
	  'xtick',[yr1:yr2],...
	  'xgrid','on',...
	  'ygrid','on',...
	  'xminortick','on',...
	  'yminortick','on');
  sll=sprintf('%s, Monthly MLD, m',nm);
  title(sll);
  
  if (nr==3); bottom_text(btx,'pwd',1); end;
end
if (nr==2); bottom_text(btx,'pwd',1); end;

  