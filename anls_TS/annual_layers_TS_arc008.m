% extract T/S fields
% calculate layer annual averaged and 
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

YR1 = 2016;
YR2 = 2016;



dday= 8; 
regn = 'ARCc0.08';
expt = 110; % no Greenland runoff  
%expt = 112;  % Greenland runoff

s_mat = 1; % =2 - skip saved time, =0 - not saved
%plr   = 1;  % layer to plot
%plr   = 16;  % layer to plot, ~100m - similar to Claudia
pfld  = 'temp';
%pfld  = 'salin';

s_fig = 0;
rg = 9806;
fprintf('Expt %3.3i, Extracting field: %s, %i-%i\n',expt,pfld,YR1,YR2);
if s_mat==0,
  fprintf(' Mat file is not saved !!!\n');
end


% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt);

meanF(1).Fields=sprintf('Layer-averaged and monthly mean %s',pfld);

%LRS = load('LRS.dat');
LRS=[0, -50; ...
     -50, -200; ...
     -200, -800];

nlrs= length(LRS); 

%meanF(1).Layers=LRS(1:nlrs,:);

% New experiments GOFS3.X use Topo 11
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


for YR=YR1:YR2
  yr = YR;
  flnm_out = sprintf('arc08_%3.3i_yr-%s_lrs_%4.4i',expt,pfld,YR);

  meanF=struct;
  for ilv=1:nlrs
    meanF(ilv).nrc=0;
    meanF(ilv).Savrg=HH*0;
    meanF(ilv).Layer(:)=LRS(ilv,:);
  end
  
  for iday=1:dday:365
    dj1=datenum(YR,1,1);
    dnmb=dj1+iday-1;
    DV=datevec(dnmb);
    
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
    if expt~=110,
      pthbin = sprintf('/nexsan/hycom/%s_%3.3i/data/%i/',regn,expt,yr);
    end

  %  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data012/%i/',yr);  % Greenland on exp
    tic; 
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end
  
    fprintf('%s, %4.4i_%2.2i_%2.2i: %s\n',pfld,DV(1:3),fina);

 % Layer thickness:
    [ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',1);
    
    [F,n,m,ll] = read_hycom(fina,finb,pfld);
    F(F>1e6)=nan;

% Depth average
% exactly over specified depth levels
    for ilv=1:nlrs
      zz1=LRS(ilv,1);
      zz2=LRS(ilv,end);
      zbtm=-9000;
      Sav = sub_zLayer_average(HH,ZZ,F,zbtm,zz1,zz2);      
%      keyboard
      if isempty(Sav), error('Sav is empty'), end;
% Update mean:
      meanF(ilv).nrc = meanF(ilv).nrc+1;
      meanF(ilv).Savrg=meanF(ilv).Savrg+Sav;
      meanF(ilv).dnmb=dnmb;

    end   % layers    
    
    fprintf('====> Processing 1 record: %6.4f min\n',toc/60);    
  end % year days

  fprintf('Averaging, end of year %i\n',YR);
  for ilv=1:nlrs
    nrc = meanF(ilv).nrc;
    A   = meanF(ilv).Savrg/nrc;
    meanF(ilv).Savrg = A;
    zz1 = LRS(ilv,1);
    zz2 = LRS(ilv,2);
    
    fprintf('::: Lr %i %5.1f - %5.1f, min S= %6.1f, max S = %6.1f\n',...
	    ilv, abs(zz1),abs(zz2),min(min(A)),max(max(A)));
  end

  fmat = sprintf('%s%s.mat',pthmat,flnm_out);
  if s_mat>0
    fprintf('Saving %s\n',fmat);
    save(fmat,'meanF');
  end
  

end;  % Years
  



