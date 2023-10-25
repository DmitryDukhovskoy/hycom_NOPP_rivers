% extract S fields
% calculate layer averaged and 
% Monthly mean S - saved by months/years
% S is integrated exactly over specified layers 0-50m etc 
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
%expt = 110; % no Greenland runoff  
expt = 112;  % Greenland runoff

s_mat = 1; % =2 - skip saved time, =0 - not saved
%plr   = 1;  % layer to plot
%plr   = 16;  % layer to plot, ~100m - similar to Claudia
%pfld  = 'temp';
pfld  = 'salin';

s_fig = 0;
rg = 9806;
fprintf('Expt %3.3i, Extracting field: %s\n',expt,pfld);
if s_mat==0,
  fprintf(' Mat file is not saved !!!\n');
end


% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt);


meanS(1).Fields='Layer-averaged and monthly mean S';
meanS(1).Units='psu';

LRS = load('LRS.dat');
nlrs= length(LRS)-1; % full depth not needed 

meanS(1).Layers=LRS(1:nlrs,:);

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
flnm_out = sprintf('arc08_%3.3i_mnthS_lrs_%4.4i',expt,YR);

for imo=1:12  
  for ilv=1:nlrs
    meanS(ilv).nrc=0;
    meanS(ilv).Savrg=HH*0;
  end

  for md=1:dday:31
    dnmb = datenum(yr,imo,md);
    DV = datevec(dnmb);
    mo=DV(2);
    mday=DV(3);
    if mo~=imo, continue; end; % month has <31 days
    iday = dnmb-datenum(yr,1,1)+1;

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
    
    [F,n,m,ll] = read_hycom(fina,finb,'salin');
    F(F>1e6)=nan;

% Depth average
% exactly over specified depth levels
    for ilv=1:nlrs
      zz1=LRS(ilv,1);
      zz2=LRS(ilv,2);
      zbtm=-9000;
      Sav = sub_zLayer_average(HH,ZZ,F,zbtm,zz1,zz2);      
%      keyboard
      if isempty(Sav), error('Sav is empty'), end;
% Update monthly mean:
      meanS(ilv).nrc = meanS(ilv).nrc+1;
      meanS(ilv).Savrg=meanS(ilv).Savrg+Sav;
      meanS(ilv).dnmb=dnmb;

    end   % layers    
    
    fprintf('====> Processing 1 record: %6.4f min\n',toc/60);    
  end % month days

  fprintf('Averaging, end of month %i\n',mo);
  for ilv=1:nlrs
    nrc = meanS(ilv).nrc;
    A   = meanS(ilv).Savrg/nrc;
    meanS(ilv).Savrg = A;
    zz1 = LRS(ilv,1);
    zz2 = LRS(ilv,2);
    
    fprintf('::: Lr %i %5.1f - %5.1f, min S= %6.1f, max S = %6.1f\n',...
	    ilv, abs(zz1),abs(zz2),min(min(A)),max(max(A)));
  end

  fmat = sprintf('%s%s%2.2i.mat',pthmat,flnm_out,mo);
  if s_mat>0
    fprintf('Saving %s\n',fmat);
    save(fmat,'meanS','-v7.3');
  end
  
end;  % months
end;  % Years
  



