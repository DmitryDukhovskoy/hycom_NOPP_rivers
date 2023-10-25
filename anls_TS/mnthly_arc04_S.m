% Monthly mean T/S
% Plot only 1 layer at a time
% if need to layer-average - modify the code
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.04';
%expt = 011;  
expt = 012;  

s_mat = 1; % =1 - save
plr   = 1;  % layer to plot
%plr   = 16;  % layer to plot, ~100m - similar to Claudia
%pfld  = 'temp';
pfld  = 'salin';

s_fig = 0;
rg = 9806;
fprintf('Extracting field: %s, Layer: %i\n',pfld,plr);

YRPLT=[];
cc=0;
for iyr=2005:2005
  for idd=1:5:366
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/data_mat/',expt);
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

fmat = sprintf('%sarc04_%3.3i_mnthly_%s_%i.mat',pthmat,expt,pfld,iyr);


meanS(1).Fields='monthly mean S in 1 layer';
meanS(1).Units='psu';

%inc1=1;
%inc2=1600;
%jnc1=1;
%jnc2=2520;
%djnc=2520;
%dinc=1600;

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

for im=1:12
  meanS(im).nrc=0;
  meanS(im).Favrg=HH*0;
end


% Plot fields:
cnc=0;
ip1=1;
mold=0;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
  pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i/',expt,yr);  
%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data012/%i/',yr);  % Greenland on exp
  
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
  if ~exist(fina,'file');
    fprintf('!!! Not found: %s, skipping\n\n',fina);
    continue;
  end
  
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);
  mo=DV(2);
  mday=DV(3);
  
  if mold==0; mold=mo; end;
  
  if mo~=mold
    fprintf('==== Averaging, end of month %i\n\n',mold);
    nrc=meanS(mold).nrc;
    A=meanS(mold).Favrg;
    meanS(mold).Favrg=A/nrc;
    mold=mo;
  end

  fprintf('%s, sfig=%i, V.Lr=%i, %4.4i_%2.2i_%2.2i: %s\n',...
	  pfld,s_fig,plr,DV(1:3),fina);
% in archm u_vel=utot; total velocity
%  keyboard
  f_getu=0;
  if f_getu==1
    ilr=1;
    [Fu,n,m,l] = read_hycom(fina,finb,'u-vel.');
    [Fv,n,m,l] = read_hycom(fina,finb,'v-vel.');
    U=squeeze(Fu(ilr,:,:));
    U(U>1e6)=nan;
    V=squeeze(Fv(ilr,:,:));
    V(V>1e6)=nan;
    S=sqrt(U.^2+V.^2);
  
  end

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

  meanS(mo).nrc=meanS(mo).nrc+1;
  meanS(mo).Favrg=meanS(mo).Favrg+AA;
  meanS(mo).dnmb=dnmb;
  

  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);
  mo=DV(2);
  mday=DV(3);

end;  % day loop

  
if ip==np
  fprintf('Averaging, end of loop %s\n',datestr(dnmb));
  nrc=meanS(mold).nrc;
  A=meanS(mold).Favrg;
  meanS(mold).Favrg=A/nrc;
end

fprintf('Saving %s\n',fmat);
save(fmat,'meanS','-v7.3');



