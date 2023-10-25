% Extract time series of monthly Laplacian of SSH
% Laplacian <0 => Convex (buldged) surface, i.e.
% local maximum
%
% in the BG to estimate strength of BG
% similar to AOO
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_mat = 0; % =0 - do not save mat file
           % =1 - save tracer and overwrite existing mat
           % =2 - load saved and plot time series

if s_mat==0,
  fprintf('==> Extraction without saving, Mat file is not created\n');
elseif s_mat == 1
  fprintf('==> Mat file will be saved, old mat file will be overridden\n');
elseif s_mat == 2
  fprintf('==> No extraction, only plotting of previously saved mat file\n');
end


rg = 9806; 
yr1=2013;
yr2=2013;

fprintf('Years to extract: %i-%i\n',yr1,yr2);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.08';
expt = 110;  
%pthfig = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig  = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmat    = sprintf('%sLapl_SSH_BG.mat',pthmat);


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);

hmsk=HH;
hmsk(HH<0)=nan;

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;


if s_mat==2,
  sub_plot_dSSH(pthfig,fmat);
%  sub_regression_dSSH(fmat);
  return;
end

% Read fields:
ip1=1;
mold = 0;
dday = 5; 
cc=0;
for iyr = yr1:yr2
  yr=iyr;
%  for iday = 337:dday:366
%  id1=336;
  id1=2;
%  if yr == 2011,
%    id1 = datenum(2011,10,1)-datenum(2011,1,1)+1;
%  end
  
  for iday = id1:dday:365
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);

    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    imo=DV(2);
    fprintf('Processing %i/%2.2i/%2.2i\n',DV(1:3));

    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    cc=cc+1;
    
    tic;
    [F,nn,mm,ll] = read_hycom(fina,finb,'srfhgt');
    toc;
    F(F>1e6)=nan;
    E=squeeze(F/9.806);
    E(HH>-500)=nan;
    E(1900:end,:)=nan;
    E=E-nanmean(nanmean(E));
    ik1=500;
    ik2=650;
    jk1=1500;
    jk2=1630;
    dmm=E;
    dmm(HH>-2000)=nan;
    AA=dmm(jk1:jk2,ik1:ik2);
    amx=max(max(AA));
    [ja,ia]=find(AA==amx,1);
    j0=jk1+ja-1;
    i0=ik1+ia-1;
%   keyboard
 % Use Laplacian - should be <0 - convex shape   
    di=80;
    e1=E(j0,i0-di);
    e2=E(j0,i0+di);
    pp=0;
    while (isnan(e1) | isnan(e2)),
      pp=pp+1;
      di=round(0.9*di);
      e1=E(j0,i0-di);
      e2=E(j0,i0+di);
      if pp>10, error('endless loop, BG max SSH ...'); end;
    end  
    e0=E(j0,i0);
    dx = sum(DX(j0,i0-di));
    d2Edx2=(e2-2*e0+e1)/(dx^2);
    
    dj=80;
    e1=E(j0-dj,i0);
    e2=E(j0+dj,i0);
    pp=0;
    while (isnan(e1) | isnan(e2)),
      pp=pp+1;
      dj=round(0.9*dj);
      e1=E(j0-dj,i0);
      e2=E(j0+dj,i0);
      if pp>10, error('endless loop dj, BG max SSH ...'); end;
    end  
    dy = sum(DY(j0-dj:j0+dj,i0));
    d2Edy2=(e2-2*e0+e1)/(dy^2);
    Le=d2Edx2+d2Edy2;
    
    f_plt=0;
    if f_plt>0
      figure(10); clf;
      pcolor(E); shading flat;
      caxis([-0.4 0.4]);
      hold on;
      contour(HH,[0 0],'k');
      colorbar
      contour(E,[0.0:0.05:0.5],'r');
    end
    
      
    
    if Le>0,
      fprintf('!!!!    Le>0 %7.4d\n',Le);
      keyboard;
    end
    
    fprintf('Lapl.E=%6.3d\n',Le);

    
    LAPLE.TM(cc,1)=dnmb;
    LAPLE.d2E_BG(cc,1)=Le;
    
  end
  
end

if s_mat>0
  fprintf('Saving %s\n',fmat);
  save(fmat,'LAPLE');
end

