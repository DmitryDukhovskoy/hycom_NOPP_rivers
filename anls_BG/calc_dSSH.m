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

s_mat = 1; % =0 - do not save mat file
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
yr1=1993;
yr2=2016;

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
[II,JJ]=meshgrid([1:nn],[1:mm]);

hmsk=HH;
hmsk(HH<0)=nan;

% Define interior Arctic Ocean
hmin = -800;
ARC = sub_arctic_domain(HH,hmin);

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;

if s_mat==2,
  fprintf('Plot results here: plot_dSSH.m');
%  sub_plot_dSSH(pthfig,fmat);
%  sub_regression_dSSH(fmat);
  return;
end

% Read fields:
ip1=1;
mold = 0;
dday = 6; 
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
    
    Inan = ARC.Inan;
    E(Inan)=nan;
    E=E-nanmean(nanmean(E));
    E(Inan) = 0;
%    E(HH>0) = nan;

% Find mean d2(SSH) within this box    
    ik1=488; % BG 
    ik2=670;
    jk1=1480;
    jk2=1665;
    iz1=948; % Eurasian basin
    iz2=1152;
    jz1=1012;
    jz2=1372;
    if ~exist('Irg','var')
      dmm = inpolygon(II,JJ,[ik1, ik1, ik2, ik2],[jk1, jk2, jk2, jk1]);
      Irg = find(dmm==1);
      nirg= length(Irg);
      dmm = inpolygon(II,JJ,[iz1, iz1, iz2, iz2],[jz1, jz2, jz2, jz1]);
      Izg = find(dmm==1);
      nizg= length(Izg);
    end
    Emn=mean(E(Irg));
    Emx=max(E(Irg));
    Ezmn=mean(E(Izg));
    Ezmin=min(E(Izg));
    
    
% Laplacian - too noisy: picking up small anomalies    
%    di=20;
%    dj=di;
%    crr=0;
%    clear LE
%    for irr=1:nirg
%      ix = Irg(irr);
%      [j0,i0]=ind2sub(size(HH),ix);
%      e1=E(j0,i0-di);
%      e2=E(j0,i0+di);
%      e0=E(j0,i0);
%      dx = sum(DX(j0,i0-di));
%      d2Edx2=(e2-2*e0+e1)/(dx^2);
%      e1=E(j0-dj,i0);
%      e2=E(j0+dj,i0);
%      dy = sum(DY(j0-dj:j0+dj,i0));
%      d2Edy2=(e2-2*e0+e1)/(dy^2);
%      crr=crr+1;
%      LE(crr)=d2Edx2+d2Edy2;
%    end   
%    dmm = HH*nan;
%    dmm(Irg) = LE;
%    pcolor(dmm); shading flat

    
    f_plt=0;
    if f_plt>0
      fprintf('Plotting ...\n');
      figure(10); clf;
      pcolor(E); shading flat;
      caxis([-0.3 0.3]);
      hold on;
      contour(HH,[0 0],'k');
      colorbar
      contour(E,[0.0:0.05:0.5],'r');
      keyboard
    end
    
      
    
    if Emn<0,
      fprintf('!!!!    BG mean(E)<0 %7.4d\n',Emn);
%      keyboard;
    end
    
%    fprintf('Lapl.E=%6.3d\n\n',Le);
    fprintf('BG max(E)=%6.3f, mean(E)=%6.3f\n', Emx,   Emn);
    fprintf('EU min(E)=%6.3f, mean(E)=%6.3f\n\n', Ezmin, Ezmn);

    IN = ARC.IN(1:4:end); % subsample
    LAPLE.Iocn               = IN;
    LAPLE.TM(cc,1)           = dnmb;
%    LAPLE.d2E_BG(cc,1) = Le;
    LAPLE.Emean_BG(cc,1)     = Emn;
    LAPLE.Emax_BG(cc,1)      = Emx;
    LAPLE.Emean_Euras(cc,1)  = Ezmn;
    LAPLE.Emin_Euras(cc,1)   = Ezmin;
    LAPLE.SSH(cc,:)          = E(IN); 
   
  end
  
end

if s_mat>0
  dmm = LAPLE.SSH;
  LAPLE.SSH = single(dmm);
  fprintf('Saving %s\n',fmat);
  save(fmat,'LAPLE','-v7.3');
end

