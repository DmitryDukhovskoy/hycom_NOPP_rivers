% Extract T&S at all layers HYCOM grid
% at specified point locations
% saved by years;
% 
% expt_110 - no Greenland runoff
% expt 112 - with Greenland runoff
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=2005;
YR2=2006;
regn = 'ARCc0.08';
%expt = 110; % no Greenland runoff  
expt = 112;  % Greenland runoff

fprintf('expt %3.3i\n\n',expt);


s_mat = 1; % =2 - load and start from last saved
%pfld  = 'salin';

s_fig = 0;
rg = 9806;
fprintf('Extracting T&S, %i-%i\n',YR1,YR2);


pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

btx='extr_TS_tser_vertlrs_pnt.m';

ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);


% Points, arc08 indices:
IJp=[490         596
         431         449
         509         401
         657         176
         712         470
         944         418
         967         611
        1017         790
        1079         546
        1151         690];
np=length(IJp);

% Plot
f_chck=0;
if f_chck==1
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-500 -500],'b');
  contour(HH,[-10000:1000:-900],'c');
  axis('equal');
  
  for ik=1:np
    i0=IJp(ik,1);
    j0=IJp(ik,2);
    plot(i0,j0,'r*');
    text(i0+5,j0,sprintf('%2.2i',ik),'Fontsize',12);
  end
  title('T/S time series locations');
  
  bottom_text(btx,'pwd',1);
end

%  keyboard
for iyr=YR1:YR2
  yr=iyr;
  dd=5;

  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
  if expt~=110,
    pthbin = sprintf('/nexsan/hycom/%s_%3.3i/data/%i/',regn,expt,yr);
  end

  cc=0;
  TS=struct;
  TS.IJ_points=IJp;
  for iday=1:dd:365
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end
    if ~exist(finb,'file');
      fprintf('Not found: %s\n\n',finb);
      continue;
    end

    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    mo=DV(2);
    mday=DV(3);

    cc=cc+1; % day counter
    fprintf('%i: %4.4i_%2.2i_%2.2i: %s\n',cc,DV(1:3),fina);

    tic;
    
   % Layer thickness:
    [F,nn,mm,ll] = read_hycom(fina,finb,'thknss');
    F=squeeze(F);
    F=F./rg;
    F(F>1e10)=0;
    F(F<1e-2)=0;
    Lthck = F;
    Iz=find(Lthck<1e-4);

  % Create Depth array of interface depths:
  % Note these are BOTTOM interfaces 
  % So Layer 1 is between interfaces 0m and ZZ(1)
    ZZb = F.*nan;
    ZZb(1,:,:) = 0;
    I = find(HH>=0);
    ZZb(1,I) = nan;
    for kk=1:ll
      ZZb(kk+1,:,:)=ZZb(kk,:,:)-Lthck(kk,:,:);
    end

  % Read T/S:    
    [F,n,m,l] = read_hycom(fina,finb,'salin');
    F(F>1e6)=nan;
    S = squeeze(F);
    S(Iz)=nan;

    [F,n,m,l] = read_hycom(fina,finb,'temp');
    F(F>1e6)=nan;
    T = squeeze(F);
    T(Iz)=nan;

    for ip=1:np
      i0=IJp(ip,1);
      j0=IJp(ip,2);
      TS.S(ip,:,cc)=squeeze(S(:,j0,i0));
      TS.T(ip,:,cc)=squeeze(T(:,j0,i0));
      TS.Z(ip,:,cc)=squeeze(ZZb(:,j0,i0));
    end
    TS.TM(cc)=dnmb;
    
    fprintf('Reading 1 rec: %8.6f min\n\n',toc/60);


  fmat = sprintf('%sarc08_%3.3i_TS_timeser_%i.mat',pthmat,expt,yr);
  if s_mat>0 & mod(cc,10)==0
    fprintf('Saving %s\n',fmat);
    save(fmat,'TS');
  end

  end;   % days
  
  if s_mat>0
    fprintf('Saving %s\n',fmat);
    save(fmat,'TS');
  end
  
end;  % year loop



  



