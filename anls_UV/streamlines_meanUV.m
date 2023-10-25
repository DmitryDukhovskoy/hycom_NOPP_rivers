% Prepare streamlines for
% U,V field (annual or monthly or any mean field)
% Plotting monthly mean UV fields
% derived in mnthly_arc08_UV.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;
f_mat = 1;
s_par = 1;

plr =1;  % U from plr layer
rg = 9806;

if s_par>0	 
  delete(gcp('nocreate'))
  if exist('spool','var'),
    delete(spool);
  end

  spool = parpool('local');
  
end

regn = 'ARCc0.08';
expt = 110;
pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

%YRPLT = [2011,2012,2013,2014,2015];
% Years to average or get monthly field:
YRS = 2005; % more than 1 year - mean over these years
np = length(YRS);
imo   = 13; % >12 - annual mean
if np>1, imo=13; end;

if imo>13
  fprintf('Streamlines for annual mean %i, mat saved 0/1: %i\n',YRS,f_mat);
end
if imo<13
  fprintf('Streamlines for month %i/%i, mat saved 0/1: %i\n',imo,YRS,f_mat);
end


ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
hmsk=HH;
hmsk(HH<0)=nan;

[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);
xlim1 = 20;
xlim2 = nn-1;
ylim1 = 100;
ylim2 = mm-100;

usm=zeros(mm,nn);
vsm=zeros(mm,nn);
cc = 0;
for ik=1:np
  iyr = YRS(ik);
  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  if imo>12
    for ik=1:12
      cc = cc+1;
      U = meanUV(ik).U;
      V = meanUV(ik).V;
      usm = usm+U;
      vsm = vsm+V;
    end
  else
    cc=1;
    usm = meanUV(ik).U;
    vsm = meanUV(ik).V;
  end
end;

U = usm./cc;
V = vsm./cc;

S = sqrt(U.^2+V.^2);
%keyboard

v_chck=0;
if v_chck>0
  figure(2); clf;
  rcmp=colormap_red(100);
  pcolor(S); shading flat;
  hold on;
  caxis([0 0.3]);
  colormap(rcmp);
  colorbar
  
  xl1=800;
  xl2=1250;
  yl1=500;
  yl2=1100;
  axis('equal');
  set(gca,'xlim',[xl1 xl2],...
	  'ylim',[yl1 yl2]);
  dd=5;
  scl=6;
  cf=0.3;
  beta=15;
  lwd=1.2;
  v_col=[0 0 0];
  dii=5;
  for ii=xl1:dii:xl2
    for jj=yl1:dii:yl2
      clear u v
      u = U(jj,ii);
      v = V(jj,ii);
      s = S(jj,ii);
      if isnan(u), continue; end;
  %    if res>0,
	u=u/s;
	v=v/s;
	scl=25;
  %    end

      x0=ii;
      y0=jj;

      x1=x0+u*scl;
      y1=y0+v*scl;
      draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
    end
  end

end


% =============================
%
% Lagrangian tracking 
%
% =============================
if ~exist('IG','var');
  [IG,JG]=meshgrid([1:nn],[1:mm]);
end

% Define time for advecting particles
% in steady-state ocean
YRPLT = [];
YR1=YRS(1);
cc=0;
dT=2;
for nyr=1:1
  for idd=1:dT:90
    cc=cc+1;
    yr = YR1+nyr-1;
    YRPLT(cc,1)=yr;       % forcing year, U,V fields
    YRPLT(cc,2)=idd;
    YRPLT(cc,3)=yr;       % same as YRPLT(1,:)
  end
end


% 
% Initialize particles:
% Set random generator 
rng('default');

PARAM.EXPERIMENT         = 'HYCOM ARCc0.08 expt_11.0';
PARAM.Outp_Freq_dt       = dT*24*3600; % model fields - output freq., sec
PARAM.Seed_Nmb_Prtcl     = 1;
PARAM.Max_Nmb_Prtcl      = 20000;
%PARAM.BG_Indx            = iBG;
PARAM.Diffusivity_RWalk  = 0;    % m2/s, eddy coef., diff, rand walk, HYCOM 1/12
dt = PARAM.Outp_Freq_dt;
%nu = PARAM.Diffusivity_RWalk;
nu = 0;

PRTCL=struct;
PRTCL.PARAM=PARAM;
PRTCL.TRACK(1).dnmb =[];
PRTCL.TRACK(1).Xcart=[];
PRTCL.TRACK(1).Ycart=[];
PRTCL.TRACK(1).I    =[];
PRTCL.TRACK(1).J    =[];

% ----------------
% Initialize:
% ----------------
if f_mat<2
  icnt=1;
  yr=YRPLT(1,1);
  iday=YRPLT(1,2);

  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);
  
  fprintf('Initializing: %4.4i_%2.2i_%2.2i \n',DV(1:3));

  PRTCL.TRACK(icnt).TM=dnmb;
%
% Seed Particles
  PRTCL = sub_seed_Greenl(HH,LON,LAT);
  UN = U;
  VN = V;
  
  f_pltg=0;
  if f_pltg>0
    figure(1); clf;
    contour(HH,[0 0],'k');
    hold on;
    contour(HH,[-200 -200],'c');
    contour(HH,[-1000 -1000],'g');

    contour(LON,[-180:20:180],'Color',[0.6 0.6 0.6]);
    contour(LON,[19 19],'Color',[0.6 0.6 0]);
    contour(LON,[5 5],'Color',[1 0.3 0]);
    contour(LAT,[50:10:88],'Color',[0.6 0.6 0.6]);
    contour(LAT,[78.5 78.5],'Color',[1 0.3 0]);

  %  II=SEED.XBar;
  %  JJ=SEED.YBar;
    plot(II,JJ,'b.-');
    II=PRTCL(1).TRACK.I;
    JJ=PRTCL(1).TRACK.J;
    plot(II,JJ,'b.-');
  end;

end;  % initialization


Npp=size(YRPLT,1);
ip1=2;
yrF = YRPLT(1,3);
fmat = sprintf('%sNAtlGreenl_particles_%i.mat',pthmat,yrF);
if f_mat>1
  fprintf('\n   Loading saved %s\n\n',fmat);
  load(fmat);
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yrF)
  [ip1,icnt,UN,VN] = sub_start_mat(PRTCL,YRPLT,expt,pthbin); % <-- need to edit
%keyboard
end

% ==================
% START
% ==================
yr_old   = YRPLT(ip1,3);
for ip=ip1:Npp
  yr   = YRPLT(ip,3);
  iday = YRPLT(ip,2);
  yrF = yr;
  
  %  Date:
  dnmb=datenum(yrF,1,1)+iday-1; % 
  DV=datevec(dnmb);
  
  fprintf('Calculating %4.4i-%2.2i-%2.2i\n',DV(1:3));
  
%  if mod(ip,10)==0
%    fmat = sprintf('%sBG_particles_%i.mat',pthmat,yr_old);
%    fprintf('=========   Saving %s\n',fmat);
%    save(fmat,'PRTCL');
%
%    clear TRACK
%    TRACK = PRTCL.TRACK(icnt); % save last record in newt array
%    icnt=1;
%    clear PRTCL
%    PRTCL = struct;
%    PRTCL.TRACK(1)=TRACK;
%    yr_old = yr;
%  end
  
  icnt=icnt+1;
  U(isnan(U))=0;
  V(isnan(V))=0;
  UN(isnan(UN))=0;
  VN(isnan(VN))=0;
  UO=UN;
  VO=VN;
  UN=U;
  VN=V;

% Initialize main variables for parallel loop:   
  Ip=PRTCL.TRACK(icnt-1).I;
  Jp=PRTCL.TRACK(icnt-1).J;
  np=find(~isnan(Ip));
  fprintf('\nActive particles: %i\n',length(np));
  
  clear II JJ TT SS ZZ
  lnp = length(Ip);
%  lnp   = length(np);
  II    = zeros(lnp,1)*nan;
  JJ    = zeros(lnp,1)*nan;
  dd    = ip;

  tic;  
  parfor ii=1:lnp
%  for ii=1:lnp
    if mod(ii,1000)==0
      fprintf('ii=%i - day %i\n',ii,dd);
    end
  
    LX = Ip(ii);
    LY = Jp(ii);
      
    [AX,AY]=runge_kutta(UO,VO,dt,...
                  UN,VN, LX, LY, DX, DY, HH, nu);
% Update particle locations:
%      LX = AX;
%      LY = AY;
%    end
%    toc

    if AX<0
      AX=nan;
      AY=nan;
    end

    II(ii)=AX;
    JJ(ii)=AY;

  end;
  fprintf('1 iteration time: %7.2f min\n',toc/60);  

  PRTCL.TRACK(icnt).TM = dnmb;
  PRTCL.TRACK(icnt).I  = II;
  PRTCL.TRACK(icnt).J  = JJ;
  
% keyboard 
  if mod(icnt,5)==0 & f_mat>0
%    fmat = sprintf('%sBG_particles_%i.mat',pthmat,yrF);
    fprintf('\n###  Saving %s\n',fmat);
    save(fmat,'PRTCL');
  end
  
end;  % end time loop


if f_mat>0
%  fmat = sprintf('%BG_particles_%i.mat',pthmat,yr);
  fprintf('\n###  Saving %s\n',fmat);
  save(fmat,'PRTCL');
end

if exist('spool','var'),
  delete(spool);
end





%keyboard

  
  




