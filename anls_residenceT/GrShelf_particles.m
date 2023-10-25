% Parallel code
% particles tracking 
% seeded at on southwestern Greenland Shelf
% and western Labrador shelf
%
% Particles are not added during the simulation
% All N prticles seeded at once at initial state
% For each particle: 
% no T or S is tracked, only time and location
% 
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

%shelf='GrSh';
shelf='WLabr'; 
% Layers:
% lr=10 ~ 45-52 m, 
% lr=15    84-92 m
% lr=23   148-156 m
% lr=31  ~400-500 m 
slr = 15;   % layer where to seed particles
Np0 = 100;   % number of initial particles
nsim = 4;   % multiple Lagr. simulations

s_par=1; % parallel session - run serial 
%f_plt=1; % plot prtcles
f_mat=2018; % =1 - start from the 1st day, save 
            % =YEAR load saved and start from the last record = YEAR, 
	    % save at the end, e.g. f_mat=1995
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
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
if strncmp(shelf,'WLabr',5)
  pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_WLabr_prt/';
end

% Repeat N years year 1994:
YRPLT=[];
cc=0;
dT=2; % model fields, outp. freq., days
YR1 = 1993; % create time array
for nyr=1:27
  iiy=0;
  for idd=1:dT:365
    cc=cc+1;
    iiy=iiy+1;
    yr = YR1+nyr-1;
    YRPLT(cc,1)=yr;       % forcing year, U,V fields
    YRPLT(cc,2)=idd;
    YRPLT(cc,3)=yr;       % same as YRPLT(1,:)
    
% Check existing fmat files to avoid overriding 
    fmat = sprintf('%s%s_part-lr%2.2i_%i-%2.2i.mat',pthmat,shelf,slr,yr,nsim);
    if exist(fmat,'file') & iiy==1 & yr>=f_mat
      fprintf('!!!  Mat file exists: %s\n',fmat);
      fprintf('!!! will be overridden, change nsim=%i ???\n',nsim);
      fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    end
    
  end
end

if f_mat == 1
  fprintf(' Saving output files from year YR1\n');
else
  fprintf('  Start saving output files from %i\n\n', f_mat);
end
fmat = sprintf('%s%s_part-lr%2.2i_XXXX-%2.2i.mat',pthmat,shelf,slr,nsim);
fprintf('The output will be saved to: %s\n',fmat);

% ------------------------
% TOPO
% ------------------------
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

% Define BG:


%% Convert geogr. -> cartesian, m
%fprintf('Converting to Cart. Grid ...\n');
%[XM,YM]=sub_dgr2m(LON,LAT);
% -----------------------------
% Calculate DX, DY of the grid
% -----------------------------
fprintf('Calculating DX, DY ...\n');
for i=1:nn-1
  lat1=LAT(:,i);
  lat3=LAT(:,i+1);

  lon1=LON(:,i);
  lon3=LON(:,i+1);

  dist = distance_spheric_coord(lat1,lon1,lat3,lon3);  % m
  DX(:,i)=dist;
end


for j=1:mm-1
  lat1=LAT(j,:);
  lat2=LAT(j+1,:);

  lon1=LON(j,:);
  lon2=LON(j+1,:);

  dist = distance_spheric_coord(lat1,lon1,lat2,lon2);  % m
  DY(j,:)=dist;
end;

DX(mm,:)=DX(mm-1,:);
DX(:,nn)=DX(:,nn-1);
DY(mm,:)=DY(mm-1,:);
DY(:,nn)=DY(:,nn-1);
if ~exist('IG','var');
  [IG,JG]=meshgrid([1:nn],[1:mm]);
end

% 
% Initialize particles:
% Set random generator 
rng('default');

PARAM.EXPERIMENT         = 'HYCOM ARCc0.08 expt_11.0';
PARAM.Outp_Freq_dt       = dT*24*3600; % model fields - output freq., sec
PARAM.Seed_Nmb_Prtcl     = 1;
PARAM.Max_Nmb_Prtcl      = 1000;
%PARAM.BG_Indx            = iBG;
PARAM.Diffusivity_RWalk  = 80;    % m2/s, eddy coef., diff, rand walk, HYCOM 1/12
dt = PARAM.Outp_Freq_dt;
%nu = PARAM.Diffusivity_RWalk;
nu = 40;

PRTCL=struct;
PRTCL.PARAM=PARAM;
PRTCL.TRACK(1).dnmb =[];
PRTCL.TRACK(1).Xcart=[];
PRTCL.TRACK(1).Ycart=[];
PRTCL.TRACK(1).I    =[];
PRTCL.TRACK(1).J    =[];
PRTCL.TRACK(1).T    =[];
PRTCL.TRACK(1).S    =[];
PRTCL.TRACK(1).Z    =[];

% ----------------
% Initialize:
% ----------------
if f_mat<2
  icnt=1;
  yr=YRPLT(1,1);
  iday=YRPLT(1,2);

  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
  if yr>2016
    expt = 112;
    pthbin = sprintf('/nexsan/hycom/ARCc0.08_112/data/%i/',yr);
  end


  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);
  
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
  end

  fprintf('Initializing: %4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

  PRTCL.TRACK(icnt).TM=dnmb;
  switch (shelf)
   case('GrSh');
    [PRTCL,UN,VN] = sub_seed_GrSh_prt(slr,HH,LON,LAT,fina,finb,Np0,nsim);
   case('WLabr');
    [PRTCL,UN,VN] = sub_seed_WLabr_prt(slr,HH,LON,LAT,fina,finb,Np0,nsim);
  end

  f_pltg=0;
  if f_pltg>0
    figure(1); clf;
    contour(HH,[0 0],'k');
    hold on;
    contour(HH,[-500 -500],'c');
    contour(HH,[-8000:1000:-1000],'Color',[0.5 0.5 0.5]);

    axis('equal');
    set(gca,'xlim',[350 700],...
	    'ylim',[300 800]);
    

  %  II=SEED.XBar;
  %  JJ=SEED.YBar;
    II=PRTCL(1).TRACK.I;
    JJ=PRTCL(1).TRACK.J;
    plot(II,JJ,'b.');
  end;

end;  % initialization
clear T S dH yr

Npp=size(YRPLT,1);
ip1=1;
yrF = YRPLT(1,3);
yr_old   = YRPLT(ip1,3);

if f_mat>1
  yrF = f_mat;
  fmat = sprintf('%s%s_part-lr%2.2i_%i-%2.2i.mat',pthmat,shelf,slr,yrF,nsim);
  fprintf('\n   Loading saved %s\n\n',fmat);
  load(fmat);
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yrF)
  if yrF>2016
    expt = 112;
    pthbin = sprintf('/nexsan/hycom/ARCc0.08_112/data/%i/',yrF);
  end
  [ip1,icnt,UN,VN] = sub_start_mat_layer(slr,PRTCL,YRPLT,expt,pthbin);
  yr_old=yrF;
%keyboard
end

% ==================
% START
% ==================
toc_av=0;
ntc=0;
for ip=ip1:Npp
  yr   = YRPLT(ip,3);
  iday = YRPLT(ip,2);
  yrF = yr;
  
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
  if yr>2016
    expt = 112;
    pthbin = sprintf('/nexsan/hycom/ARCc0.08_112/data/%i/',yr);
  end
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
  end

  %  Date:
  dnmb=datenum(yrF,1,1)+iday-1; % 
  DV=datevec(dnmb);
  
  fprintf('Layer %i, Sim %i, %s, Calculating %4.4i-%2.2i-%2.2i\n',...
           slr,nsim, shelf, DV(1:3));
%keyboard  
% Save if end of year: 
  if (yr_old ~= yr) 
   if ip>ip1 % do not need to save if start from mat file
    tm1=PRTCL.TRACK(end).TM;
    dvE=datevec(tm1);
    iyE=dvE(1);
    fmat = sprintf('%s%s_part-lr%2.2i_%i-%2.2i.mat',pthmat,shelf,slr,iyE,nsim);
    fprintf('========= SAVING: Last record in PRTCL=%s\n',datestr(tm1));
    fprintf('========= END OF YEAR  Saving %s\n',fmat);
    save(fmat,'PRTCL');
   end
   
    clear TRACK
    TRACK = PRTCL.TRACK(icnt); % save last record in newt array
    icnt=1;
    clear PRTCL
    PRTCL = struct;
    PRTCL.TRACK(1)=TRACK;
    yr_old = yr;
  end
  
  icnt=icnt+1;
  UO=UN;
  VO=VN;
  [F,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',slr);
  F(F>1e6)=nan;
  UN=squeeze(F);
  [F,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',slr);
  F(F>1e6)=nan;
  VN=squeeze(F);

% Initialize the main variables for parallel loop:   
  Ip=PRTCL.TRACK(icnt-1).I;
  Jp=PRTCL.TRACK(icnt-1).J;
  np=find(~isnan(Ip));
  fprintf('\nActive particles: %i\n',length(np));
  
% Time intervals for Runge-Kutta:  
%  tspan=[0:dt/dT:dt]; % 1 day
%  ntstp = length(tspan)-1;
  
  clear II JJ TT SS ZZ
  lnp = length(Ip);
%  lnp   = length(np);
  II    = zeros(lnp,1)*nan;
  JJ    = zeros(lnp,1)*nan;
  dd    = ip;
%keyboard
  tic;  
  parfor ii=1:lnp
%  for ii=1:lnp
    if s_par==0 & mod(ii,20)==0
      fprintf('ii=%i - day %i\n',ii,dd);
    end
  
    LX = Ip(ii);
    LY = Jp(ii);
%    tic
%    LX = Ip;
%    LY = Jp;
%    for itm = 1:ntstp
      
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

  amm=toc/60;
  toc_av=toc_av+amm;
  ntc=ntc+1;
  nitr=Npp-ip;
  htime=nitr*(toc_av/ntc)/60;

  fprintf('1 iteration time: %7.2f min\n',amm);  
  fprintf('Average time 1 iter: %7.2f min, Estimated time to end: %6.2f hrs\n\n',...
	  toc_av/ntc,htime);

  PRTCL.TRACK(icnt).TM = dnmb;
  PRTCL.TRACK(icnt).I  = II;
  PRTCL.TRACK(icnt).J  = JJ;
  
% keyboard 
  if mod(icnt,20)==0 & f_mat>0
    tm1=PRTCL.TRACK(end).TM;
    dvE=datevec(tm1);
    iyE=dvE(1);
    fmat = sprintf('%s%s_part-lr%2.2i_%i-%2.2i.mat',pthmat,shelf,slr,iyE,nsim);
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

