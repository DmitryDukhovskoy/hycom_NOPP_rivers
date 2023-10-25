% Calculate monthly mean fields from daily mean
% U, V, S, tracer conc, dP
% For the upper Nlev 
% averaging time period = Nav days
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat=1;
Nlev=20;  % save for Nlev
rg = 9806;
hgg= 1e20; % 
Nav=30;   % averaging time, days

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_monthly_mean/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);
tmpmat = sprintf('%sMonth_tmp.mat',pthmat);

% Need to modify code to include existing days
%if ~exist(flxmat,'file');

%figure out dates that needed for averaging
%yr=YRPLT(1,1);
%iday=YRPLT(1,2);
%icyc=YRPLT(1,3);
%dnmb=datenum(yr,1,1)+iday-1;
%DV=datevec(dnmb);

% Averaging time period:
% When # of records = Nav, the mean is 
% is saved at time d2-Nav/2, d2- end of averaging
YRPLT=[];
cc=0;
%if expt=='061', yend=2008; end;
for icyc=2:2
  for iyr=2005:2005
    for imo=1:12
      cc=cc+1;
      YRPLT(cc,1)=iyr;
      YRPLT(cc,2)=imo;
      YRPLT(cc,3)=icyc;
    end
  end
end
% Can restart from any month
% Start from:
ip1=1; 
fprintf('Starting from Month = %i\n',ip1);

%ftopo = sprintf('%sdepth_ARCc0.08_09.nc',pthtopo); % 
ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
% Collocate sub-regions, find index offset
INDa = smaller_domain_indices('NorthAtl');
inc1=INDa.i1;
inc2=INDa.i2;
jnc1=INDa.j1;
jnc2=INDa.j2;
djnc=INDa.dj;
dinc=INDa.di;
%keyboard

HH=HH(jnc1:jnc2,inc1:inc2);
LON=LON(jnc1:jnc2,inc1:inc2);
LAT=LAT(jnc1:jnc2,inc1:inc2);
[mm,nn]=size(LON);


nrec=0;  % # of saved mean U,V averaged over Nav days
np=length(YRPLT);

for ip=ip1:np
  yr=YRPLT(ip,1);
  imo=YRPLT(ip,2);
  icyc=YRPLT(ip,3);
  if icyc==1,
    expt='060';
  else
    expt='061';
  end
  pthbin = sprintf(...
  '/Net/yearly/ddmitry/hycom/ARCc0.08/%s/%4.4i/bin_outp/',expt,yr);

  d1=datenum(yr,imo,1);
  d2=d1+35;
  dv=datevec(d2);
  d2=datenum(dv(1),dv(2),1);
  mday=d2-d1;
  dJ1=datenum(yr,1,1);
  cnc=0;
  fmatout = sprintf('%smonth_UVCSdP_%4.4i%2.2i-%i.mat',...
		    pthmat,yr,imo,icyc);
  for imday=1:mday
    dnmb=datenum(yr,imo,imday);
    iday=dnmb-dJ1+1;

    fina = sprintf('%s%s_archm.%4.4i_%3.3i_12nAtl.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%s_archm.%4.4i_%3.3i_12nAtl.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
    else
      cnc=cnc+1;
      dnmb=datenum(yr,1,1)+iday-1;
      DV=datevec(dnmb);
      TM(cnc,1)=dnmb;
      CYC(cnc,1)=icyc;

      fprintf('%s: %4.4i_%2.2i_%2.2i: %s\n',expt,DV(1:3),fina);

      tic;
      [F,n,m,l] = read_hycom(fina,finb,'u-vel.');
  %  F=F(:,jnc1:jnc2,inc1:inc2);
      F(F>hgg)=0;
      UN=squeeze(F(1:Nlev,:,:));

      [F,n,m,l] = read_hycom(fina,finb,'v-vel.');
      F(F>hgg)=0;
      VN=squeeze(F(1:Nlev,:,:));

      [F,n,m,l] = read_hycom(fina,finb,'tracer');
      F(F>hgg)=0;
      TR=squeeze(F(1:Nlev,:,:));

      [F,n,m,l] = read_hycom(fina,finb,'salin');
      F(F>hgg)=0;
      SN=squeeze(F(1:Nlev,:,:));

      [F,n,m,l] = read_hycom(fina,finb,'thknss');
      F=F./rg;
      F(F>hgg)=0;
      dH=squeeze(F(1:Nlev,:,:)); 

      UN(dH<0.001)=nan;
      VN(dH<0.001)=nan;
      TR(dH<0.001)=nan;
      SN(dH<0.001)=nan;

      if cnc==1
	Uav  = zeros(Nlev,m,n);
	Vav  = zeros(Nlev,m,n);
	TRav = zeros(Nlev,m,n);
	Sav  = zeros(Nlev,m,n);
	dHav = zeros(Nlev,m,n);
      end
      Uav = Uav + UN.*dH;
      Vav = Vav + VN.*dH;
      TRav= TRav+ TR.*dH;
      Sav = Sav + SN.*dH;
      dHav= dHav+ dH;

      fprintf('=======   1 day %8.3f sec\n\n',toc);
      
    end  % if file exists
    
% Update:
    if imday==mday
      Uav = Uav./dHav;
      Vav = Vav./dHav;
      TRav= TRav./dHav;
      Sav = Sav./dHav;
      dHav= dHav/cnc;

      FMEAN.code     = 'monthly_mean_flds.m';
      FMEAN.TM_avrg  = TM;
      FMEAN.CYC      = CYC;
      FMEAN.Nlev     = Nlev;
      FMEAN.Nav_days = Nav;
      FMEAN.Uav      = Uav;
      FMEAN.Vav      = Vav;
      FMEAN.Trac_av  = TRav;
      FMEAN.Sav      = Sav;
      FMEAN.dHav     = dHav;

      fprintf('Monthly means: Uav=%8.2f %8.2f\n',...
	      min(min(Uav(1,:,:))),max(max(Uav(1,:,:))));
      fprintf('Monthly means: Vav=%8.2f %8.2f\n',...
	      min(min(Vav(1,:,:))),max(max(Vav(1,:,:))));
      fprintf('Monthly means: Trav=%8.2f %8.2f\n',...
	      min(min(TRav(1,:,:))),max(max(TRav(1,:,:))));
      fprintf('Monthly means: Sav=%8.2f %8.2f\n',...
	      min(min(Sav(1,:,:))),max(max(Sav(1,:,:))));
      fprintf('Monthly means: dH=%8.2f %8.2f\n',...
	      min(min(dHav(1,:,:))),max(max(dHav(1,:,:))));

      fprintf('####  Saving %s\n\n',fmatout);
      save(fmatout,'FMEAN');

    end;
  
  end  % days  
end

exit
