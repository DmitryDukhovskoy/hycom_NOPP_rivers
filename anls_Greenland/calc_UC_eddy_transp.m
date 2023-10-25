% Calculate edy transport of tracer concentration 
% Following Nagai et al., 2015, Dominant role of eddies ..., JGR
% The methodology is:
% 1) Calculate 1-month mean fields monthly_mean_flds.m
% 2) Cacl. anomalies: subtract means from instanteneous (daily mean flds)
% 3) Average fluxes for 1 mo
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat=0;
Nlev=10;  % in the upper Nlev
%Zmn = 100; % average over the top Zmn m
rg = 9806;
hgg= 1e20; % 
Nav=30;   % averaging time, days

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_monthly_mean/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);
YRPLT=[];
cc=0;

%if expt=='061', yend=2008; end;
for icyc=2:2
  for iyr=2007:2007
    for imo=1:12
      cc=cc+1;
      YRPLT(cc,1)=iyr;
      YRPLT(cc,2)=imo;
      YRPLT(cc,3)=icyc;
    end
  end
end
% Start from:
ip1=1;

%ftopo = sprintf('%sdepth_ARCc0.08_09.nc',pthtopo); % 
ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % in new GOFS3. expts
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
  fmean = sprintf('%smonth_UVCSdP_%4.4i%2.2i-%i.mat',...
		    pthmat,yr,imo,icyc);
  fprintf('Loading mean: %s\n',fmean);
  load(fmean);
  Umn = squeeze(FMEAN.Uav(1:Nlev,:,:));
  Vmn = squeeze(FMEAN.Vav(1:Nlev,:,:));
  TRmn= squeeze(FMEAN.Trac_av(1:Nlev,:,:));
  dHmn= squeeze(FMEAN.dHav(1:Nlev,:,:));
  Smn = squeeze(FMEAN.Sav(1:Nlev,:,:));
  clear FMEAN;
  TRmn(TRmn<0)=0;
  

  cnc=0;
  fout = sprintf('%smonth_UC_eddy_%4.4i%2.2i-%i.mat',pthmat,yr,imo,icyc);
  
  for imday=1:mday
    dnmb=datenum(yr,imo,imday);
    iday=dnmb-dJ1+1;

    fina = sprintf('%s%s_archm.%4.4i_%3.3i_12nAtl.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%s_archm.%4.4i_%3.3i_12nAtl.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end
    
    cnc=cnc+1;
%    dnmb=datenum(yr,1,1)+iday-1;
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
    SS=squeeze(F(1:Nlev,:,:));

    [F,n,m,l] = read_hycom(fina,finb,'thknss');
    F=F./rg;
    F(F>hgg)=0;
    dH=squeeze(F(1:Nlev,:,:)); 

    UN(dH<0.001)=nan;
    VN(dH<0.001)=nan;
    TR(dH<0.001)=nan;
    SS(dH<0.001)=nan;

% Anomalies:
% To get tracer flux U*C in kg/s, where U=u*dH (m2/s)
% need to multiply Upr*TRpr (kg/[s*m]) by dx or dy - grid size length
% 
    Upr=UN.*dH-Umn.*dHmn; % m2/s
    Vpr=VN.*dH-Vmn.*dHmn;
    TRpr=TR-TRmn;   % kg/m3
    Spr=SS-Smn;

%      UpCp0=squeeze(Upr(1,:,:).*TRpr(1,:,:)); % surf. layer
%      VpCp0=squeeze(Vpr(1,:,:).*TRpr(1,:,:)); % surf. layer
    UpCpZ=squeeze(sum(Upr.*TRpr,1)); % intgr over N layer
    VpCpZ=squeeze(sum(Vpr.*TRpr,1)); % intgr over N layer
    UpSpZ=squeeze(sum(Upr.*Spr,1)); % intgr over N layer
    VpSpZ=squeeze(sum(Vpr.*Spr,1)); % intgr over N layer

    UCPR.TM(cnc,1)           = dnmb;
    UCPR.Nlev                = Nlev;
    UCPR.UpTRp_Nlev(cnc,:,:) = UpCpZ;
    UCPR.VpTRp_Nlev(cnc,:,:) = VpCpZ;
    UCPR.UpSp_Nlev(cnc,:,:)  = UpSpZ;
    UCPR.VpSp_Nlev(cnc,:,:)  = VpSpZ;

    fprintf('=======   1 day %8.3f sec\n\n',toc);


    fprintf('Monthly means: Uav=%8.2f %8.2f\n',...
	    min(min(Umn(1,:,:))),max(max(Umn(1,:,:))));
    fprintf('Monthly means: Vav=%8.2f %8.2f\n',...
	    min(min(Vmn(1,:,:))),max(max(Vmn(1,:,:))));
    fprintf('Monthly means: Trav=%8.2f %8.2f\n',...
	    min(min(TRmn(1,:,:))),max(max(TRmn(1,:,:))));
    fprintf('Monthly means: Sav=%8.2f %8.2f\n',...
	    min(min(Smn(1,:,:))),max(max(Smn(1,:,:))));

    if mod(cnc,7)==0 
      fprintf('####  Saving %s\n\n',fout);
      save(fout,'UCPR');
    end
  end  % days  
%keyboard
  fprintf('####  Saving %s\n\n',fout);
  save(fout,'UCPR');

end





