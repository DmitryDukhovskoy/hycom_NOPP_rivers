% 0.04 HHYCOM-CICE
%
% Analyze eddy fluxes of U,T,S
% across Greenland contour
% fields are prepared in:
% mean <u>, <t>, <s> in greenl_meanUTS_Zlvls.m
% perturbations: greenl_pertUTS_Zlvls.m: <u'*T'>
% 
% Also estimate eddy flux using indirect method:
% <u'*T'>=<u*T>-<u_bar>*<T_bar>

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat = 0; % =0 - load saved and plot; =2 - do not save

Tref=0; % if change Tref to not 0, need to recalculate <u*T> terms
Cp=4200; % J/kg*C
rhow=1027;  % water density, kg/m3

YR1=2005;
YR2=2007;
YRS = [YR1:YR2];
%YRS=2005;
nyr=length(YRS);
rg  = 9806;  % convert pressure to depth, m
hgg = 1e20; 

plr=0; % highlight this interface
btx = 'eddy_flux_arc04_Greenl.m';

fprintf('Calculate daily perturb U,T,S, Greenland Contour, on Zlevels, %i\n',YRS);

regn = 'ARCc0.04';
expt = 011; % experiment without runoff
%expt = 012;  % epxeriment with Greenland runoff and monthly Arctic rivers

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.04/data_mat2/';

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath_arc04(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section
np = length(Hs);


fmatout = sprintf('%sarc004_%3.3i_eddyTVolflx_%4.4i-%4.4i.mat',...
		  pthmat,expt,min(YRS),max(YRS));

if s_mat>0
  cc=0;
  kkf=0;
  for iyr=1:nyr
    yr=YRS(iyr);

  % Mean fields:  
    fmatu  = sprintf('%sarc04_expt%3.3i_greenl_contr_mnthUzlv_%i.mat',...
		     pthmat,expt,yr);
    fmatt  = sprintf('%sarc04_expt%3.3i_greenl_contr_mnthTzlv_%i.mat',...
		     pthmat,expt,yr);
    fmats  = sprintf('%sarc04_expt%3.3i_greenl_contr_mnthSzlv_%i.mat',...
		     pthmat,expt,yr);
    
 % Total mean flux: <u*T>   
    fmatut = sprintf('%sarc04_expt%3.3i_greenl_contr_meanUT_%i.mat',...
		     pthmat,expt,yr);

  % Perturbations
    fmatup = sprintf('%sarc04_expt%3.3i_greenl_contr_prtrbUzlv_%i.mat',...
		     pthmat,expt,yr);
    fmattp = sprintf('%sarc04_expt%3.3i_greenl_contr_prtrbTzlv_%i.mat',...
		     pthmat,expt,yr);
    fmatsp = sprintf('%sarc04_expt%3.3i_greenl_contr_prtrbSzlv_%i.mat',...
		     pthmat,expt,yr);

  % Heat fluxes: <u><t>, <u'*t'>
    fprintf('Loading flux files %s \n',fmatu);
    fprintf('Loading flux files %s \n',fmatup);
    fprintf('Loading flux files %s \n',fmatt);
    fprintf('Loading flux files %s \n',fmattp);
    fprintf('Loading flux files %s \n',fmatut);
    load(fmatu); 
    load(fmatup);
    load(fmatt);
    load(fmattp);
    load(fmatut); % <u*T>

    if ~exist('dL','var');
      dL=UZGR.Distance_m*1e-3; % km
      ZZ=UZGR.ZZlevels;
      dZ=abs(diff(ZZ));
      ddL=diff(dL(:));
      ddL=[ddL;ddL(end)];
    end

    Umn=UZGR.U;
%    Upr=UpZGR.U;
    Tmn=TZGR.T;
    TpUp=TpZGR.UpTp; % <u'*t'> - averaged by months
    UT=UTZGR.UmTm; % <u*T> - averaged by months

    for im=1:12
      um_bar=squeeze(Umn(im,:,:)); % u_bar - reference mean u
%      up=squeeze(Upr(im,:,:));
      tm_bar=squeeze(Tmn(im,:,:));  % t_bar
      
      if im==1
	avum=um_bar*0;
	avtm=tm_bar*0;
      end
      avum=avum+um_bar;
      avtm=avtm+tm_bar;
      
      tpup=squeeze(TpUp(im,:,:));
      uTtot=squeeze(UT(im,:,:));

      cc=cc+1;
      if cc==1
	avVFm=ddL*0;
	avHFp=ddL*0;
	avHFtot=ddL*0;
	uTvrt=dZ*0;  % vertical profile of Cp*rho*<uT>
	hfp_vrt=dZ*0;  % vert prof of Cp*rho*<u'*T'>
      end

  % Perturb: Cp*rho*<u'*T'>
      hfp=Cp*rhow*tpup;
      Iocn=find(~isnan(hfp(1,:)));
      hfp(isnan(hfp))=0;
      HFp=hfp'*dZ;
      intgr_hfp=hfp*ddL/sum(ddL(Iocn)); % along contour-mean flx, W/m2
      hfp_vrt=hfp_vrt+intgr_hfp; % W/m2

  % Mean Vol Flux m3/s per 1 m of segment length:
      vf=um_bar;
      vf(isnan(um_bar))=0;
      VFm=vf'*dZ;

      avHFp=avHFp+HFp;
      avVFm=avVFm+VFm;  % <---- not balanced, overall vol flux not 0
      
% average monthly mean u*T terms (total flux)
% for indirect calculation of <u'*T'>
% <u*(T-Tref)>=...=<uT>-<u>*Tref, the first
% term is saved in greenl_pertUTS_Zlvsl.m
      uTtot(isnan(uTtot))=0;
      hftot=Cp*rhow*(uTtot-Tref*um_bar);
      hftot(isnan(um_bar))=0;
      avHFtot=avHFtot+(hftot'*dZ); % Cp*rho*<u*T> intgr over depth
      
% Save vertical mean profile of <u*T>
      intgr_hftot=hftot*ddL/sum(ddL(Iocn)); % along contour-mean flux, W/m2
      uTvrt=uTvrt+intgr_hftot;  % W/m2
      
% Check cov(u,T), upper ocean
      uu=squeeze(Umn(im,1,:));
      tt=squeeze(Tmn(im,1,:));
      
      ums=nanmean(uu); % mean along contour
      tms=nanmean(tt);
      
      crrut(cc,1)=nanmean((uu-ums).*(tt-tms))/(nanstd(uu)*nanstd(tt));

% Plot mean U_bar
f_plt=0;
      if f_plt==1
	ZM=ZZ(2:end);
	pcolor(dL,ZM,um_bar); shading flat
	set(gca,'ylim',[-1500 0])
	title('<u_{bar}>');
        bottom_text(btx,'pwd',1);
      end
      
    end % months
% Monthly mean <u_bar>*<t_bar>
% Mean H Flux W per 1 m of segment length::
    aTm=avtm./im;
    aUm=avum./im;
    hf=Cp*rhow*(aTm-Tref).*aUm;
    hf(isnan(hf))=0;
    HFm=hf'*dZ; 

    kkf=kkf+1; % years
    if kkf==1
      avHFm=HFm*0; % time averaged heat flux - mean component
    end
    avHFm=avHFm+HFm;

    
  end

  % Save time-averaged fluxes:
  % Cp*rho*<u>*<t>:
  HVFLX.Title='Greenl Contour, Mean and perturb fluxes averaged over time';
  HVFLX.Tref=Tref;
  HVFLX.HeatFlx_mean_W=avHFm/kkf;     % Cp*rho*<T_bar-Tref>*<u_bar>
  HVFLX.HeatFlx_prtrb_W=avHFp/cc;     % Cp*rho*<u'*T'>
  HVFLX.VolFlx_mean_m3_sec=avVFm/cc;  % <u*u>
  HVFLX.uT_total=avHFtot/cc;          % Cp*rho*<u*T>
  HVFLX.Dist_contoru_km=dL;
  HVFLX.corr_U_T=crrut;
  HVFLX.uT_vertical=uTvrt;              % vertical distr of Cp*rho*<u*T>
  HVFLX.HeatFlx_prtrb_vertical=hfp_vrt; % vert of Cp*rho*<u'*T'>

  if s_mat==1
    fprintf('Saving %s\n',fmatout);
    save(fmatout,'HVFLX');
  end
  
else
  fprintf('Loading %s\n',fmatout);
  load(fmatout);
  
  fmatu  = sprintf('%sarc04_expt%3.3i_greenl_contr_mnthUzlv_%i.mat',...
		     pthmat,expt,YR1);
  load(fmatu); 
  
  if ~exist('ZZ','var');
%    dL=UZGR.Distance_m*1e-3; % km
    ZZ=UZGR.ZZlevels;
    dZ=abs(diff(ZZ));
  end
end

dL=HVFLX.Dist_contoru_km;
hfm=HVFLX.HeatFlx_mean_W;       % Cp*rho*<T_bar-Tref>*<u_bar> - mean heat flx
hfp=HVFLX.HeatFlx_prtrb_W;      % Cp*rho*<u'*T'>
vfm=HVFLX.VolFlx_mean_m3_sec;   % <u*u>
uTm=HVFLX.uT_total;             % Cp*rho*<u*T> - total heat flux

% ===========================
% Bug - use wrong starting contour pnt
% in Greenland contour isobath
ftmp=0;
if ftmp==1
 idx=1474; % index from the end 
  nll=length(uTm);
%  a1=dL(1:nll-idx);
%  a2=dL(nll-idx+1:end);
%  a1=a1(:);
%  a2=a2(:);
%  dL=[a2;a1];

  a1=hfm(1:nll-idx);
  a2=hfm(nll-idx+1:end);
  hfm=[a2;a1];

  a1=hfp(1:nll-idx);
  a2=hfp(nll-idx+1:end);
  hfp=[a2;a1];

  a1=vfm(1:nll-idx);
  a2=vfm(nll-idx+1:end);
  vfm=[a2;a1];

  a1=uTm(1:nll-idx);
  a2=uTm(nll-idx+1:end);
  uTm=[a2;a1];
end
% ===========================

ddL=diff(dL(:));
ddL=[ddL;ddL(end)];
I=find(ddL==0);
if ~isempty(I),
  error('Repeated nodes along the section ...');
end
%mdl=mean(ddL(I))*1e3; %mean segm dist, m

% Butterworth filter
% 1 mo cutoff = 1/12 -> in Matlab Wn=1/6
% 
nl=round(length(ddL)/2); 
iid=60; % # pnts for filtering
Wn = iid/nl;
avwdth=round(dL(iid)); % averaging window, inv of cutoff freq
[Bf,Af] = butter(9,Wn,'low');
[Bfh,Afh] = butter(9,1/20,'low');
yy = filtfilt(Bf,Af,hfm); % 
hfmf=yy;  % filtered mean heat flux, J/s per 1 m segm
yy = filtfilt(Bf,Af,vfm); % 
vfmf=yy;  % filtered vol flux, m3/s per 1 m of segment
%vfmf=yy;  % filtered vol flux, m3/s 
yy = filtfilt(Bf,Af,hfp); % 
hfpf=yy;  % filtered heat flux perturb, J/s per 1 m of segm

yy = filtfilt(Bf,Af,uTm);
uTmf=yy;



% Integrate heat flux along the contour:
% Note that directly estimated total heat flux
% is for Tref=0C
iuTm=sum(uTm.*ddL*1e3);
iuTm2=sum((hfm+hfp).*ddL*1e3);
fprintf('==== Directly estimated integrated total heat flux =%5.4f PW ==\n',iuTm*1e-15);
fprintf('==== Inirectly estimated integrated total heat flux=%5.4f PW ==\n',iuTm2*1e-15);



figure(1); clf;
axes('Position',[0.08 0.72 0.88 0.2]);
plot(dL,vfmf,'linewidth',1.6);
stl=sprintf('$$004-%3.3i~VolFlx, m3/(s*m),~%i-%i,~ButterwFltr~dw=%5.4f~km^{-1}$$',...
	    expt,min(YRS),max(YRS),(1/avwdth));
title(stl,'Fontsize',12,'interpreter','latex');
set(gca,'tickdir','out',...
	'xlim',[100 9050],...
	'xtick',[0:500:9500],...
	'xgrid','on',...
	'ygrid','on');
set(gcf,'position',[1277 374 1104 967]);
xlabel('Distance along Greenland contour, km');

f_corr=0;
if f_corr==1
  axes('Position',[0.08 0.42 0.88 0.2]);
  plot(HVFLX.corr_U_T);
  stl=sprintf('$$004-%3.3i~Corr(u,t)$$',expt);
  title(stl,'Fontsize',12,'interpreter','latex');
  set(gca,'tickdir','out',...
	  'xlim',[1 length(HVFLX.corr_U_T)],...
	  'xtick',[0:5:100],...
	  'xgrid','on',...
	  'ygrid','on');
  set(gcf,'position',[1277 374 1104 967]);
  xlabel('Time, mo');
end

% Plot mean total heat flux: Cp*rho*<u*T>,
% should be same as <u_bar>*<T_bar>+<u'*T'>
axes('Position',[0.08 0.42 0.88 0.2]);
plot(dL,uTmf,'linewidth',1.6);
stl=sprintf('$$004-%3.3i~C_p\\rho<u*t>,~W,~%i-%i,~BttrwFltr~dw=%5.4f~km^{-1}$$',...
	    expt,min(YRS),max(YRS),(1/avwdth));
title(stl,'Fontsize',12,'interpreter','latex');
set(gca,'tickdir','out',...
	'xlim',[100 9050],...
	'xtick',[0:500:9500],...
	'xgrid','on',...
	'ygrid','on');
set(gcf,'position',[1277 374 1104 967]);
xlabel('Distance along Greenland contour, km');

% Indirect estimate of Cp*rho*<u'*T'>=Cp*rho*<u*T>-Cp*rho*<u_bar>*<T_bar>
hfp2=uTm-hfm;
yy = filtfilt(Bf,Af,hfp); %
hfp2=yy;

axes('Position',[0.08 0.09 0.88 0.2]);
plot(dL,hfp2,'linewidth',1.6);
stl=sprintf('$$004-%3.3i~C_p\\rho<u''*t''>=C_p*\\rho(<uT>-<u><T>)$$',...
	    expt);
title(stl,'Fontsize',12,'interpreter','latex');
set(gca,'tickdir','out',...
	'xlim',[100 9050],...
	'xtick',[0:500:9500],...
	'xgrid','on',...
	'ygrid','on');

set(gcf,'position',[1277 374 1104 967]);
xlabel('Distance along Greenland contour, km');
bottom_text(btx,'pwd',1);
set(gcf,'position',[1277 374 1104 967]);





figure(2); clf;
axes('Position',[0.08 0.72 0.88 0.2]);
plot(dL,hfmf+hfpf,'linewidth',1.6);
stl=sprintf('$$004-%3.3i~<u><t>+<u''*t''>~W,%i-%i,~ButtrFltr~dw=%5.4f~km^{-1}$$',...
	    expt,min(YRS),max(YRS),(1/avwdth));
title(stl,'Fontsize',12,'interpreter','latex');
set(gca,'tickdir','out',...
	'xlim',[100 9050],...
	'xtick',[0:500:9500],...
	'xgrid','on',...
	'ygrid','on');

axes('Position',[0.08 0.42 0.88 0.2]);
plot(dL,hfmf,'linewidth',1.6);
stl=sprintf('$$004-%3.3i~C_p\\rho<u><t>~W,~ButtrFltr~dw=%5.4f~km^{-1}$$',...
	    expt,(1/avwdth));
title(stl,'Fontsize',12,'interpreter','latex');
set(gca,'tickdir','out',...
	'xlim',[100 9050],...
	'xtick',[0:500:9500],...
	'xgrid','on',...
	'ygrid','on');

axes('Position',[0.08 0.09 0.88 0.2]);
plot(dL,hfpf,'linewidth',1.6);
stl=sprintf('$$004-%3.3i~C_p\\rho<u''*t''>~W,~ButterwFltr~dw=%5.4f~km^{-1}$$',...
	    expt,(1/avwdth));
title(stl,'Fontsize',12,'interpreter','latex');
set(gca,'tickdir','out',...
	'xlim',[100 9050],...
	'xtick',[0:500:9500],...
	'xgrid','on',...
	'ygrid','on');

set(gcf,'position',[1277 374 1104 967]);
xlabel('Distance along Greenland contour, km');
bottom_text(btx,'pwd',1);
set(gcf,'position',[1277 374 1104 967]);


% Plot vertical distr of fluxes
% These are depth-line integrated, units W
nl=length(dZ);
ZM=ZZ(1:nl)-dZ/2;
hfp_vrt=HVFLX.HeatFlx_prtrb_vertical;
uTvrt=HVFLX.uT_vertical;

figure(3); clf;
axes('Position',[0.08 0.2 0.35 0.7]);
xl1=min(hfp_vrt);
xl2=1.05*max(hfp_vrt);
plot(hfp_vrt,ZM,'linewidth',1.6);
set(gca,'tickdir','out',...
	'xlim',[xl1 xl2],...
	'ylim',[-1000 0],...
	'ytick',[-1000:100:0],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',12);
stl=sprintf('$$004-%3.3i~C_p\\rho<u''*t''>,~W/m^2,~%i-%i$$',...
	    expt,YR1,YR2);
title(stl,'Fontsize',12,'interpreter','latex');


axes('Position',[0.58 0.2 0.35 0.7]);
xl1=min(uTvrt);
xl2=1.05*max(uTvrt);
plot(uTvrt,ZM,'linewidth',1.6);
hold
plot(hfp_vrt,ZM,'Color',[0.8 0.3 0]);

set(gca,'tickdir','out',...
	'xlim',[xl1 xl2],...
	'ylim',[-1000 0],...
	'ytick',[-1000:100:0],...
	'xgrid','on',...
	'ygrid','on');
stl=sprintf('$$004-%3.3i~C_p\\rho<u*t>,~W/m^2,~%i-%i$$',...
	    expt,YR1,YR2);
title(stl,'Fontsize',12,'interpreter','latex');
legend('<u*t>','<u''*t''>','Location','southeast');

bottom_text(btx,'pwd',1);
set(gcf,'position',[1149 386 946 946]);





