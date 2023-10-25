% Calculate and Plot Eddy flux across sections
% Eddy fluxes are calculated based on 
% Nagai et al., JGR, 2015, "Dominant role of eddies ..."
%
% The mean and perturabations (anomalies) are 
% calcualted in monthly_mean_flds.m and calc_UC_eddy_transp.m
% 
% Note this is different approach from
% tracer_eddy_flux.m
%
% 
% Note that HYCOM files are subset files: smaller region
%  and fewer variables (look for nAtl in the name)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig=0;
s_mat=0; % =0 - plot fluxes;
          %=1 - save mat; =2 - calculate fluxes, do not save

pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_monthly_mean/';
pthmat2  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig  = '/Net/tholia/ddmitry/hycom/ARCc0.08/fig_TrEddy/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);
txtbtm='/hycom_arc08/Greenland_expt/tracer_eddy_sectionsNagai.m';


YRPLT=[];
cc=0;
for iyr=2004:2008
  for im=1:12
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=im;
  end
end
npp=length(YRPLT);
yr1=YRPLT(1,1);
yr2=YRPLT(end,1);

fmatFlx = sprintf('%stracEddyFlx_sect_Nagai-%i%i.mat',pthmat,yr1,yr2);


ftopo = sprintf('%sdepth_ARCc0.08_09.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

vrs=3; 
fmat = sprintf('%sGrCoast_sections.mat',pthmat2); % Greenland coast sections

fprintf('\nLoading sections %s\n\n',fmat);
load(fmat);
inc1=IND.i1;
inc2=IND.i2;
jnc1=IND.j1;
jnc2=IND.j2;
djnc=IND.dj;
dinc=IND.di;

% Collocate sub-regions, find index offset
INDa = smaller_domain_indices('NorthAtl');
inc1a=INDa.i1;
inc2a=INDa.i2;
jnc1a=INDa.j1;
jnc2a=INDa.j2;
djnca=INDa.dj;
dinca=INDa.di;

if inc1~=inc1a
  error('Need to modify code: add index offset for different regions');
elseif jnc1~=jnc1a
  error('Need to modify code: add index offset for different regions');
end
%keyboard

HH=HH(jnc1:jnc2,inc1:inc2);
LON=LON(jnc1:jnc2,inc1:inc2);
LAT=LAT(jnc1:jnc2,inc1:inc2);
[mm,nn]=size(LON);

%if s_mat<0
%  sub_plot_eddy_flux(flxmat);
%  return
%end

% -------------
% Plot sections
% --------------
f_sct=0;
if f_sct>0
  figure(10); clf;
  hold on;
  contour(HH,[0 0],'Color',[0.4 0.4 0.4]);
  contour(HH,[-8000:500:-100],'Color',[0.8 0.8 0.8]);
  contour(HH,[-1000:100:-10],'Color',[0.7 0.7 0.7]);
  axis('equal');
  set(gca,'xlim',[100 700],'ylim',[250 1000]);

  ns=length(SEGM);
  for ir=1:ns
    nm=SEGM(ir).Name;
    I=SEGM(ir).I;
    J=SEGM(ir).J;

% Plot all points:    
%    plot(I,J,'b.','linewidth',2);
%    plot(I,J,'b.');
% Plot segments:
    plot([I(1) I(end)],[J(1) J(end)],'b-','linewidth',2);
    plot(I(1),J(1),'bo','linewidth',2);
    plot(I(end),J(end),'bo','linewidth',2);
    ttx=sprintf('%2.2i',ir);
    x0=mean(I);
    y0=mean(J);
%    text(x0,y0,ttx,'Fontsize',12);
  end
  set(gca,'xtick',[],'ytick',[]);
end

% Distances along the sections:
ns=length(SEGM);
for ir=1:ns
  nm=SEGM(ir).Name;
  X=SEGM(ir).LON;
  Y=SEGM(ir).LAT;
  ll=length(X);
  clear D
  for l=1:ll-1
    x1=X(l);
    y1=Y(l);
    x2=X(l+1);
    y2=Y(l+1);
    D(l,1)=distance_spheric_coord(y1,x1,y2,x2);
  end
  SEGM(ir).Dist=D;
end

% Plot fluxes:
if s_mat == 0,
  sub_plot_coast_eddyFlx(fmatFlx, SEGM, s_fig, pthfig);
  return
end





cnc=0;
ip1=1;


for ir=1:ns
  AV(ir).Ua  = [];
  AV(ir).UCa = [];
  AV(ir).Ca  = [];
  AV(ir).Ua0 = [];
  AV(ir).UCa0= [];
  AV(ir).Ca0 = [];
end;


% Calculate fluxes across secitons:
cc=0;
for ip=1:npp
  yr=YRPLT(ip,1);
  im=YRPLT(ip,2);
  icyc=2;
% Eddy fluxes by U and V components
% integrated over Nlev:
  fout = sprintf('%smonth_UC_eddy_%4.4i%2.2i-%i.mat',pthmat,yr,im,icyc);
  fprintf('Loading %s\n',fout);
  load(fout);
  
  TM=UCPR.TM;
  nd=length(TM);
  Nlev=UCPR.Nlev;

% Load precalculated monthly means:
  fmean = sprintf('%smonth_UVCSdP_%4.4i%2.2i-%i.mat',...
		    pthmat,yr,im,icyc);
  fprintf('Loading mean: %s\n',fmean);
  load(fmean);
  Umn = squeeze(FMEAN.Uav(1:Nlev,:,:));
  Vmn = squeeze(FMEAN.Vav(1:Nlev,:,:));
  TRmn= squeeze(FMEAN.Trac_av(1:Nlev,:,:));
  dHmn= squeeze(FMEAN.dHav(1:Nlev,:,:));
  Smn = squeeze(FMEAN.Sav(1:Nlev,:,:));
  clear FMEAN;
  TRmn(TRmn<0)=0;
  UmTRm=squeeze(sum(Umn.*dHmn.*TRmn));
  VmTRm=squeeze(sum(Vmn.*dHmn.*TRmn));

% Anomalies/perturbations averaged over the upper Nlev
  UpTRp=UCPR.UpTRp_Nlev;
  VpTRp=UCPR.VpTRp_Nlev;


  [nd,a1,a2]=size(UpTRp);

  
% Calculate transport
% Sign convention:
% U+ - positive X (right)
% V+ - positive Y (up)
%keyboard
% Monthly average eddy fluxes
%  for iday=1:nd
  iday=15;
    tic;    
    dnmb=datenum(yr,im,iday);
    cc=cc+1;
    
    ns=length(SEGM);
    for ir=1:ns    % all sections
      fprintf('Calculating fluxes, section %i\n',ir);
      nm = SEGM(ir).Name;
      I  = SEGM(ir).I;
      J  = SEGM(ir).J;
      DX = SEGM(ir).Dist;
      nx = SEGM(ir).Norm(1); % positive towards Greenl. coast 
      ny = SEGM(ir).Norm(2); %

  % For flux calculation in boxes
  % important to keep right direction of segments
  % to close contours, i and j should increase
      [b1,b2]=size(I);
      if b1<b2, I=I'; end;
      [b1,b2]=size(J);
      if b1<b2, J=J'; end;

      if J(end)==J(1),
	if I(end)<I(1)
	  I=flipud(I);
	  J=flipud(J);
	end
      end
      if I(end)==I(1)
	if J(end)<J(1)
	  I=flipud(I);
	  J=flipud(J);
	end
      end

      ll=length(I);
      clear FLX
      UpCp=[];
      UmCm=[];
      U =[];
      C =[];
      DH=[];
  %   keyboard
  %    parfor l=1:ll-1   % little segments of a section
      for l=1:ll-1   % little segments of a section
	i1=I(l);
	i2=I(l+1);
	j1=J(l);
	j2=J(l+1);
	dd=DX(l);  % segment length, m
  % Find Normal:
  % Sign convention:
	di=abs(i2-i1);
	dj=abs(j2-j1);
	if di>0 & dj>0
	  error('Not step-like section ...');
	end

  % Check norms:
	if (I(1)-I(end))==0 & ny~=0 % vertical straight section
	  error('Check norms: ny not 0, di=0');
	elseif (J(1)-J(end))==0 & nx~=0
	  error('Check norms: nx not 0, dj=0');
	end

	flx = []; 
	if di==0  % Y-section, V*norm=0
	  UTp   = squeeze(nanmean(UpTRp(:,j1,i1))); % time average	
	  flx   = nx*UTp*dd; % UTp=[kg/m3*m*m/s]*dd=[m]=kg/s->flux over Nlev depth
	  flxMn = nx*UmTRm(j1,i1)*dd; % mean non-eddying flux
	elseif dj==0
	  VTp   = squeeze(nanmean(VpTRp(:,j1,i1)));
	  flx   = ny*VTp*dd; % UTp=[kg/m3*m*m/s]*dd=[m]=kg/s->flux over Nlev depth
	  flxMn = ny*VmTRm(j1,i1)*dd; % mean non-eddying flux
	end

	UpCp(l,1) = flx; % eddy flux across a grid-cell side
	UmCm(l,1) = flxMn; 
	DH(l,1)   = HH(j1,i1);
      end;  % little segments along a section loop

  % Sum over the section and only over ocean cells
%	Ih=find(DH<0);
	UmCm_intgr = nansum(UmCm); % intgr_section(<U>*<C>) - regular flux
	UpCp_intgr = nansum(UpCp); % eddy flux<U'*C'>, intger over section, z=0

  %    keyboard
      if ip==1
	TRFLX(ir).Info='Tracer eddy fluxes';
	TRFLX(ir).Code='tracer_eddy_sectionsNagai.m';
	TRFLX(ir).Name=nm;
	TRFLX(ir).Norm=[nx,ny];
	TRFLX(ir).Nlev_avrg=Nlev;
      end;
     
      TRFLX(ir).Nav_averaged_days= 31; 
      TRFLX(ir).TM(cc)           = dnmb; % 
      TRFLX(ir).CYC(cc)          = icyc;
      TRFLX(ir).RegFlux_UmCm(cc)  = UmCm_intgr; % <U>*<C>
      TRFLX(ir).EddyFlux_UpCp(cc) = UpCp_intgr; % <U'>*<C'>

      fprintf('Section %i, Reg.Flux=%7.3d EddyFlux=%7.3d\n',...
	      ir,UmCm_intgr,UpCp_intgr);
%      keyboard
    end;    % x-sections

    fprintf('1 day: %7.1f sec\n',toc);
%  end;  % days  
  
  if s_mat==1
    fprintf('Saving %s\n',fmatFlx);
    save(fmatFlx,'TRFLX');
  end
  
  
end % months






