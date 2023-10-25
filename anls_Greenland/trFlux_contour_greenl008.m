% Calculate ocean FW Tracer flux across a closed contour
% around Greenland
% isobath around Greenland
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2007;
YR2 = 2007;
Sref = 34.8;

expt = 110;

s_mat = 1; % =0 - extract data, no save; =1 extract & save, =2 - load saved

rg=9806;  % convert pressure to depth, m
hgg=1e20; 

plr=0; % highlight this interface
btx = 'fwflux_contour_greenl008.m';

fprintf('Tracer Flux, Greenland Section, %i-%i\n',YR1,YR2);


regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthfig=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_green_xsct/',...
		  regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section

% ================================================
% Plot Greenland map and the contour
% ================================================
f_pltgr=0;
if f_pltgr==1
  Hg = HH;

  Hg(1:380,:)=nan;
  Hg(1100:end,:)=nan;
  Hg(:,1050:end)=nan;
  Hg(:,1:450)=nan;

  fprintf('Plotting Greenland heat contour ...\n');
  fn=2;
  domname = '0';
  sub_plot_bath(Hg,LON,LAT,fn,domname);
  contour(Hg,[-100:10:-5],'Color',[0.85 0.85 0.85]);
  contour(Hg,[-3500:100:-10],'Color',[0.6 0.6 0.6]);
  contour(Hg,[-4000:500:-50],'Color',[0.25 0.25 0.25]);
  IIs = GC.cntr_Iindx;
  JJs = GC.cntr_Jindx;
  x   = GC.Distance_m*1e-3; % m->km

  plot(IIs,JJs,'b-','Linewidth',2);
  for km=0:500:max(x)
    d=abs(x-km);
    i0=find(d==min(d));
    if km==0
       plot(IIs(i0),JJs(i0),'r.','Markersize',14);
       plot(IIs(i0),JJs(i0),'rd','Markersize',6);
    else
      plot(IIs(i0),JJs(i0),'r.','Markersize',11);
    end
  end

  set(gca,'xlim',[450 1050],...
	  'ylim',[380 1100],...
	  'xtick',[],...
	  'ytick',[]);
  
  title('Contour (~800m) for heat flux calculation');

  bottom_text(btx,'pwd',1);
end

% ================================================
FWFLX = struct;
FWFLX.Info = 'Tracer Fluxes across 800m isobath around Greenland';
FWFLX.GrCntr_II = GC.cntr_Iindx;
FWFLX.GrCntr_JJ = GC.cntr_Jindx;
FWFLX.Hbottom   = GC.Hbottom;
FWFLX.DistCntr  = GC.Distance_m;

mold  = 0;
yrold = YR1;
dday  = 6; 
for iyr=YR1:YR2
  cc=0;
  TM = [];
  yr = iyr;
  
  for iday = 1:dday:365
    pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
    cc   = cc+1;
    j1d  = datenum(yr,1,1);
    dnmb = j1d+iday-1;
    DV   = datevec(dnmb);
    imo  = DV(2);
    TM(cc,1) = dnmb;

% ==================
% Save monthly means
% ==================
    if mold~=imo
      fmat = sprintf('%s%3.3i_Greenl_TrFlx_%i.mat',...
		   pthmat,expt,yrold);
      nlr = 41;
      FWFLX = sub_io_fwflx_month(FWFLX, s_mat, mold, imo, fmat, nlr);
      mold = imo;
      yrold = yr;
    end
    
    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);
    
    tic;
%    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
%    F(F>hgg)=nan;
%    T=F;

%    [F,n,m,nlr] = read_hycom(fina,finb,'salin');
%    F(F>hgg)=nan;
%    S=F;
%
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F(F>hgg)=0;
    U=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>hgg)=0;
    V=F;

    fld='thknss';
    [F,n,m,l] = read_hycom(fina,finb,fld);
%    [F,n,m,l] = read_hycom(fina,finb,fld,'r_layer',34);
    F(F>1e18)=0;
    F=F/rg;
    F(F<1e-2)=0;
    dH = F;

    nTr = 1;
    [F,nn,mm,ll] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
    F(F>1e6)=nan;
    Tr = F;
    
%    [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
%    ZZ(isnan(ZZ))=100;
%    ZM(isnan(ZM))=100;

%    Vct = sub_crrct_vflx(GC,U,V,dH,DX,DY);
    Vct = 0;
% Process segments
    II = GC.cntr_Iindx;
    JJ = GC.cntr_Jindx;
    np = length(II);
    Hflx = zeros(nlr,np)*nan;
    ZZ = zeros(nlr,np)*nan;
%    Vflx = zeros(1,np); % volume flux for checking
    for ig=1:np-1
      i1=II(ig);
      i2=II(ig+1);
      j1=JJ(ig);
      j2=JJ(ig+1);
      
      if HH(j1,i1)>=0 & HH(j2,i2)>=0, 
	continue; 
      end;
      
      if Hs(ig)>0, continue; end;
      
      xsct = 1;
      if i2==i1 % y section
	xsct = 0;
      end
% Reorder i1/j1 and i2/j2 so that i2>i1 or j1>j2
      if i2<i1, 
	imm=i2;
	i2=i1;
	i1=imm;
      end
      if j2<j1
	jmm=j2;
	j2=j1;
	j1=jmm;
      end

      nrm = GC.Norm_in(ig,:);
      nx  = nrm(1);
      ny  = nrm(2);
      dZ = 0;
      
      clear dh*

% Fluxes calculated:
%            v(i,j+1)        v(i+1,j+1)
%      |--------|-------|-------|---------|
%      |                |                 |
%      |                |                 |
%      |     V1         |       V2        |  Flx = Cp*rho*(T-Tref)V1*dH1*dX1/2+
%      -        *===============*         |        Cp*rho*(T-Tref)V2*dH2*dX2/2+
%u(i,j)|     T,S,dH     |                 |   where V1 and V2 are flux-averaged
%      |                |                 |   v(i,j+1)&v(i,j) and v(i+1,j+1)&v(i+1,j)
%      |                |                 |
%      |--------|-------|-------|---------|
%      |      v(i,j)         v(i+1,j)
%      |   
      
      if xsct==0 % i2=i1
	if nx==0, error('Check norm x comp., it is 0!'); end;
% U Fluxes in grid cell (j1,i1):	
	u1 = squeeze(U(:,j1,i1))*nx;
	dh1= 0.5*squeeze(dH(:,j1,i1-1)+dH(:,j1,i1));
	u1(dh1<1e-3)=nan;
        u1  = sub_chck_uv(u1,dh1,0);
	
	u2 = squeeze(U(:,j1,i1+1))*nx;
	dh2= 0.5*squeeze(dH(:,j1,i1)+dH(:,j1,i1+1));
	u2(dh2<1e-3) = nan;
        u2  = sub_chck_uv(u2,dh2,0);	

% Values in the center of the grid cell	
	u12 = (u1.*dh1+u2.*dh2)./(dh1+dh2);
	dh12= squeeze(dH(:,j1,i1));
	t12 = squeeze(Tr(:,j1,i1));
	dy12 = DY(j1,i1);
	hf12 = t12.*u12.*dh12*dy12; % kg/m3*m/2*m2 => kg/s
	
% U Fluxes in grid cell (j2,i2):	
	u3  = squeeze(U(:,j2,i2))*nx;
	dh3 = 0.5*squeeze(dH(:,j2,i2-1)+dH(:,j2,i2));
        u3  = sub_chck_uv(u3,dh3,0);
	
	u4  = squeeze(U(:,j2,i2+1))*nx;
	dh4 = 0.5*squeeze(dH(:,j2,i2)+dH(:,j2,i2+1));
        u4  = sub_chck_uv(u4,dh4,0);
	
% Values in the center of the grid cell	
	u34 = (u3.*dh3+u4.*dh4)./(dh3+dh4);
	dh34= squeeze(dH(:,j2,i2));
	t34 = squeeze(Tr(:,j2,i2));
	dy34 = DY(j2,i2);
	hf34 = t34.*u34.*dh34*dy34; % W=J/s
	
	hflx = 0.5*(hf12+hf34);
	dZ   = 0.5*squeeze(dH(:,j1,i1)+dH(:,j2,i2));
	
      else  % xsection
	if ny==0, error('Check norm y comp., it is 0!'); end;
	v1  = squeeze(V(:,j1,i1))*ny;
	dh1 = 0.5*(squeeze(dH(:,j1-1,i1)+dH(:,j1,i1)));
	v1(dh1<1e-3) = nan;
        v1  = sub_chck_uv(v1,dh1,0);	
	
	v2  = squeeze(V(:,j1+1,i1))*ny;
	dh2 = 0.5*(squeeze(dH(:,j1,i1)+dH(:,j1+1,i1)));
	v2(dh2<1e-3) = nan;
        v2  = sub_chck_uv(v2,dh2,0);	

% Values in the center of the grid cell	
	v12   = (v1.*dh1+v2.*dh2)./(dh1+dh2);
        dh12  = squeeze(dH(:,j1,i1));	
	t12   = squeeze(Tr(:,j1,i1));
	dx12  = DX(j1,i1);
	hf12  = t12.*v12.*dh12*dx12;

	v3  = squeeze(V(:,j2,i2))*ny;
	dh3 = 0.5*(squeeze(dH(:,j2,i2)+dH(:,j2-1,i2)));
	v3(dh3<1e-3) = nan;
        v3  = sub_chck_uv(v3,dh3,0);	

	v4  = squeeze(V(:,j2+1,i2))*ny;
	dh4 = 0.5*(squeeze(dH(:,j2,i2)+dH(:,j2+1,i2)));
	v4(dh4<1e-3) = nan;
        v4  = sub_chck_uv(v4,dh4,0);	
	
% Values in the center of the grid cell	
        v34 = (v3.*dh3+v4.*dh4)./(dh3+dh4);
	dh34= squeeze(dH(:,j2,i2));	
	t34  = squeeze(Tr(:,j2,i2));
	dx34 = DX(j2+1,i2);
	hf34 = t34.*v34.*dh34*dx34;
	
	hflx = 0.5*(hf12+hf34);
	dZ   = 0.5*squeeze(dH(:,j1,i1)+dH(:,j2,i2));

      end
      
      Hflx(:,ig) = hflx;
      zZ         = -cumsum(dZ);
      zZ(dZ==0)  = nan;
      ZZ(:,ig)   = zZ;
      
    end  % ig - segment

    fprintf('Tracer Flux %5.3d ton/s\n',...
	    nansum(nansum(Hflx))*1e-3);
    
    fprintf('1 day processed %6.3f min\n\n',toc/60);
    
    FWFLX(imo).nrec         = FWFLX(imo).nrec+1;
    dmm                     = FWFLX(imo).TrFlux_kg_s;
    FWFLX(imo).TrFlux_kg_s  = dmm+Hflx;
    dmm                     = FWFLX(imo).ZZ;
    FWFLX(imo).ZZ           = dmm+ZZ;
    FWFLX(imo).TM           = dnmb;

  end  % day
    
end   % year

fmat = sprintf('%s%3.3i_Greenl_TrFlx_%i.mat',...
		   pthmat,expt,yrold);
FWFLX = sub_io_fwflx_month(FWFLX, s_mat, mold, imo, fmat, nlr);


f_pflx = 0;
if f_pflx==1
  fnmb = 10;
  dx = GC.Distance_m*1e-3; % m->km
  Hs = GC.Hbottom;
 % plot(dx,Hs); % plot Bottom profile along contour
 % set(gca,'xlim',[0 max(dx)],'xtick',[0:500:max(dx)])
 % title('Bottom profile along the heat contour');
% Convert W -> W/m2
  ar = ZZ*0;
  [Dst,dmb] = meshgrid(dx,[1:nlr]);
  ar = abs(Dst*1e3.*ZZ); % m2
  cff = 1e-3;  % kg -> tons 
%      Wm = Hflx./ar*cff; % W -> W/m2 -> MW/m2
%      tstr=sprintf('HFlx, W/m2*%3.1d, %s',cff,datestr(dnmb));
  Wm = Hflx*cff; %
  tstr=sprintf('TrFlux, kg/s*%3.1d, %s',cff,datestr(dnmb));
% For plotting add surface layer      
  ZZp = [ZZ(1,:);ZZ];
  ZZp(1,:)=0;
  Wmp  = [Wm(1,:); Wm];
  Dstp = [Dst(1,:);Dst];
  sub_plot_hflxZcntr(Wmp,ZZp,Dstp,fnmb,tstr,Hs);

end




