% Calculate ocean heat flux to Greenland
% across specified contour -
% isobath around Greenland
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2008;
YR2 = 2010;
Cp = 4200; % J/kg K
%Tref= -273.15; % Ref T to calc. H flux
Tref= -1.8; % Ref T to calc. H flux

expt = 110;

s_mat = 1; % =0 - extract data, no save; =1 extract & save, =2 - load saved

rg=9806;  % convert pressure to depth, m
hgg=1e20; 

plr=0; % highlight this interface
btx = 'ocn_hflx_greenl008.m';


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
HFLX = struct;
HFLX.Info = 'Heat & Vol. Fluxes across 800m isobath around Greenland';
HFLX.Tref      = Tref; 
HFLX.GrCntr_II = GC.cntr_Iindx;
HFLX.GrCntr_JJ = GC.cntr_Jindx;
HFLX.Hbottom   = GC.Hbottom;
HFLX.DistCntr  = GC.Distance_m;

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
      fmat = sprintf('%s%3.3i_Greenl_HVflx_%i.mat',...
		   pthmat,expt,yrold);
      nlr = 41;
      HFLX = sub_io_hflx_month(HFLX, s_mat, mold, imo, fmat, nlr);
      mold = imo;
      yrold = yr;
    end
    
    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);
    
    tic;
    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
  %  F=F(:,jnc1:jnc2,inc1:inc2);
    F(F>hgg)=nan;
    T=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'salin');
    F(F>hgg)=nan;
    S=F;
%
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F(F>hgg)=0;
    U=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>hgg)=0;
    V=F;
%    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.','r_layer',34);
%    F(F>hgg)=0;
%    vmm=squeeze(F);

    fld='thknss';
    [F,n,m,l] = read_hycom(fina,finb,fld);
%    [F,n,m,l] = read_hycom(fina,finb,fld,'r_layer',34);
    F(F>1e18)=0;
    F=F/rg;
    F(F<1e-2)=0;
    dH = F;
    
    
%    [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
%    ZZ(isnan(ZZ))=100;
%    ZM(isnan(ZM))=100;

    Vct = sub_crrct_vflx(GC,U,V,dH,DX,DY);
%    Vct = 0;
% Process segments
    II = GC.cntr_Iindx;
    JJ = GC.cntr_Jindx;
    np = length(II);
    Hflx = zeros(nlr,np)*nan;
    ZZ = zeros(nlr,np)*nan;
    Vflx = zeros(1,np); % volume flux for checking
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
      clear hf* vf*
% collocate dH, T, U, V
Need to fix this:
u0=(u1*dh1+u2*dh2+...)/(dh1+dh2+...)
not 0.25*(u1*dh1+u2*dh2+...)
!!!!!!! - may be ok , double check?

      if xsct==0 % i2=i1
	if nx==0, error('Check norm x comp., it is 0!'); end;
	u1 = squeeze(U(:,j1,i1))*nx;
	dh1= 0.5*squeeze(dH(:,j1,i1-1)+dH(:,j1,i1));
	t1 = 0.5*squeeze( T(:,j1,i1-1)+ T(:,j1,i1));
	s1 = 0.5*squeeze( S(:,j1,i1-1)+ S(:,j1,i1));
	rho1= sw_dens0(s1,t1);
	dy1 = DY(j1,i1);
	u1(dh1<1e-3)=nan;
        u1  = sub_chck_uv(u1,0);
	u1  = u1+Vct;  % correction to have - net vol flux
	hf1 = Cp*rho1.*(t1-Tref).*u1.*dh1*dy1; % W=J/s
	vf1 = u1.*dh1.*dy1;                    % m3/s - volume flux
	
	u2 = squeeze(U(:,j2,i2))*nx;
	dh2= 0.5*squeeze(dH(:,j2,i2-1)+dH(:,j2,i2));
	t2 = 0.5*squeeze( T(:,j2,i2-1)+ T(:,j2,i2));
	s2 = 0.5*squeeze( S(:,j2,i2-1)+ S(:,j2,i2));
	rho2= sw_dens0(s2,t2);
	dy2 = DY(j2,i2);
	u2(dh2<1e-3) = nan;
        u2  = sub_chck_uv(u2,0);	
	u2  = u2+Vct;  % correction to have - net vol flux
	hf2 = Cp*rho2.*(t2-Tref).*u2.*dh2*dy2;
	vf2 = u2.*dh2.*dy2;                    % m3/s - volume flux
	
	u3 = squeeze(U(:,j2,i2+1))*nx;
	dh3= 0.5*squeeze(dH(:,j2,i2)+dH(:,j2,i2+1));
	t3 = 0.5*squeeze( T(:,j2,i2)+ T(:,j2,i2+1));
	s3 = 0.5*squeeze( S(:,j2,i2)+ S(:,j2,i2+1));
	rho3= sw_dens0(s3,t3);
	dy3 = DY(j2,i2+1);
	u3(dh3<1e-3) = nan;
        u3  = sub_chck_uv(u3,0);
	u3  = u3+Vct;  % correction to have - net vol flux
	hf3 = Cp*rho3.*(t3-Tref).*u3.*dh3*dy3;
	vf3 = u3.*dh3.*dy3;                    % m3/s - volume flux
	
	u4 = squeeze(U(:,j1,i1+1))*nx;
	dh4= 0.5*squeeze(dH(:,j1,i1)+dH(:,j1,i1+1));
	t4 = 0.5*squeeze( T(:,j1,i1)+ T(:,j1,i1+1));
	s4 = 0.5*squeeze( S(:,j1,i1)+ S(:,j1,i1+1));
	rho4= sw_dens0(s4,t4);
	dy4 = DY(j1,i1+1);
	u4(dh4<1e-3) = nan;
        u4  = sub_chck_uv(u4,0);	
	u4  = u4+Vct;  % correction to have - net vol flux
	hf4 = Cp*rho4.*(t4-Tref).*u4.*dh4*dy4;
	vf4 = u4.*dh4.*dy4;                    % m3/s - volume flux
	
	
	hflx = 0.25*(hf1+hf2+hf3+hf4);
	dZ   = 0.25*(dh1+dh2+dh3+dh4);
	vflx = 0.25*(vf1+vf2+vf3+vf4);
	
      else  % xsection
	if ny==0, error('Check norm y comp., it is 0!'); end;
	v1  = squeeze(V(:,j1,i1))*ny;
	dh1 = 0.5*(squeeze(dH(:,j1-1,i1)+dH(:,j1,i1)));
	t1  = 0.5*(squeeze( T(:,j1-1,i1)+ T(:,j1,i1)));
	s1  = 0.5*(squeeze( S(:,j1-1,i1)+ S(:,j1,i1)));
	rho1= sw_dens0(s1,t1);
	dx1 = DX(j1,i1);
	v1(dh1<1e-3) = nan;
        v1  = sub_chck_uv(v1,0);	
	v1  = v1+Vct;  % correction to have - net vol flux
	hf1 = Cp*rho1.*(t1-Tref).*v1.*dh1*dx1;
	vf1 = v1.*dh1.*dx1;                    % m3/s - volume flux
	
	v2  = squeeze(V(:,j1+1,i1))*ny;
	dh2 = 0.5*(squeeze(dH(:,j1,i1)+dH(:,j1+1,i1)));
	t2  = 0.5*(squeeze( T(:,j1,i1)+ T(:,j1+1,i1)));
	s2  = 0.5*(squeeze( S(:,j1,i1)+ S(:,j1+1,i1)));
	rho2= sw_dens0(s2,t2);
	dx2 = DX(j1+1,i1);
	v2(dh2<1e-3) = nan;
        v2  = sub_chck_uv(v2,0);	
	v2  = v2+Vct;  % correction to have - net vol flux
	hf2 = Cp*rho2.*(t2-Tref).*v2.*dh2*dx2;
	vf2 = v2.*dh2.*dx2;                    % m3/s - volume flux
	
	v3  = squeeze(V(:,j2+1,i2))*ny;
	dh3 = 0.5*(squeeze(dH(:,j2,i2)+dH(:,j2+1,i2)));
	t3  = 0.5*(squeeze( T(:,j2,i2)+ T(:,j2+1,i2)));
	s3  = 0.5*(squeeze( S(:,j2,i2)+ S(:,j2+1,i2)));
	rho3= sw_dens0(s3,t3);
	dx3 = DX(j2+1,i2);
	v3(dh3<1e-3) = nan;
        v3  = sub_chck_uv(v3,0);	
	v3  = v3+Vct;  % correction to have - net vol flux
	hf3 = Cp*rho3.*(t3-Tref).*v3.*dh3*dx3;
	vf3 = v3.*dh3.*dx3;                    % m3/s - volume flux
	
	v4  = squeeze(V(:,j2,  i2))*ny;
	dh4 = 0.5*(squeeze(dH(:,j2-1,i2)+dH(:,j2,i2)));
	t4  = 0.5*(squeeze( T(:,j2-1,i2)+ T(:,j2,i2)));
	s4  = 0.5*(squeeze( S(:,j2-1,i2)+ S(:,j2,i2)));
	rho4= sw_dens0(s4,t4);
	dx4 = DX(j2,i2);
	v4(dh4<1e-3) = nan;
        v4  = sub_chck_uv(v4,0);	
	v4  = v4+Vct;  % correction to have - net vol flux
	hf4 = Cp*rho4.*(t4-Tref).*v4.*dh4*dx4; % 
	vf4 = v4.*dh4.*dx4;                    % m3/s - volume flux
	
	hflx = 0.25*(hf1+hf2+hf3+hf4);
	dZ   = 0.25*(dh1+dh2+dh3+dh4);
	vflx = 0.25*(vf1+vf2+vf3+vf4);
	
      end
      
      Hflx(:,ig) = hflx;
      zZ         = -cumsum(dZ);
      zZ(dZ==0)  = nan;
      ZZ(:,ig)   = zZ;
      Vflx(ig) = nansum(vflx); % depth-integrated vol flux, m3/s
      
    end  % ig - segment

    fprintf('VFlux %4.2f Sv, HFlux %5.3d W\n',...
	    nansum(Vflx)*1e-6, nansum(nansum(Hflx)));
    
    fprintf('1 day processed %6.3f min\n\n',toc/60);
    
    HFLX(imo).nrec         = HFLX(imo).nrec+1;
    dmm                    = HFLX(imo).Vol_flux_m3s;
    HFLX(imo).Vol_flux_m3s = dmm+Vflx;
    dmm                    = HFLX(imo).Hflux_W;
    HFLX(imo).Hflux_W      = dmm+Hflx;
    dmm                    = HFLX(imo).ZZ;
    HFLX(imo).ZZ           = dmm+ZZ;
    HFLX(imo).TM           = dnmb;

  end  % day
    
end   % year

fmat = sprintf('%s%3.3i_Greenl_HVflx_%i.mat',...
		   pthmat,expt,yrold);
HFLX = sub_io_hflx_month(HFLX, s_mat, mold, imo, fmat, nlr);



    
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
      cff = 1e-12; 
%      Wm = Hflx./ar*cff; % W -> W/m2 -> MW/m2
%      tstr=sprintf('HFlx, W/m2*%3.1d, %s',cff,datestr(dnmb));
      Wm = Hflx*cff; % stay in Watts
      tstr=sprintf('HFlx, W*%3.1d, %s',cff,datestr(dnmb));
% For plotting add surface layer      
      ZZp = [ZZ(1,:);ZZ];
      ZZp(1,:)=0;
      Wmp  = [Wm(1,:); Wm];
      Dstp = [Dst(1,:);Dst];
      sub_plot_hflxZcntr(Wmp,ZZp,Dstp,fnmb,tstr,Hs);
      
    end
    



