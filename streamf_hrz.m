% calculate streamfunctions
% for time-averaged fields
% Use barotropic flows, otherwise there can be 
% errors in streamfunctions
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;
f_extr = 0; % extract/calculate mean transport
% Barotropic or baroclinic U,V for calculating 
% transports and stremfunctions:
% use nlr to specify # of vert. layers for averaging
% if nlr > 100 - the code uses barotropic u
%  if u-barotropic is not saved in *.a file (check *.b)
% then will need to do averaging yourself:
% specify nlr=32 - max number of v. layers
nlr = 32; % strmfunc for nlr vert. layers (averaged); 
            % set some large number to calculate barotropic flow
rg  = 9806;  % convert pressure to depth, m


figure(1); clf;
%set(gcf,'Visible','off');

%YPLT=[2004,299;2005,351;2006,295;2007,354;2008,286;2009,236];
% Animation:
YRPLT=[];
cc=0;
for iyr=2004:2004
  for idd=1:7:365
%for iyr=2008:2008
%  for idd=1:7:365
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);


regn = 'ARCc0.08';
expt = '060';
%pthfig  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/fig_Greenland_river/';
pthfig  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/fig_BaffinBay/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat/';

% to get new transport:
fmt=sprintf('%smean_transp_Lr05_nrec0053.mat',pthmat);
ie=exist(fmt,'file');
if ~ie & f_extr==0
  fprintf('File does not exist %s\n',fmt);
  error('Need to extract data');
end

%IND = smaller_domain_indices('Green');
IND = smaller_domain_indices('NorthAtl');

inc1=IND.i1;
inc2=IND.i2;
jnc1=IND.j1;
jnc2=IND.j2;
djnc=IND.dj;
dinc=IND.di;

% Streamlines are calculated upto:
imx=1100;
jmx=1200;



ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

HH=HH(jnc1:jnc2,inc1:inc2);
LON=LON(jnc1:jnc2,inc1:inc2);
LAT=LAT(jnc1:jnc2,inc1:inc2);
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

if ~exist('DX','var')
  DX=zeros(mm,nn);
  DY=zeros(mm,nn);
  for i=1:nn-1
    dx=distance_spheric_coord(LAT(:,i),LON(:,i),LAT(:,i+1),LON(:,i+1));
    DX(:,i)=dx;
  end
  DX(:,nn)=dx;
  for j=1:mm-1
    dy=distance_spheric_coord(LAT(j,:),LON(j,:),LAT(j+1,:),LON(j+1,:));
    DY(j,:)=dy;
  end
  DY(mm,:)=dy;
end

% Get HYCOM fields:
% calculate transport/m: U*dH (not * by DX, DY): m2/s
cnc=0;
ip1=1;
Txmn=[];
Tymn=[];
ihh=find(HH>=0);

if f_extr>0 
  for ip=ip1:np
    yr=YRPLT(ip,1);
    iday=YRPLT(ip,2);
  %yr=2009;
  %iday=236;
    if yr<2008
      pthbin = sprintf(...
      '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/060/%4.4i/bin_outp/',yr);
    else
      pthbin = sprintf(...
      '/Net/yearly/ddmitry/hycom/ARCc0.08/060/%4.4i/bin_outp/',yr);
    end

  % Experiments:
  % 060 - CCMP winds, run 2004 -
  %       Greenland runoff on
  %       restart from old hycom expt 01.0 - 2006 jan.

  %  fina = sprintf('%s%s_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
  %  finb = sprintf('%s%s_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
    fina = sprintf('%s%s_archm.%4.4i_%3.3i_12r.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%s_archm.%4.4i_%3.3i_12r.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end

    cnc=cnc+1;
    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);

    fprintf('%4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

    if nlr>100  % barotropic flow
      [F,n,m] = read_hycom(fina,finb,'u_btrop');
      F=F(jnc1:jnc2,inc1:inc2);
      F(F>1e6)=nan;
      U=squeeze(F);

      [F,n,m] = read_hycom(fina,finb,'v_btrop');
      F=F(jnc1:jnc2,inc1:inc2);
      F(F>1e6)=nan;
      V=squeeze(F);
    else  % baroclinic velocities
      tic;
      [F,n,m,l] = read_hycom(fina,finb,'u-vel.');
      toc;
      F=F(:,jnc1:jnc2,inc1:inc2);
      F(F>1e6)=nan;
      U=squeeze(F);
      [F,n,m,l] = read_hycom(fina,finb,'v-vel.');
      F=F(:,jnc1:jnc2,inc1:inc2);
      F(F>1e6)=nan;
      V=squeeze(F);
      
      [F,n,m,l] = read_hycom(fina,finb,'thknss');
      F=F(:,jnc1:jnc2,inc1:inc2)/rg;
      F(F>1e10)=nan;
      dH=squeeze(F); % note this is not U,V thickness, center grid
		     % there are thck-u & thck-v 
      ih=find(dH<1e-1);
      U(ih)=nan;
      V(ih)=nan;
      
 % Check spurious U,V:
      I=find(abs(U)>2);
      U(I)=0; V(I)=0;
      I=find(abs(V)>2);
      U(I)=0; V(I)=0;
      
    end      
  % ------------
  % Check if u-vel. is
  % the total baroclinic velocity or baroclinic anomalies:
  % in instanteneous archive files
  % u-vel. = utot - ubtrop
  %
  % in archm u_vel=utot;
  % ------------
    f_chckU=0;
    if f_chckU>0
      [F,n,ml] = read_hycom(fina,finb,'u_btrop');
      F=F(1,jnc1:jnc2,inc1:inc2);
      F(F>1e6)=nan;
      Ub=squeeze(F);
      i=700;
      j=350;
      uu=squeeze(U(:,j,i));
      dh=squeeze(dH(:,j,i));
      uua=nansum(uu.*dh)./nansum(dh)
      ub=Ub(j,i)
    % uua and ub should be close:
      if abs(uua-ub)/abs(ub)<0.01;
	fprintf('u_vel. and v_vel. are total baroclinic U\n');
      else
	fprintf('u_vel. and v_vel. are barocl. anomlaies U\n');
      end
    end
  % ------------
  % Calculate transport:
  % in archm, archv: u, v anre not collocated 
  % with H, dH:
  % u(i,j) is between h(i-1,j) and h(i,j)
  % To check look at near land points
  % such that HH(j-1,i) and HH(j,i-1) is land and HH(j,i) is not
  % if V(j,i) is nan  than V is not collocated with H
  % if U(j,i) is nan than U is not collocated with H
    if nlr>100  % barotropic
      Tx=U*0;
      Ty=V*0;
      uu=U(:,2:imx);
      dh1=HH(:,1:imx-1);
      dh2=HH(:,2:imx);
      dh1(ihh)=0;
      dh2(ihh)=0;
      dh=0.5*(dh1+dh2);
      Tx=uu.*dh;
      
      vv=V(2:jmx,:);
      dh1=HH(1:jmx-1,:);
      dh2=HH(2:jmx,:);
      dh1(ihh)=0;
      dh2(ihh)=0;
      dh=0.5*(dh1+dh2);
      Ty=vv.*dh;
	
      Uav=U;
      Vav=V;
    else  % average over nlr layers:
      [ll,mm,nn]=size(U);
      Tx=U*0;
      Ty=V*0;
      for ii=2:imx
	uu=squeeze(U(:,:,ii));
	dh1=squeeze(dH(:,:,ii-1));
	dh2=squeeze(dH(:,:,ii));
	dh=0.5*(dh1+dh2);
	Tx(:,:,ii)=uu.*dh;
      end

      for jj=2:jmx
	vv=squeeze(V(:,jj,:));
	dh1=squeeze(dH(:,jj-1,:));
	dh2=squeeze(dH(:,jj,:));
	dh=0.5*(dh1+dh2);
	Ty(:,jj,:)=vv.*dh;
      end

      Txav=squeeze(nansum(Tx(1:nlr,:,:),1));
      Tyav=squeeze(nansum(Ty(1:nlr,:,:),1));
      Txav=Txav(1:jmx,1:imx);
      Tyav=Tyav(1:jmx,1:imx);
% Depth-mean U and V
      uh=U.*dH;
      vh=V.*dH;
      sdH=squeeze(nansum(dH(1:nlr,:,:),1));
      uh(ih)=nan;
      vh(ih)=nan;
      Uav=squeeze(nansum(uh(1:nlr,:,:),1))./sdH;
      Vav=squeeze(nansum(vh(1:nlr,:,:),1))./sdH;
      Uav=Uav(1:jmx,1:imx);
      Vav=Vav(1:jmx,1:imx);
    end

    
 % Time-mean fields:   
    if isempty(Txmn)
      Txmn=Txav*0;
      Tymn=Tyav*0;
      Umn=Txav*0;
      Vmn=Txav*0;
    end
    Txmn = Txmn+Txav;
    Tymn = Tymn+Tyav;
    Umn  = Umn+Uav;
    Vmn  = Vmn+Vav;
  end;
  Txmn=Txmn./cnc;
  Tymn=Tymn./cnc;

% save mean transport:
  Nlrs_av=nlr;
  fmt=sprintf('%smean_transp_Lr%2.2i_nrec%4.4i.mat',pthmat,nlr,cnc);
  fprintf('Saving %s\n',fmt);
  save(fmt,'Txmn','Tymn','Umn','Vmn','Nlrs_av');

else
  fprintf('Loading transport %s\n',fmt);
  load(fmt);
end

HHa = HH(1:jmx,1:imx);
DXa = DX(1:jmx,1:imx);
DYa = DY(1:jmx,1:imx);
LM  = HHa;
LM(LM>=0) = 0;
LM(LM<0)  = 1;

% Box-average:
db=0;
if db>0
  Tx0=Txmn*0;
  Ty0=Tymn*0;
  fprintf('Spacial filtering, # pnts %i\n',db*db);
  Txmn0=Txmn;
  Tymn0=Tymn;

  for ii=db+1:imx-db
    ni=ii/imx*100;
    if mod(ii,50)==0;
      fprintf('Filtering %i%% done...\n',round(ni));
    end
    for jj=db+1:jmx-db
      A=squeeze(Txmn(jj-db:jj+db,ii-db:ii+db));
      c=nanmean(nanmean(A));
      Tx0(jj,ii)=c;
      A=squeeze(Tymn(jj-db:jj+db,ii-db:ii+db));
      c=nanmean(nanmean(A));
      Ty0(jj,ii)=c;
    end
  end
  Txmn=Tx0;
  Tymn=Ty0;
end

i0=1; j0=1;
Psi = stream_fn(Txmn,Tymn,DXa,DYa,i0,j0,'simps',LM);
Psi(HHa>=0)=0;
db=5;
if db>0
  Psi0 = Psi;
  Psi  = sub_boxav(Psi,db);
end;

%keyboard

f_pltu=0;
if f_pltu,
%  Txav=squeeze(nansum(Tx(1:nlr,:,:),1));
%  Tyav=squeeze(nansum(Ty(1:nlr,:,:),1));
%  Tx=squeeze(nansum(U(1:nlr,:,:),1));
%  Ty=squeeze(nansum(V(1:nlr,:,:),1));
  Tx=zeros(mm,nn);
  Ty=zeros(mm,nn);
  Tx(1:jmx,1:imx)=Txmn;
  Ty(1:jmx,1:imx)=Tymn;
  res=1;
  stl='Mass Flux';
  nf=3;
  sub_plotu(LON,LAT,nf,HH,Tx,Ty,stl,res)
end

f_chck=0;
if f_chck
  A=Psi;
%  Psi=zeros(mm,nn);
%  Psi(1:jmx,1:imx)=A;
  Psi(HHa>0)=nan;

  figure(1); clf;
  contour(Psi,20,'b');
  hold on;
  contour(HH,[0 0],'k');
  Txav(HHa>0)=nan;
  Tyav(HHa>0)=nan;
  dm=10;
  x=[1:dm:nn];
  y=[1:dm:mm];

  scl=0.05;
  cf=0.25;
  beta=20;
  lwd=1.2;
  v_col=[0 0 0];
  dii=6;
  res=1;  % >0 have only unit vectors
  for ii=10:dii:500
    for jj=50:dii:850
      clear u v
      u = Txav(jj,ii);
      v = Tyav(jj,ii);
      s = sqrt(u*u+v*v);
      if isnan(u), continue; end;
      if res>0,
	u=u/s;
	v=v/s;
	scl=4;
      end

      x0=ii;
      y0=jj;

      x1=x0+u*scl;
      y1=y0+v*scl;
      draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
    end
  end


  axis('equal');
  set(gca,'xlim',[1 1000],'ylim',[1 1100]);

%  keyboard
end;



  
%  ffg=sprintf('%stracers_conc_%4.4i_%2.2i_%2.2i-layer%2.2i',pthfig,DV(1:3),plr);
%ffg=sprintf('%s060divU_layer%2.2i-%4.4i',pthfig,plr,cnc);

%  keyboard

if s_fig>0
  fprintf('Saving %s\n\n',ffg);
  print('-djpeg','-r200',ffg);
end


%keyboard
% Plot mean divU
for ik=1:12
  U=[];
  dU=DivU(ik).DIV;
  coeff=DivU(ik).coeff;
  nf=ik;
  DV=datevec(DivU(ik).TM(end));
  stl=sprintf('divU*%d %4.4i/%2.2i/%2.2i, Layer %i',coeff,DV(1:3),plr);
  sub_plot_divu(dU,LON,LAT,nf,HH,U,V,stl); 
  
  if s_fig
    ffg=sprintf('%s060divU_layer%2.2i-mo%2.2i',pthfig,plr,ik);
    fprintf('Saving %s\n',ffg);
    print('-djpeg','-r200',ffg);
  end
    
end





