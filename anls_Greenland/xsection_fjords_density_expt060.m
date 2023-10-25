% Plot vertical disctribution
% of densty (brunt-vias. freq) in fjords
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat = 0; % =0 - extract data, no save; =1 extract & save, =2 - load saved
s_fig = 1;
f_tracloc = 0; % =1 - add locations of tracer release to transect maps
rg=9806;  % convert pressure to depth, m
fld0='brunt'; %'salin' - background field
plr=0; % highlight this interface
txtb = 'ddmitry@mars:hycom_NOPP_rivers/anls_Greenland/xsection_fjords_density.m';


regn = 'ARCc0.08';
expt = '060';
pthfig  = '/Net/tholia/ddmitry/hycom/ARCc0.08/fig_Grtrcr_fjord/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat/';
ntopo   = 9;
fmat=sprintf('%svert_vel_layers.mat',pthmat);

PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/relax/';
PTH.topo=pthtopo;
PTH.river='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.mat = pthmat;
PTH.fig = pthfig;


domname = 'NorthAtl';
%IND = smaller_domain_indices('Green');
IND = smaller_domain_indices(domname);
%IND = smaller_domain_indices('Baffin');

inc1=IND.i1;
inc2=IND.i2;
jnc1=IND.j1;
jnc2=IND.j2;
djnc=IND.dj;
dinc=IND.di;

YRPLT=[];
cc=0;
iyr2=2010;
for icyc=1:1
  for iyr=2008:2008
%    for idd=85:84:365
    for idd=253:84:253
      if idd==1, idd=2; end;
      cc=cc+1;
      YRPLT(cc,1)=iyr;
      YRPLT(cc,2)=idd;
      YRPLT(cc,3)=icyc;
    end
  end
end

%YRPLT=[2008,50];
np=size(YRPLT,1);



ftopo = sprintf('%s/depth_%s_%2.2i.nc',pthtopo,regn,ntopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

HH=HH(jnc1:jnc2,inc1:inc2);
LON=LON(jnc1:jnc2,inc1:inc2);
LAT=LAT(jnc1:jnc2,inc1:inc2);
[mm,nn]=size(LON);

%fn=1;
%sub_plot_bath(HH,LON,LAT,fn,domname);

% Specify segments of the x-section
xname='uummq'; % Fjord Uummannaq - western coast
%xname='scoresby'; % Scoresby Sound - eastern coast
xdr=1;
switch(xname);
 case('uummq');
  IJs=[266, 639;
      244, 660;
      196, 678];
  xdr=-1; % reverse x dir
 case('scoresby');
  IJs=[   481   562
   500   563
   549   545];
end
IIa=IJs(:,1);
JJa=IJs(:,2);

nij=size(IJs,1);
IIs=[];
JJs=[];
for ii=1:nij-1
  i1=IJs(ii,1);
  i2=IJs(ii+1,1);
  j1=IJs(ii,2);
  j2=IJs(ii+1,2);
  [I,J]=sub_xsct_indx(i1,j1,i2,j2);
  if size(I,1)==1;
    I=I';
    J=J';
  end
  
  IIs=[IIs;I];
  JJs=[JJs;J];
end;

IJs=[IIs,JJs];

nS=length(IIs);
clear Xl Yl
for ii=1:nS
  i0=IJs(ii,1);
  j0=IJs(ii,2);
  Xl(ii,1)=LON(j0,i0);
  Yl(ii,1)=LAT(j0,i0);
end;
INDs=sub2ind(size(HH),JJs,IIs);

f_map=0;
if f_map>0
  figure(10); clf;
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-5000:500:-100],'Color',[0.7 0.7 0.7]);
  plot(IIs,JJs,'b.-');
  set(gca,'xlim',[100 600],...
	  'ylim',[200 1000]);
end

fprintf('Section: %s, Saving fig: %i\n',xname,s_fig);
%keyboard

% Plot fields:
cnc=0;  % 
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  icyc=YRPLT(ip,3);
  switch (icyc),
   case(1)
    expt='060';
   case(2),
    expt='061';
  end
%yr=2009;
%iday=236;
    pthbin = sprintf(...
    '/Net/yearly/ddmitry/hycom/ARCc0.08/%s/%4.4i/bin_outp/',expt,yr);


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

  tic;
  [F,n,nlev] = read_hycom(fina,finb,'temp');
  toc;
  F=F(:,jnc1:jnc2,inc1:inc2);
  F(F>1e6)=nan;

  T=squeeze(F(:,INDs));

  [F,n,nlev] = read_hycom(fina,finb,'salin');
  F=F(:,jnc1:jnc2,inc1:inc2);
  F(F>1e6)=nan;

  S=squeeze(F(:,INDs));

  Rho = sw_dens0(S,T);
  
%  keyboard
%  lTr=log(Tr);
% Prepare for plotting - add extra bogus layer
% at the bottom, Matlab won't show it
  [a1,a2]=size(Rho);
  
% Prepare vertical thicknesses for plotting  
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F=F(:,jnc1:jnc2,inc1:inc2);
  F(F>1e10)=nan;
  F(F<0.1)=nan;
  F=F/rg;
  Dsec=squeeze(F(:,INDs));
  Dsec(Dsec==0)=nan;
% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
  clear ZZb
  Dsec(isnan(Dsec))=0;
  ZZb(1,:)=-Dsec(1,:);
  for kk=2:l
    ZZb(kk,:)=ZZb(kk-1,:)-Dsec(kk,:);
  end
% For plotting need to have #of layers +1 interfaces
% add surface layer, otherwise all values
% will be shifted 1 layer down
  [nl,npb]=size(ZZb);
  ZZ=zeros(nl+1,npb);
  ZZ(2:nl+1,:)=ZZb;
  ZZav=ZZ;
  
% Depths of middle of the cells:
  ZM(1,:)=0.5*ZZ(1,:);
  for kk=1:l
    ZM(kk,:)=0.5*(ZZ(kk+1,:)+ZZ(kk,:));
  end

% Interpolate into reg. vert. grid:
  ZZi=[0:-5:-1000]';
  nli=length(ZZi);
  Rhoi=zeros(nli,npb)*nan;
  for ipp=1:npb
    Ib=min(find(Dsec(:,ipp)==0));
    if Ib==1, continue; end;
    rho=Rho(:,ipp);
    zz=ZM(:,ipp);
    hb=-sum(Dsec(:,ipp));
    ibz = max(find(ZZi>=hb));
    for kl=Ib:nl
      zz(kl)=zz(kl-1)-0.1;
    end;
    zz=[0;zz];
    rho=[rho(1);rho];
    if zz(nl)>ZZi(end)
      zz(nl)=ZZi(end);
    end
    rhoi = interp1(zz,rho,ZZi,'cubic');
    rhoi(ibz+1:end)=nan;
    Rhoi(:,ipp)=rhoi;
  end
  
  
% Calculate Brunt-Vaisala freq.  
  g=9.8;
  rho0=1027;
  for k=2:nli-1
    drho = Rhoi(k-1,:)-Rhoi(k+1,:);
    drho(drho>=0)=-1e-8; % unstable strat. or mixed
    dz   = ZZi(k-1,:)-ZZi(k+1,:);
    dz(dz==0)=nan;
    nfr=-g/rho0.*drho./dz;
    NFR(k,:)=nfr;
  end
  NFR(1,:)=NFR(2,:);
  NFR(nli,:)=NFR(nli-1,:);
 
  if ~exist('XL','var')
    nn=length(Xl);
    x1=Xl(1);
    y1=Yl(1);
    for ii=1:nn
      x2=Xl(ii);
      y2=Yl(ii);
      dXX(ii)=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
      x1=x2;
      y1=y2;
    end
    XL=cumsum(dXX);
% Correct distances: dXX is distance for "zigzag" sections
% that follow the straight lines, thus the total
% distancs is longer than the straight line distances
% need to adjust it:
%    xtot=7922; % that what it should be, km, from Myers
%    xmx=max(max(XL));
%    crct=xtot/xmx;
%    XL=XL*crct;
  end   

% Plot:
  Ing=find(NFR<=0);
  if ~isempty(NFR);
    fprintf('!!!!  Warning:  N freq<0 found, replaced with 1e-9 !!!\n');
    NFR(Ing)=1e-9;
  end
  
  lTr=log(sqrt(NFR)); % sec^-2 -> sec^-1

  ln1=LON(j1,i1);
  ln2=LON(j2,i1);
  cyc=YRPLT(1,3);
  stl=sprintf('%s, ln(sec-1), %4.4i/%2.2i/%2.2i, cyc %i, %s',...
	      fld0,DV(1:3),cyc,xname);
  nf=1;
  sub_plot_xsectZZ(nf,XL,ZZi,lTr,stl,fld0,xdr);

  ffg=sprintf('%ssgm_%s-%s_%3.3i',pthfig,xname,fld0,ip);

  if s_fig>0
    fprintf('Saving %s\n\n',ffg);
    print('-dpng','-r200',ffg);
  end
end

% --------------------
% Plot the cross-section
% delete unneeded topo:
HH(1:250,:)=nan;
HH(:,700:end)=nan;

fn=2;
sub_plot_bath(HH,LON,LAT,fn,domname);
contour(HH,[-3500:100:-10],'Color',[0.8 0.8 0.8]);
contour(HH,[-4000:500:-50],'Color',[0.4 0.4 0.4]);
plot(IIa,JJa,'k-','linewidth',2);
x=XL(1,:);
for kk=1:50:max(x)
  d=abs(x-kk);
  i0=find(d==min(d));
  plot(IIs(i0),JJs(i0),'ro','linewidth',2);
end

tt = ' ';
if f_tracloc>0
  ftrca = sprintf('%srelax_trcr_Greenland.a',PTH.data);
  ftrcb = sprintf('%srelax_trcr_Greenland.b',PTH.data);
  imonth = 7;
  FTr = sub_read_tracer_relax(ftrca,ftrcb,IND,imonth,nlev);
  lFTr = log(FTr); % kg/m3
  
  c1=-15;
  c2=5;
  nint=100;
  cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
  Rr=-1;                        % colors of max intensity, clockwise rotation 
  C0=[-1,0,1];                 % starting point, red dimension is Z axis
  Cend=[1, 1, 0.6];
  cmpT = colormap_spiral(nint,'C0',C0,'Rr',Rr,'Cend',Cend);
  
  pcolor(lFTr); shading flat;
  colormap(cmpT);

  hght=[];
  lngth=[];
  mint=10;
  mbx=mint;
  fsz=12;
  bxc='k';
  posc=[0.85 0.1 0.8 0.06];
  aend=0;
  colorbar_vert(cmpT,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);
  tt=sprintf(', Tracer Rlx, kg/m3, mo=%i',imonth);

end


axis('equal');
switch (xname)
 case('uummq')
  set(gca,'xlim',[170 300],'ylim',[600 720]); 
  stl = sprintf('Fjord Uummannaq %s',tt);
 case('scoresby');
  set(gca,'xlim',[440 580],'ylim',[500 630]); 
  stl = sprintf('Scoresby Sound%s',tt);
end  
set(gca,'xtick',[],'ytick',[]);
title(stl,'Fontsize',16);

bottom_text(txtb);

ffg=sprintf('%ssgm_%s-map',pthfig,xname);
if s_fig>0
  fprintf('Saving %s\n\n',ffg);
  print('-dpng','-r200',ffg);
end

