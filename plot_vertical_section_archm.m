% Vertical sections output archive files:
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;
%pfld  = 'temp';
pfld  = 'salin';
fld0  = pfld;
%xname = 'BerNatl'; % Bering - N.Atl
%xname = 'BeaufNatl'; % Beauf. Shelf - South Norway
xname = 'LaptSea'; % Laptev Sea

hg    = 2^100;
rg    = 9806;
% Change Tv and nlev for different ARCc fields
%Tv    = 7; % bathym version: 07 or 11
%nlev  = 32; % 32 or 41
Tv    = 11; % bathym version: 07 or 11
nlev  = 41; % 32 or 41

% ARCc0.08 110 - tracers, river climatology, no Greenland
% ARCc0.08-112 - tracers, river monthly data, Greenland on
expt = 112;  % ARCc HYCOM experiment #

regn = 'ARCc0.08';
%pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_xsect/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,Tv); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
% 
%domname = 'NorthAtl';
%IND = smaller_domain_indices('Green');
%IND = smaller_domain_indices(domname);
%IND = smaller_domain_indices('Baffin');

inc1=1;
inc2=nn;
jnc1=1;
jnc2=mm;
djnc=mm;
dinc=nn;


% 
YRPLT = [];
cc=0;
for iyr=2005:2005
  for idd=210:210
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

% Specify segments of the x-section
switch(xname);
 case('BerNatl');
  IJs=[645, 1916;
       934, 1242;
       1038, 855;
       1020, 495;
       815, 179];
 case('BeaufNatl');
  IJs=[479, 1687;
       928, 1250
       1144, 272]
 case('LaptSea');
  IJs=[1096 1446
       1170 1546
       1170 1784];
  xl1=650;
  xl2=1420;
  yl1=-100;
  yl2=0;
  cs1=10;
  cs2=34;  
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

f_map=1;
if f_map>0
  figure(10); clf;
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7]);
  plot(IIs,JJs,'b.-');
%  set(gca,'xlim',[100 600],...
%	  'ylim',[200 1000]);
end

fprintf('Section: %s, Flag fig: %i\n',xname,s_fig);

%keyboard


np=size(YRPLT,1);
% Plot fields:
cnc=0;
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
%  pthbin = sprintf(...
%    '/nexsan/people/ddmitry/hycom/ARCc0.08/%s/bin_outp/',expt);
%  pthbin = sprintf('/nexsan/archive/%s_%3.3i/1997/',regn,expt)
  pthbin = sprintf('/nexsan/hycom/ARCc0.08_%3.3i/data/%i/',expt,yr);

  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

%fprintf('Plotting: Tracer #%i V.Layer=%i\n\n',nTr,plr);
  fprintf(':: %4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);
% in archm u_vel=utot; total velocity
%  keyboard
  tic;
  
  [F,n,nlev] = read_hycom(fina,finb,pfld);
  toc;
  F=F(:,jnc1:jnc2,inc1:inc2);
  F(F>1e6)=nan;

  T=squeeze(F(:,INDs));
  [a1,a2]=size(T);
% Prepare for plotting - add extra bogus layer
% at the bottom, Matlab won't show it
  T(end+1,1:a2)=T(end,:);
  AA = T; 

%  [F,n,nlev] = read_hycom(fina,finb,'salin');
%  F=F(:,jnc1:jnc2,inc1:inc2);
%  F(F>1e6)=nan;
%
%  S=squeeze(F(:,INDs));

%  Rho = sw_dens0(S,T);
  
  
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
  f_intrp2zz=0;
  if f_intrp2zz==1
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
  end  

  Nintf=size(ZZ,1);
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
      XL(1:Nintf,ii)=sum(dXX(1:ii));
    end
% Note correct distances: dXX is distance for "zigzag" sections
  end   
  
  
  ln1=LON(j1,i1);
  ln2=LON(j2,i1);
  stl=sprintf('%s-%3.3i, %s, %4.4i/%2.2i/%2.2i, %s',...
	      regn,expt,fld0,DV(1:3),xname);
  nf=1;
  xdr = 1;
  f_layer=0;
  [ma,na]=size(AA);
  if strncmp(pfld,'salin',4)
    c1=cs1;
    c2=cs2;
  else
    c1=ct1;
    c2=ct2;
  end
  
%  sub_plot_xsection(nf,XL,ZZ,AA,stl,fld0,f_layer);
  sub_plot_xsection2(nf,XL,ZZ,AA,stl,fld0,f_layer,...
		     xl1,xl2,yl1,yl2,c1,c2);
  if f_intrp2zz  
    sub_plot_xsectZZ(nf,XL,ZZi,lTr,stl,fld0,xdr);
  end

  ffg=sprintf('%s%s-%s_%2.2i%2.2i%4.4i',...
	      pthfig,xname,fld0,DV(3),DV(2),DV(1));
  txtb = 'plot_vertical_section_archm.m';
  bottom_text(txtb,'pwd',1,'Fontsize',10);
  if s_fig>0
    fprintf('Saving %s\n\n',ffg);
    print('-dpng','-r200',ffg);
  end
end





 
