% Vertical sections output archive files:
% for most recent sensit. experiments with CICE5
% old simulations 
%
% Fields are long-term average over specified time
% 
% Vertical transects of T and S
% Across the Arctic, Bering-to-North Pole-to-Nordic Seas along 

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
f_save = 2;  % =0 - no save; = 1 - save rewrite; =2 - continue saved

YR1 = 2017;
YR2 = YR1;
regn = 'ARCc0.04';
ires = 0.04;
expt = 23;
expt_nm = 'dltEddTmltEAPJRA';

%xname = 'BerNatl'; % Bering - N.Atl
%xname = 'BeaufNatl'; % Beaufort Shelf - South Norway
%xname = 'LaptSea'; % Laptev Sea
xname = 'BeaufIcel';  % Beaufort Sea to Iceland via N Pole
%xname = 'BafnFram'; % section around Greenland from Baffin to Fram str 


hg    = 2^100;
rg    = 9806;
% Change Tv and nlev for different ARCc fields
Tv    = 17; % bathym v11
nlev  = 41; % 41

pthout = '/nexsan/people/ddmitry/hycom/data_extracted/';

% ARCc0.04 010 similar to ARCc0.08 110 - tracers, river climatology, no Greenland
% ARCc0.04 012 similar to ARCc0.08-112 - tracers, river monthly data, Greenland on
% 022 - GOFS3.5 CICEv5

if f_save == 0
  fprintf(' !!!!!!!!!!!!!!!!! \n');
  fprintf(' Output not saved f_save=%i \n\n',f_save);
end

%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
% 
% 
YRPLT = [];
cc=0;
for iyr=YR1:YR2
  for idd=1:7:365
    if idd==1, idd=2; end;
    dJ1  = datenum(iyr,1,1);
    cc=cc+1;
    dnmb = dJ1+idd-1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
    YRPLT(cc,3)=dnmb;
  end
end

fprintf('Extracting %s %i-%i 0.04 HYCOM-CICE expt%3.3i\n',pfld,YR1,YR2,expt);

% Specify segments of the x-section
XSCT = sub_xsections(xname,ires);
IJs  = XSCT.IJs;
IIs  = IJs(:,1);
JJs  = IJs(:,2);

nS=length(IIs);
clear Xl Yl
for ii=1:nS
  i0=IJs(ii,1);
  j0=IJs(ii,2);
  Xl(ii,1)=LON(j0,i0);
  Yl(ii,1)=LAT(j0,i0);
  Hb(ii,1)=HH(j0,i0);
end;
INDs=sub2ind(size(HH),JJs,IIs);

btx = 'extr_vertTS_mean_arc004.m';
f_map=0;
if f_map>0
  figure(10); clf;
  set(gcf,'Position',[1000         483         752         833]);
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7]);
  plot(IIs,JJs,'b.-');
  axis('equal');
  title(xname);
%
% Plot distance marks along the section
% Distance from point 0
 nn=length(Xl);
 x1=Xl(1);
 y1=Yl(1);
 for ii=1:nn
  x2=Xl(ii);
  y2=Yl(ii);
  dXX(ii)=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
%  x1=x2;
%  y1=y2;
 end
%  Lsct = cumsum(dXX);
  Lsct = dXX;

  dL = 1000;
  for ll=0:dL:max(Lsct)
    dmm = abs(Lsct-ll);
    ii=find(dmm==min(dmm),1);
    i0=IIs(ii);
    j0=JJs(ii);

    plot(i0,j0,'r.','Markersize',12);
    lts=sprintf('%i',ll);
%    text(i0+5,j0,lts,'Fontsize',12);
  end

  bottom_text(btx,'pwd',1);
  keyboard
end

XSCT.resolution = ires;
XSCT.expt_nmb   = expt;
XSCT.expt_nm    = expt_nm;
XSCT.Sect_I     = IIs;
XSCT.Sect_J     = JJs;
XSCT.Sect_lon   = Xl;
XSCT.Sect_lat   = Yl; 
XSCT.Bott_depth = Hb;

fmat_out = sprintf('%shycom%3.3i_%3.3i_xsect_%s_%s_%i-%i.mat',...
                    pthout,ires*100,expt,xname,pfld,YR1,YR2);



np=size(YRPLT,1);
% Plot fields:
cnc=0;
ip1=1;
AAsum = [];
dlast = 0;

if f_save == 2
  fprintf('Start from saved %s\n',fmat_out);
  load(fmat_out);
  nrc = XSCT.Nmb_rcrds;
  AAsum = XSCT.Mean_Field * nrc;
  TM = XSCT.TM;
  dlast = XSCT.TM(nrc);

  fprintf('Last saved record #%i %s\n',nrc, datestr(dlast)); 
end

expt0 = expt;
for ip=ip1:np
  tic;
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);
  dnmbE = datenum(2017,03,31);


  if dnmb <= dnmbE
    pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%4.4i_BL99Tfrz/',yr);
    fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,22,yr,iday);
    finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,22,yr,iday);
  else
    pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%4.4i_%s/',...
                    expt,yr,expt_nm);
    fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
  end

  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;

  if dnmb <= dlast
    fprintf('%i %s <= last saved %s, skipping\n',cnc,datestr(dnmb),datestr(dlast));
    continue;
  end

  fprintf(':: %4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);
  
  [F,n,nlev] = read_hycom(fina,finb,pfld);
  F(F>1e6)=nan;

  T=squeeze(F(:,INDs));
  [a1,a2]=size(T);

%
% Prepare vertical thicknesses 
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e10)=nan;
  F(F<0.1)=nan;
  F=F/rg;
  Dsec=squeeze(F(:,INDs));
  Dsec(Dsec==0)=1.0e-16;  % for interpolation, avoid nans
% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
  clear ZZb
  Dsec(isnan(Dsec))=0;
  ZZb(1,:)=-Dsec(1,:);
  for kk=2:l
    ZZb(kk,:)=ZZb(kk-1,:)-Dsec(kk,:);
  end

  [nl,npb]=size(ZZb);
  ZZ=zeros(nl+1,npb);
  ZZ(2:nl+1,:)=ZZb;
 
% For interpolation - add surface =0 values:
%  AA = zeros(nl+1,npb);
%  AA(1,:) = T(1,:);
%  AA(2:nl+1,:) = T;
  AA = T;
 
% Depths of middle of the cells:
  ZM(1,:)=0.5*ZZ(1,:);
  for kk=1:l
    ZM(kk,:)=0.5*(ZZ(kk+1,:)+ZZ(kk,:));
  end

% Interpolate into reg. vert. grid:
  fprintf('Interpolating into Z-levels ...\n');
  npb = size(AA,2);
  ZZi=[0:-1:-100,-105:-5:-500, -510:-10:-2000, -2020:-20:-6000]';
  nli=length(ZZi);
  AAi=zeros(nli,npb)*nan;
  for ipp=1:npb
    Ib=min(find(Dsec(:,ipp)<=1.e-9));
    if Ib==1, continue; end;
    tt=AA(:,ipp);
    zz=ZM(:,ipp);
    hb=-sum(Dsec(:,ipp));
    ibz = max(find(ZZi>=hb));
    for kl=Ib:nl
      zz(kl)=zz(kl-1)-0.1;
    end;
    zz=[0;zz];
    tt=[tt(1);tt];
    if zz(nl)>ZZi(end)
      zz(nl)=ZZi(end);
    end
    aai = interp1(zz,tt,ZZi,'pchip');
    aai(ibz+1:end)=nan;
    AAi(:,ipp)=aai;
  end

  Nintf=size(ZZ,1);
  if ~exist('XL','var')
    nn=length(Xl);
    x1=Xl(1);
    y1=Yl(1);
    for ii=1:nn
      x2=Xl(ii);
      y2=Yl(ii);
% Distance from point 0:
      dxx=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
      XL(ii) = dxx;
%      dxx=max([dxx,0.01]);
%      dXX(ii)=dxx;
%      x1=x2;
%      y1=y2;
%      XL(ii)=sum(dXX(1:ii));
    end
% Note correct distances: dXX is distance for "zigzag" sections
    XSCT.Depths    = ZZi;
    XSCT.Sect_dist = XL;
  end   
  
%keyboard
 
  f_plt = 0;
  if f_plt==1 
    stl=sprintf('%4.2f, %s %s, %4.4i/%2.2i/%2.2i, 170W-10E ',...
         ires,expt_nm,pfld,DV(1:3));
    nf=1;
    xdr = 1;
    f_layer=0;
    [ma,na]=size(AA);
    s0=[];
    sub_plot_xsectZ(nf,XL,ZZi,AAi,stl,fld0,xdr,Hb);
    bottom_text(btx,'pwd',1);

    keyboard
  end

% Sum fields  
  if cnc == 1
    AAsum = AAi;
  else
    AAsum = AAsum + AAi;
  end

  TM(cnc,1) = dnmb; 
  AA = AAsum / cnc;
  XSCT.Mean_Field = AA;
  XSCT.Nmb_rcrds = cnc;
  XSCT.TM = TM;
  if mod(cnc,10) == 0 & f_save>0
    fprintf('Saving output ---> %s',fmat_out);
    save(fmat_out,'XSCT');
  end

  fprintf('Processed 1 record %6.2f min\n',toc/60);        
end

if f_save>0
  fprintf('Saving output ---> %s\n',fmat_out);
  save(fmat_out,'XSCT');
end





 
