% Plot vertical disctribution
% of T & S  in fjords
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

pfld = 'salin';
s_mat = 0; % =0 - extract data, no save; =1 extract & save, =2 - load saved
s_fig = 0;
s_map = 0; % plot map with segments
f_temp =0; % plot T section with S

%f_tracloc = 0; % =1 - add locations of tracer release to transect maps
rg=9806;  % convert pressure to depth, m
%fld0='salin'; %'temp' - background field

plr=0; % highlight this interface
btx = 'xsct_ts_frames_fjords008.m';

regn = 'ARCc0.08';
%expt = 110; % experiment without runoff
expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthfig=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_green_xsctFR/',...
		  regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


YRPLT=[];
cc=0;
for iyr=2015:2015
%  for idd=226:226
  for idd=45:45
%  for idd=165:30:365
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np   = size(YRPLT,1);
SGM  = sub_fjord_sections(HH,LON,LAT);
%A = SGM(3); % do only 2 sections
%A(2)=SGM(4);
%SGM=A; 
nsgm = length(SGM);

f_map_check=0; % quickly see the segments
if f_map_check>0
  figure(10); clf;
%  fn = 10;
%  sub_plot_bath(HH,LON,LAT,fn,domname);
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-1000:100:-100],'Color',[0.8 0.8 0.8]);
  contour(HH,[-5000:1000:-1000],'Color',[0.9 0.9 0.9]);
  for kk=1:nsgm
    IIs = SGM(kk).IIs;
    JJs = SGM(kk).JJs;
    plot(IIs,JJs,'b.-');
    stx=sprintf('Sect %2.2i %s',kk,SGM(kk).Name);
    text(min(IIs),min(JJs),stx);
  end
  axis('equal');
  set(gca,'xlim',[500 960],...
	  'ylim',[400 1080]);
  bottom_text(btx,'pwd',1);
end

%fprintf('Section: %s, Saving fig: %i\n',xname,s_fig);
%keyboard

% Plot fields:
cnc=0;  % 
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)
  
  if expt==112
    pthbin = sprintf('/nexsan/hycom/ARCc0.08_112/data/%i/',yr);
  end
    
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

  fprintf('%4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

  tic;
  [F,n,m,nlev] = read_hycom(fina,finb,'temp');
  F(F>1e6)=nan;
  TT = F;

  [F,n,m,nlev] = read_hycom(fina,finb,'salin');
  F(F>1e6)=nan;
  SS = F;

  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e18)=nan;
  F=F/rg;
  F(F<1e-2)=0;
  dH = F;
  
  for kk=8:8
%  for kk=1:nsgm
    INDs = SGM(kk).INDs;
    xnm  = SGM(kk).Name;
    T=squeeze(TT(:,INDs));
    S=squeeze(SS(:,INDs));
%    Rho = sw_dens0(S,T);
  
% Prepare for plotting - add extra bogus layer
% at the bottom, Matlab won't show it
    [a1,a2]=size(T);
  
% Prepare vertical thicknesses for plotting  
    Dsec=squeeze(dH(:,INDs));
    Dsec(Dsec==0)=nan;
% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
    clear ZZb
    Dsec(isnan(Dsec))=0;
    ZZb(1,:)=-Dsec(1,:);
    for kl=2:l
      ZZb(kl,:)=ZZb(kl-1,:)-Dsec(kl,:);
    end
% For plotting need to have #of layers +1 interfaces
% add surface layer, otherwise all values
% will be shifted 1 layer down
    [nl,npb]=size(ZZb);
    ZZ=zeros(nl+1,npb);
    ZZ(2:nl+1,:)=ZZb;
    ZZav=ZZ;
  
% Depths of middle of the cells:
    clear ZM
    ZM(1,:)=0.5*ZZ(1,:);
    for kl=1:l
      ZM(kl,:)=0.5*(ZZ(kl+1,:)+ZZ(kl,:));
    end

% Interpolate into reg. vert. grid:
%    hb=min(min(ZZb));
    ZZi=[0:-5:-1800]';
    nli=length(ZZi);
    Ti=zeros(nli,npb)*nan;
    Si=zeros(nli,npb)*nan;
    Hs=HH(INDs);
    for ipp=1:npb
      Ib=min(find(Dsec(:,ipp)==0));
      if Ib==1, continue; end;
      t=T(:,ipp);
      s=S(:,ipp);
      zz=ZM(:,ipp);
      hb=-sum(Dsec(:,ipp));
      ibz = max(find(ZZi>=hb));
      for kl=Ib:nl
	zz(kl)=zz(kl-1)-0.1;
      end;
      zz=[0;zz];
      t=[t(1);t];
      s=[s(1);s];
      if zz(nl)>ZZi(end)
	zz(nl)=ZZi(end);
      end
      si = interp1(zz,s,ZZi,'pchip');
%      si(ibz+1:end) = nan; % nan for bottom
      Si(:,ipp) = si;
      ti = interp1(zz,t,ZZi,'pchip');
%      ti(ibz+1:end) = nan; % na for bottom
      Ti(:,ipp) = ti;
    end

    XL = SGM(kk).dist_m*1e-3; % distance along segment, km  
    if strncmp(xnm,'WGr',3),
      xdr=-1;
    else
      xdr=1;
    end

% Plot:
    fld0='salin';
    stl=sprintf('%i, sct %2.2i, %s: %3.1f-%3.1f, %4.4i/%2.2i/%2.2i',...
	      expt, kk,fld0,min(min(Si)),max(max(Si)),DV(1:3));
    nf1=kk;
    sub_plot_xsectZZ(nf1,XL,ZZi,Si,stl,fld0,xdr,Hs);
      bottom_text(btx,'pwd',1,'position',[0.01 0.2 0.5 0.05]);
    set(gcf,'Position',[144 173 1975 1123]);

%    keyboard
    if s_fig>0
%      bottom_text(btx,'pwd',1,'position',[0.01 0.2 0.5 0.05]);
      ffg=sprintf('%s%3.3i_sgm%2.2i_%s_%4.4i',...
		pthfig,expt,kk,fld0,cnc);
      fprintf('Saving %s\n\n',ffg);
      print('-dpng','-r200',ffg);
    end
   
    if f_temp==1
    fld0='temp';
    stl=sprintf('%i, sct %2.2i, %s: %3.1f-%3.1f, %4.4i/%2.2i/%2.2i',...
	      expt, kk,fld0,min(min(Ti)),max(max(Ti)),DV(1:3));
    nf=nf1+20;
    sub_plot_xsectZZ(nf,XL,ZZi,Ti,stl,fld0,xdr,Hs);
      bottom_text(btx,'pwd',1,'position',[0.01 0.2 0.5 0.05]);
    set(gcf,'Position',[144 173 1975 1123]);
    
%    keyboard
    if s_fig>0
%      bottom_text(btx,'pwd',1,'position',[0.01 0.2 0.5 0.05]);
      ffg=sprintf('%s%3.3i_sgm%2.2i_%s_%4.4i',...
		pthfig,expt,kk,fld0,cnc);
      fprintf('Saving %s\n\n',ffg);
      print('-dpng','-r200',ffg);
    end
    end

  end % segments

  fprintf('1 time processed %6.1fmin\n\n',toc/60);
  
end



