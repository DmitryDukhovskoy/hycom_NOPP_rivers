% Calculate volume transport across straits/sections
% from ARCc0.04 simulations
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.04';
%expt = 011;
expt = 012;
TV   = 17;  % topo version
%segm = 'BeringS';
%segm = 'FramS';
%segm = 'BarOp';
%segm = 'DavisS';

SEGM{1} = 'BeringS';
SEGM{2} = 'FramS';
SEGM{3} = 'BarOp';
SEGM{4} = 'DavisS';
nsgm = length(SEGM);


s_mat=1;

%Zmn = 100; % average over the top Zmn m
rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.04/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);
YRPLT=[];
cc=0;

%if expt=='061', yend=2008; end;
dd = 5;
for iyr=2005:2005
  nday=datenum(iyr,12,31)-datenum(iyr,1,1)+1;
%  for iday=243:dd:nday
  for iday=243:dd:247
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
  end
end
fmat = sprintf('%s%3.3i_TS_straits_%4.4i_%4.4i.mat',...
	       pthmat,expt,YRPLT(1,1),YRPLT(end,1));

if s_mat==0,
  fprintf('No mat files saved ... \n');
else
  fprintf('Data will be saved -> %s\n',fmat);
end

% Start from:
%ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

[mm,nn]=size(LON);
SGMTS = struct;
for ig=1:nsgm
  segm=SEGM{ig};
  switch(segm)
   case('BeringS');
    is1=634*2;
    js1=1919*2;
    is2=658*2;
    js2=js1;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('FramS');
%    is1=935; too far north - ~80N
%    js1=959;
%    is2=1070;
%    js2=js1;
    is1=908*2; % this ~78.5
    is2=1106*2;
    js1=915*2;
    js2=915*2;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('BarOp');
    is1=1166*2;
    js1=456*2;
    is2=is1;
    js2=929*2;
    X=LON(js1:js2,is1);
    Y=LAT(js1:js2,is1);
    J=[js1:js2]';
    I=ones(size(J))*is1;
   case('DavisS');
    is1=478*2;
    js1=680*2;
    is2=568*2;
    js2=680*2;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
  end

  p_map=0;
  if p_map>0
    figure(11); clf;
    contour(HH,[0 0],'k');
    hold on
    plot(I,J,'r.-');
    keyboard
  end



  SGMTS(ig).Name = segm;
  SGMTS(ig).X    = X;
  SGMTS(ig).Y    = Y;
  SGMTS(ig).I    = I;
  SGMTS(ig).J    = J;
  ll=length(I);
  clear D
  for k=1:ll-1
    x1=X(k);
    y1=Y(k);
    x2=X(k+1);
    y2=Y(k+1);
    D(k,1)=distance_spheric_coord(y1,x1,y2,x2);
  end
  D(k+1,1)=D(k);
  SGMTS(ig).Dist = D;

end


nrec=0;  % # of saved mean U,V averaged over Nav days
np=size(YRPLT,1);
cnc=0;
for ip=1:np
  yr   = YRPLT(ip,1);
  iday = YRPLT(ip,2);
  j1d  = datenum(yr,1,1);
  dnmb = j1d+iday-1;
  DV   = datevec(dnmb);
  
  if expt==11;
    pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data/%i/',yr);  
  else
    pthbin = sprintf('/nexsan/hycom/ARCc0.04_012/data/%i/',yr);  
  end
  
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  if ~exist(fina,'file') | ~exist(finb,'file');
    fprintf('Not found *a or *b: %s\n',fina);
    fprintf('                     %s\n',finb);
    FLXV(cnc,1)=nan;
    continue;
  end
    
  cnc=cnc+1;
  TM(cnc,1)=dnmb;
  
  fprintf('Getting data expt %3.3i: %4.4i_%2.2i_%2.2i: %s\n',expt,DV(1:3),fina);


  tic;
  [F,n,m,nlr] = read_hycom(fina,finb,'temp');
%  F=F(:,jnc1:jnc2,inc1:inc2);
  F(F>hgg)=nan;
  T=F;

  [F,n,m,nlr] = read_hycom(fina,finb,'salin');
  F(F>hgg)=nan;
  S=F;

  [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
  ZZ(isnan(ZZ))=100;
  ZM(isnan(ZM))=100;
  

% Process segments
  for ig = 1:nsgm
    nm = SGMTS(ig).Name;
    I  = SGMTS(ig).I;
    J  = SGMTS(ig).J;
    DX = SGMTS(ig).Dist;

    [b1,b2]=size(I);
    if b1<b2, I=I'; end;
    [b1,b2]=size(J);
    if b1<b2, J=J'; end;

    xsct = 0;
    if J(end)==J(1), % X-section
      if I(end)<I(1)
	I=flipud(I);
	J=flipud(J);
      end
      xsct = 1; 
    end
    if I(end)==I(1) % Y-section
      if J(end)<J(1)
	I=flipud(I);
	J=flipud(J);
      end
    end

    if xsct==1
      i1=min(I);
      i2=max(I);
      j1=J(1);
      tsgm = squeeze(T(:,j1,i1:i2));
      ssgm = squeeze(S(:,j1,i1:i2));
      zz   = squeeze(ZZ(:,j1,i1:i2));
      zm   = squeeze(ZM(:,j1,i1:i2));
    else
      j1=min(J);
      j2=max(J);
      i1=I(1);
      tsgm = squeeze(T(:,j1:j2,i1));
      ssgm = squeeze(S(:,j1:j2,i1));
      zz   = squeeze(ZZ(:,j1:j2,i1));
      zm   = squeeze(ZM(:,j1:j2,i1));
    end
    
    SGMTS(ig).TM=TM;
    SGMTS(ig).Temp(cnc,:,:) = tsgm;
    SGMTS(ig).Saln(cnc,:,:) = ssgm;
    SGMTS(ig).ZZ(cnc,:,:) = zz;
    SGMTS(ig).ZM(cnc,:,:) = zm;
    
  end;  % segments loop

    
    
    
%  keyboard
  plt_sec=0;
  if plt_sec>0
    ig = 1;
    it = 1;
    tt = squeeze(SGMTS(ig).Temp(it,:,:));
    ss = squeeze(SGMTS(ig).Saln(it,:,:));
    zz = squeeze(SGMTS(ig).ZZ(it,:,:));
    zm = squeeze(SGMTS(ig).ZM(it,:,:));
    dd = SGMTS(ig).Dist;
    snm= SGMTS(ig).Name;
    zz(zz>0)=0;
%    zz(isnan(tt))=nan;
    
    tt(end+1,:)=tt(end);
    ss(end+1,:)=ss(end);
    xx = (cumsum(dd)-dd(1))*1e-3;
    
    figure(10); clf;
%    pcolor(xx,zz,tt); shading flat;
    pcolor(xx,zz,ss); shading flat;
    set(gca,'xlim',[xx(1) xx(end)]);
    stl = sprintf('%s',snm);
    title(stl);

    ZZ=cumsum(DZ,1);
    ZZ(II)=nan;
    XX=cumsum(D)*1e-3;
    pcolor(XX,ZZ,FLX); shading flat
    caxis([-12000 12000]);
    colorbar
    stt=sprintf('%s, Volume Flux, m3/s',nm);
    title(stt);
  end

  fprintf('1 day %6.3f min\n',toc/60);
  fprintf('== %i/%2.2i/%2.2i,  %s Tmin: %6.2f Tmax=%6.2f\n\n',...
	  DV(1:3),nm,min(min(tsgm)),max(max(tsgm)));

  if s_mat>0 & mod(ip,10)==0
    fprintf('####  Saving %s\n\n',fmat);
    save(fmat,'SGMTS');
  end


end

if s_mat>0
  fprintf('####  Saving %s\n\n',fmat);
  save(fmat,'SGMTS');
end

%save(fout,'UCPR');

%exit




