% Calculate volume transport across straits/sections
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt = '110';
TV   = 11;  % topo version
%segm = 'BeringS';
%segm = 'FramS';
%segm = 'BarOp';
%segm = 'DavisS';

SEGM{1} = 'BeringS';
SEGM{2} = 'FramS';
SEGM{3} = 'BarOp';
SEGM{4} = 'DavisS';
nsgm = length(SEGM);

s_mat=2; % =2 - start from last saved record

rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);
YRPLT=[];
cc=0;

%if expt=='061', yend=2008; end;
dd = 5;
for iyr=1993:2016
  nday=datenum(iyr,12,31)-datenum(iyr,1,1)+1;
  for iday=1:dd:nday
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
  end
end
fmat = sprintf('%s%s_UV_straits_%4.4i_%4.4i.mat',...
	       pthmat,expt,YRPLT(1,1),YRPLT(end,1));

if s_mat==0,
  fprintf('No mat files saved ... \n');
else
  fprintf('Data will be saved -> %s\n',fmat);
end


ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

[mm,nn]=size(LON);
UV = struct;

for ig=1:nsgm
  segm=SEGM{ig};
  switch(segm)
   case('BeringS');
    is1=634;
    js1=1919;
    is2=658;
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
    is1=908; % this ~78.5
    is2=1106;
    js1=915;
    js2=915;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('BarOp');
    is1=1166;
    js1=456;
    is2=is1;
    js2=929;
    X=LON(js1:js2,is1);
    Y=LAT(js1:js2,is1);
    J=[js1:js2]';
    I=ones(size(J))*is1;
   case('DavisS');
    is1=478;
    js1=680;
    is2=568;
    js2=680;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
  end
  
  p_map=0;
  if p_map>0
%    ig = 2;
%    I = UV(ig).I;
%    J = UV(ig).J;
    figure(11); clf;
    contour(HH,[0 0],'k');
    hold on
    plot(I,J,'r.-');
    keyboard
  end

  UV(ig).Name = segm;
  UV(ig).X    = X;
  UV(ig).Y    = Y;
  UV(ig).I    = I;
  UV(ig).J    = J;
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
  UV(ig).Dist = D;

end


np=length(YRPLT);
cnc=0;
ip1=1;

if s_mat == 2
  clear UV
  fprintf('Loading saved fmat to restart %s\n',fmat);
  load(fmat);
  
  TM=UV(1).TM;
  nrc=length(TM);
  dv0=datevec(TM(nrc));
  fprintf('Start from the last saved record %i/%2.2i/%2.2i\n',dv0(1:3));
  ip1=nrc+1;
  cnc=nrc;
  fprintf('# saved flds: %i\n',cnc);
  yr   = YRPLT(ip1,1);
  iday = YRPLT(ip1,2);
  j1d  = datenum(yr,1,1);
  dnmb = j1d+iday-1;
  DV   = datevec(dnmb);
  fprintf('Next record to process %i/%2.2i/%2.2i\n',DV(1:3));
%  keyboard
end  


for ip=ip1:np
  yr   = YRPLT(ip,1);
  iday = YRPLT(ip,2);
  j1d  = datenum(yr,1,1);
  dnmb = j1d+iday-1;
  DV   = datevec(dnmb);
  
  pthbin = sprintf('/nexsan/archive/ARCc0.08_%s/data/%4.4i/',expt,yr);
  
  fina = sprintf('%s%s_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%s_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  if ~exist(fina,'file') | ~exist(finb,'file');
    fprintf('Not found *a or *b: %s\n',fina);
    fprintf('                     %s\n',finb);
    FLXV(cnc,1)=nan;
    continue;
  end
    
  cnc=cnc+1;
  TM(cnc,1)=dnmb;
  
  fprintf('Getting data expt %s: %4.4i_%2.2i_%2.2i: %s\n',expt,DV(1:3),fina);


  
% Calculate transport:
% in archm, archv: u, v anre not collocated 
% with H, dH:
% u(i,j) is between h(i-1,j) and h(i,j)
% To check look at near land points
% such that HH(j-1,i) and HH(j,i-1) is land and HH(j,i) is not
% if V(j,i) is nan  than V is not collocated with H
% if U(j,i) is nan than U is not collocated with H
%need barotropic U?

  tic;
  [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
%  F=F(:,jnc1:jnc2,inc1:inc2);
  U1=squeeze(F(1,:,:));
  F(F>hgg)=0;
  UN=F;

  [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
  V1=squeeze(F(1,:,:));
  F(F>hgg)=0;
  VN=F;

  if ~exist('colloc','var'),
    colloc = check_collocatedU(HH,U1,V1);
    if colloc
      fprintf('U and V are collocated with H points\n');
      fprintf('The code is written for not collocated U,V\n');
      error(' Use different code for collocated U/V ...');
    end
    fprintf('         U, V Not collocated \n');
  end

  [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
  ZZ(isnan(ZZ))=100;
  ZM(isnan(ZM))=100;
  

% Process segments
  for ig = 1:nsgm
    nm = UV(ig).Name;
    I  = UV(ig).I;
    J  = UV(ig).J;
    DX = UV(ig).Dist;

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
      UVsgm = squeeze(VN(:,j1,i1:i2)); % save normal component
      zz   = squeeze(ZZ(:,j1,i1:i2));
    else
      j1=min(J);
      j2=max(J);
      i1=I(1);
      UVsgm = squeeze(UN(:,j1:j2,i1));
      zz   = squeeze(ZZ(:,j1:j2,i1));
    end
    
    UV(ig).TM = TM;
    UV(ig).UV_normal(cnc,:,:) = UVsgm;
    UV(ig).ZZ(cnc,:,:) = zz;

    umx = max(max(UVsgm));
    umn = min(min(UVsgm));
    fprintf(':: %s Umax=%4.2f  Umin=%4.2f m/s\n',nm,umx,umn);
    
  end;  % segments loop
    
  
%keyboard
  plt_sec=0;
  if plt_sec>0
    ig = 3;
    it = 1;
    vv = squeeze(UV(ig).UV_normal(it,:,:));
    zz = squeeze(UV(ig).ZZ(it,:,:));
    dd = UV(ig).Dist;
    snm= UV(ig).Name;
    zz(zz>0)=0;
    
    vv(end+1,:)=vv(end);
    xx = (cumsum(dd)-dd(1))*1e-3;
    
    figure(10); clf;
    pcolor(xx,zz,vv); shading flat
    caxis([-0.2 0.2]);
    colorbar
    stt=sprintf('%s, V normal',snm);
    title(stt);
  end
  
%  FLXZ=nansum(FLX);

  fprintf(' ===  %i/%2.2i/%2.2i, 1 day %6.3f min\n\n',DV(1:3),toc/60);

  if s_mat>0 & mod(ip,50)==0
    fprintf('####  Saving %s\n\n',fmat);
    save(fmat,'UV');
  end


end

if s_mat>0
  fprintf('####  Saving %s\n\n',fmat);
  save(fmat,'UV');
end

%save(fout,'UCPR');

%exit




