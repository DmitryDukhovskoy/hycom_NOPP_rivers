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
SEGM{5} = 'DenmarkS';
SEGM{6} = 'IclNorw';
SEGM{7} = 'FaroeShtl';
SEGM{8} = 'NorthEGr';
SEGM{9} = 'CntrEGr';
SEGM{10}= 'SouthEGr';
SEGM{11}= 'SouthWGr';
SEGM{12}= 'CntrWGr';
SEGM{13}= 'NarresS';
SEGM{14}= 'LancasterS';
nsgm = length(SEGM);


s_mat=1;

%Zmn = 100; % average over the top Zmn m
rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);
YRPLT=[];
cc=0;

%if expt=='061', yend=2008; end;
dd = 10;
for iyr=1993:2016
  nday=datenum(iyr,12,31)-datenum(iyr,1,1)+1;
  for iday=1:dd:nday
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
  end
end
%fmat = sprintf('%s%s_TS_straits_%4.4i_%4.4i.mat',...
%	       pthmat,expt,YRPLT(1,1),YRPLT(end,1));
fmat = sprintf('%s%s_TS_straits_%4.4i_%4.4iv2.mat',...
	       pthmat,expt,YRPLT(1,1),YRPLT(end,1));

if s_mat==0,
  fprintf('No mat files saved ... \n');
else
  fprintf('Data will be saved -> %s\n',fmat);
end

% Start from:
ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

[mm,nn]=size(LON);
SGMTS = struct;

SGMTS = sub_set_xsct_v2(HH,LON,LAT,SEGM);


nrec = 0;  % # of saved mean U,V averaged over Nav days
np   = length(YRPLT);
cnc  = 0;
ip1  = 1;
if s_mat == 2
  clear SGMTS
  fprintf('Loading saved fmat to restart %s\n',fmat);
  load(fmat);
  
  TM=SGMTS(1).TM;
  nrc=length(TM);
  dv0=datevec(TM(nrc));
  fprintf('===  Start from the last saved record %i/%2.2i/%2.2i\n',dv0(1:3));
  ip1=nrc+1;
  cnc=nrc;
  fprintf('# saved flds: %i\n',cnc);
  yr   = YRPLT(ip1,1);
  iday = YRPLT(ip1,2);
  j1d  = datenum(yr,1,1);
  dnmb = j1d+iday-1;
  DV   = datevec(dnmb);
  fprintf('===  Next record to process %i/%2.2i/%2.2i\n\n',DV(1:3));
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
    ig = 2;
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

  if s_mat>0 & mod(ip,50)==0
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




