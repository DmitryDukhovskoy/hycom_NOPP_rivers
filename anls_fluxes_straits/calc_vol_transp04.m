% Calculate volume transport across straits/sections
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear


YR1=2017;
YR2=YR1;
dday=5;

%segm = 'BeringS';
%segm = 'FramS';
%segm = 'BarOp';
%segm = 'DavisS';
segm = 'DenmarkS';

% This code is for archv output fields
% CHeck if archm - do not add U/V barotropic !!! - use extr_*mean_ 
s_mat=2;  % =2 - load saved and start from the last record

ixx    = 9; % experiment name and dir - check with EXPT - expt 023
%ixx    = 6;  % expt 022 original 
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;


%Zmn = 100; % average over the top Zmn m
rg  = 9806;
hgg = 1e20; % 


pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/strait_fluxes/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

YRPLT=[];
cc=0;

%if expt=='061', yend=2008; end;
dd = 5;
for iyr=YR1:YR2
  nday=datenum(iyr,12,31)-datenum(iyr,1,1)+1;
  for iday=2:dd:nday
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
    j1d  = datenum(iyr,1,1);
    dnmb = j1d+iday-1;
    YRPLT(cc,3)=dnmb;
  end
end
fmat = sprintf('%s%3.3i_Vflux_%s_%4.4i_%4.4i.mat',...
	       pthmat,expt,segm,YRPLT(1,1),YRPLT(end,1));

if s_mat==0,
  fprintf('No mat files saved ... \n');
else
  fprintf('Data will be saved -> %s\n',fmat);
end

% Start from:
ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


switch(segm)
 case('BeringS');
  is1=634;
  js1=1919;
  is2=658;
  js2=js1;
  is1=2*is1;
  is2=2*is2;
  js1=2*js1;
  js2=2*js2;
  X=LON(js1,is1:is2);
  Y=LAT(js1,is1:is2);
  I=[is1:is2]';
  J=ones(size(I))*js1;
 case('FramS');
  is1=935;
  js1=959;
  is2=1070;
  js2=js1;
  is1=2*is1;
  is2=2*is2;
  js1=2*js1;
  js2=2*js2;
  X=LON(js1,is1:is2);
  Y=LAT(js1,is1:is2);
  I=[is1:is2]';
  J=ones(size(I))*js1;
 case('BarOp');
  is1=1166;
  js1=456;
  is2=is1;
  js2=929;
  is1=2*is1;
  is2=2*is2;
  js1=2*js1;
  js2=2*js2;
  X=LON(js1:js2,is1);
  Y=LAT(js1:js2,is1);
  J=[js1:js2]';
  I=ones(size(J))*is1;
 case('DavisS');
  is1=478;
  js1=680;
  is2=568;
  js2=680;
  is1=2*is1;
  is2=2*is2;
  js1=2*js1;
  js2=2*js2;
  X=LON(js1,is1:is2);
  Y=LAT(js1,is1:is2);
  I=[is1:is2]';
  J=ones(size(I))*js1;
 case('DenmarkS');
  is1=1413; % already on 0.04 grid !
  js1=1103;
  is2=1714;
  js2=1103;
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



SEGM.Name=segm;
SEGM.X=X;
SEGM.Y=Y;
SEGM.I=I;
SEGM.J=J;
ll=length(I);
for k=1:ll-1
  x1=X(k);
  y1=Y(k);
  x2=X(k+1);
  y2=Y(k+1);
  D(k,1)=distance_spheric_coord(y1,x1,y2,x2);
end
SEGM.Dist=D;


nrec=0;  % # of saved mean U,V averaged over Nav days
np=length(YRPLT);
cnc=0;
ip1=1;

if s_mat==2
  load(fmat);
  cnc = length(TM);
  dlst = TM(cnc);
  ip1 = find(YRPLT(:,3)==dlst)+1;
  if isempty(ip1),
% should not happen unless YRPLT is changed during the code exectution
    fprintf('Couldnt find exact date in YRPLT %s\n',datestr(dlst)); 
    ip1=max(find(YRPLT(:,3)<=dlst))+1;
    fprintf('Use next closest date to start: %s\n\n',datestr(YRPLT(ip1,3))); 
  end
  
end


for ip=ip1:np
  yr   = YRPLT(ip,1);
  iday = YRPLT(ip,2);
  dnmb = YRPLT(ip,3);
  DV   = datevec(dnmb);
 
	pthbin=sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%4.4i_%s/',...
									expt,yr,texpt);
	fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
	finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
 

  if yr==2017 & DV(2)<4
    pthbin='/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_BL99Tfrz/';
    fina = sprintf('%s022_archv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
    finb = sprintf('%s022_archv.%4.4i_%3.3i_00.b',pthbin,yr,iday);
  end
    

  cnc=cnc+1;
  TM(cnc,1)=dnmb;

  if ~exist(fina,'file') | ~exist(finb,'file');
    fprintf('Not found *a or *b: %s\n',fina);
    fprintf('                     %s\n',finb);
    FLXV(cnc,1)=nan;
    continue;
  end
    
  
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

	[F,n,m,llr] = read_hycom(fina,finb,'u_btrop');
	F(F>hgg)=0;
	Ub=squeeze(F);

% Add barotropic U,V for archv files:
	for ilr=1:nlr
		UN(ilr,:,:)=squeeze(UN(ilr,:,:))+Ub;
	end

  [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
  V1=squeeze(F(1,:,:));
  F(F>hgg)=0;
  VN=F;

	[F,n,m,llr] = read_hycom(fina,finb,'v_btrop');
	F(F>hgg)=0;
	Vb=squeeze(F);

	for ilr=1:nlr
		VN(ilr,:,:)=squeeze(VN(ilr,:,:))+Vb;
	end

  if ~exist('colloc','var'),
    colloc = check_collocatedU(HH,U1,V1);
    if colloc
      fprintf('U and V are collocated with H points\n');
      fprintf('The code is written for not collocated U,V\n');
      error(' Use different code for collocated U/V ...');
    end
    fprintf('         U, V Not collocated \n');
  end

  [F,n,m,nlr] = read_hycom(fina,finb,'thknss');
  F=F./rg;
  F(F>hgg)=0;
  dH=F; 

  UN(dH<0.001)=nan;
  VN(dH<0.001)=nan;

% Calculate transport
% Sign convention:
% U+ - positive X (right)
% V+ - positive Y (up)
  nm = SEGM.Name;
  I  = SEGM.I;
  J  = SEGM.J;
  DX = SEGM.Dist;

% This flux calculation is for segments
% not for boxes where it is
% important to keep right direction of segments
% to close contours (i, j increasing)
  [b1,b2]=size(I);
  if b1<b2, I=I'; end;
  [b1,b2]=size(J);
  if b1<b2, J=J'; end;

  if J(end)==J(1),
    if I(end)<I(1)
      I=flipud(I);
      J=flipud(J);
    end
  end
  if I(end)==I(1)
    if J(end)<J(1)
      I=flipud(I);
      J=flipud(J);
    end
  end
  
  ll=length(I);
  clear FLX DZ
  FLX = zeros(nlr,ll-1);
  DZ  = zeros(nlr,ll-1);
  for l=1:ll-1   % little segments of a section
    i1=I(l);
    i2=I(l+1);
    j1=J(l);
    j2=J(l+1);
    dd=DX(l);  % segment length, m
% Find Normal:
% Sign convention:
    di=abs(i2-i1);
    dj=abs(j2-j1);
    if di>0 & dj>0
      error('Not step-like section ...');
    end
    if di==0 & dj==0
      error('Check I,J indices of the section ...');
    end

    H1  = HH(j1,i1);
    H2  = HH(j2,i2);
    flx = zeros(nlr,1); 
    U0  = [];
    V0  = [];
    if H1==0 & H2==0; FLX(:,l)=flx; continue; end;
%  Collocate U,V
% interpolate into the middle of the grid
    if di==0  % Y-section, V*norm=0
      U0    = squeeze(UN(:,j1,i1));
      dHj1  = squeeze(dH(:,j1,i1-1));
      dHj2  = squeeze(dH(:,j1,i1));
      U0(dHj1<1e-3) = 0;
      U0(dHj2<1e-3) = 0;
      U0(isnan(U0)) = 0;
      dH0   = 0.5*(dHj1+dHj2);
      flx   = U0.*dH0.*dd; % m3/s
    elseif dj==0
      V0    = squeeze(VN(:,j1,i1));
      dHi1  = squeeze(dH(:,j1-1,i1));
      dHi2  = squeeze(dH(:,j1,i1));
      dH0   = 0.5*(dHi1+dHi2);
      V0(dHi1<1e-3)=0;
      V0(dHi2<1e-3)=0;
      V0(isnan(V0))=0;
      flx   = V0.*dH0.*dd; % m3/s
    end

    FLX(:,l) = flx;
    DZ(:,l)  = -dH0;

  end;  % little segments along a section loop
  
%keyboard
  plt_sec=0;
  if plt_sec>0
    figure(10); clf;
    II=find(abs(DZ)<1e-6);
    ZZ=cumsum(DZ,1);
    ZZ(II)=nan;
    XX=cumsum(D)*1e-3;
    pcolor(XX,ZZ,FLX); shading flat
    caxis([-12000 12000]);
    colorbar
    stt=sprintf('%s, Volume Flux, m3/s',nm);
    title(stt);
  end
  
%  FLXZ=nansum(FLX);
  FlxTot = nansum(nansum(FLX))*1e-6; % total  flux

  fprintf('1 day %6.3f min\n',toc/60);
  fprintf('== %i/%2.2i/%2.2i, %s Flux: %6.2f Sv\n\n',DV(1:3),nm,FlxTot);
  FLXV(cnc,1)=FlxTot;

  if s_mat>0 & mod(ip,5)==0
    fprintf('####  Saving %s\n\n',fmat);
    save(fmat,'FLXV','TM','SEGM');
  end


end

if s_mat>0
  fprintf('####  Saving %s\n\n',fmat);
  save(fmat,'FLXV','TM','SEGM');
end

%save(fout,'UCPR');

%exit




