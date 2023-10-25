% Calculate S flux - monthly means
% across specified sections
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat = 1; 

expt = 110;
TV   = 11;  % topo version
rg  = 9806;
hgg = 1e20; % 

segm = 'AtlOB';

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);
dd  = 15;
YR1 = 1993;
YR2 = 2016;

fmat = sprintf('%s%3.3i_SFlux_%s_%4.4i_%4.4i.mat',...
	       pthmat,expt,segm,YR1,YR2);

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


switch(segm)
 case('AtlOB');
  is1=233;
  js1=53;
  is2=1019;
  js2=js1;
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
D(ll)=D(ll-1);
SEGM.Dist=D;

nrec=0;  % # of saved mean U,V averaged over Nav days
cnc=0;
ip1=1;

if s_mat==2
  load(fmat);
  cnc = length(TM);
  dlst = TM(cnc);
  ip1 = find(YRPLT(:,3)==dlst)+1;
  if isempty(ip1),
    error('Couldnt find date in YRPLT %s',datestr(dlst));
  end
  
end

dday = 15;
for yr=YR1:YR2
  pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
  SF = 0;
  VF = 0;
  for im=1:12
    nmm=0;
    for md=1:dday:31
      tic;
      dnmb = datenum(yr,im,md);
      DV = datevec(dnmb);
      mo=DV(2);
      mday=DV(3);
      if mo~=im, continue; end;
      iday = dnmb-datenum(yr,1,1)+1;

      fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
      finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

      if ~exist(fina,'file');
        fprintf('Not found: %s\n\n',fina);
        continue;
      end
  
%      nrec=nrec+1;
      fprintf('%4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

% Layer thickness:
%      [ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',1);
      [F,n,m,nlr] = read_hycom(fina,finb,'thknss');
      F=F./rg;
      F(F>hgg)=0;
      dH=F; 
      
      if ~exist('Dst','var')
	for kk=1:nlr
	  Dst(kk,:)=D;
	end
      end
      
    
% Calculate transport:
% in archm, archv: u, v anre not collocated 
% with H, dH:
% u(i,j) is between h(i-1,j) and h(i,j)
% To check look at near land points
% such that HH(j-1,i) and HH(j,i-1) is land and HH(j,i) is not
% if V(j,i) is nan  than V is not collocated with H
% if U(j,i) is nan than U is not collocated with H
%need barotropic U?
      [F,n,m,ll] = read_hycom(fina,finb,'salin');
      F(F>1e6)=nan;
      S2=squeeze(F(:,js1,is1:is2));
      S1=squeeze(F(:,js1-1,is1:is2));
      
      [F,n,m,ll] = read_hycom(fina,finb,'v-vel.');
      F(F>1e6)=nan;
      V=squeeze(F(:,js1,is1:is2));
      
      dHi1  = squeeze(dH(:,js1-1,is1:is2));
      dHi2  = squeeze(dH(:,js1,is1:is2));
      dH0   = 0.5*(dHi1+dHi2);
      
      S = (S1.*dHi1+S2.*dHi2)./(dHi1+dHi2);
      
      V(dHi1<1e-3)=0;
      V(dHi2<1e-3)=0;
      I=find(isnan(V));
      if ~isempty(I)
        V(I)=0;
      end
      Vflx = V.*dH0.*Dst; % m3/s
      vf = nansum(nansum(Vflx)); % total vol flux m3/s
      VF = VF+vf;
      
% Salt flux calculated as mass flux
% where salinity 1 = 0.001 kg/m3 of salt
      Sflx = S.*Vflx;  
      sf = nansum(nansum(Sflx)); %kg/s
      SF = SF+sf;
      
      nmm=nmm+1;
      
      fprintf('Vol Flux = %8.4f Sv, Sflux = %8.4f kg\n',VF/nmm*1e-6,SF/nmm);
      
      
      fprintf('1 day %6.3f min\n',toc/60);
      
    end  % days
    
    nrec=nrec+1;

    SEGM.TM(nrec,1)    = dnmb;
    SEGM.Vol_Flx(nrec) = VF/nmm;
    SEGM.S_Flx(nrec)   = SF/nmm;
      
  end
  if s_mat==1;
    fprintf('Saving %s\n',fmat);
    save(fmat,'SEGM');
  end
end

if s_mat==1;
  fprintf('Saving %s\n',fmat);
  save(fmat,'SEGM');
end

      
%       
      
