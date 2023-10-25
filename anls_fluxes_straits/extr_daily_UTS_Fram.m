% Extract time series of U, T, S
% along xsection with +/- Nrows  <--- Not implemented yet
%     just 1 section
% This is needed to cover location
% of moorings in the strait, for instance
%  1 year at a time
% Fields are on HYCOM grid and 
% NOT interpolated onto z-grid
%
%
% -------------------------------------
%              |
%            u -    *(i+1,j+1)
%              |
%              | 
%              *----|-----*--------*
%              |    v
%              |
%              |
% (i,j)      u -    * (i+1,j)
%              |
%              |
% ---|---------*----|-----*--------*
%   v               v
%
% Xsection is follow the edges of the grid cells
% i.e. goes across u,v-points 
% of the grid cells
% All Variables are interpolated into the 
% middle of the segments (betw i and i+1,  etc)
% ready for calculating fluxes across the small segments
% Variables: T, S, normal V, dH
% 

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt = 112;
TV   = 11;  % topo version
sctnm='Fram';

s_mat=1;
YR1=2004;
YR2=YR1;

rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

fprintf('Extracting U,T,S %s %i-%i\n\n',sctnm,YR1,YR2);

fmat = sprintf('%sarc08-%3.3i_UTS_daily_%s_%4.4i.mat',...
	       pthmat,expt,sctnm,YR1);

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
%[DX,DY]=sub_dx_dy(LON,LAT);
HH(HH>=0)=0;

UTS = struct;


% 79N:
IJ=[908         904
    915         904
         922         904
         929         904
         938         904
         945         905
         955         905
         966         906
         975         907
         985         909
         997         911
        1009         913
        1021         916
        1033         919
        1043         923
        1055         927
        1065         930
        1074         935
        1083         939];
nij=length(IJ);
IIs=[];
JJs=[];
for k=1:nij-1
  is=IJ(k,1);
  js=IJ(k,2);
  ie=IJ(k+1,1);
  je=IJ(k+1,2);
  [II,JJ]=sub_xsct_indx(is,js,ie,je);
  II=II(:);
  JJ=JJ(:);
  if k<nij-1
    IIs=[IIs;II(1:end-1)];
    JJs=[JJs;JJ(1:end-1)];
  else
    IIs=[IIs;II];
    JJs=[JJs;JJ];
  end
end

SCT.I=IIs;
SCT.J=JJs;
npt=length(IIs);

for ik=1:npt
  j=JJs(ik);
  i=IIs(ik);
  SCT.X(ik)=LON(j,i);
  SCT.Y(ik)=LAT(j,i);
end

IN=sub2ind([mm,nn],JJs,IIs);
[Im,Jm]=meshgrid([1:nn],[1:mm]);
SCT.Ind=IN;
Hb=HH(IN);
SCT.Hbottom=Hb;

% Depth, coordinates in the mid points, segment length
for ik=1:npt-1
  j1=JJs(ik);
  i1=IIs(ik);
  j2=JJs(ik+1);
  i2=IIs(ik+1);
  
  ln1=LON(j1,i1);
  lt1=LAT(j1,i1);
  ln2=LON(j2,i2);
  lt2=LAT(j2,i2);
  
% Segment length, m  
  dx=distance_spheric_coord(lt1,ln1,lt2,ln2);
  SCT.SGM_DX(ik)=dx;
  
% coord of the middle of the segment
  SCT.SGM_Lat(ik)=0.5*(lt1+lt2); 
  SCT.SGM_Lon(ik)=0.5*(ln1+ln2);
  
% depth at the middle point  
  hb=0.5*(Hb(ik)+Hb(ik+1));
  SCT.SGM_Hb(ik)=hb; 
    
% Norm, u or v
  if i1==i2
    SCT.Nrm(ik,:)=[1,0]; % u comp
  else
    SCT.Nrm(ik,:)=[0,1]; % v comp
  end
  
end
DX=SCT.SGM_DX;

% Interpolate onto fixed z-grid
% Specify fixed ZZ (interf) and ZM levels:
%ZZf = [(0:-1:-200)';...
%       (-202:-2:-1000)';(-2025:-25:-3000)';(-3050:-50:-5000)'];
%kzz = length(ZZf);


%dZf=diff(ZZf);
%ZMf = [];
%for ik=1:kzz-1
%  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
%end
%SCT.ZZ=ZZf;
%SCT.ZM=ZMf;


yr=YR1;
cc=0;
for iday=1:365
  dj1=datenum(yr,1,1);
  dnmb = dj1+iday-1; 
  DV   = datevec(dnmb);

%  pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
  pthbin = sprintf('/nexsan/hycom/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);

  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  if ~exist(fina,'file') | ~exist(finb,'file')
    fprintf('Not found %s or %s\n\n',fina,finb);
    continue;
  end

  cc=cc+1;

  fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);

  tic;
  [F,n,m,nlr] = read_hycom(fina,finb,'temp');
  F(F>hgg)=nan;
  T=F;

  [F,n,m,nlr] = read_hycom(fina,finb,'salin');
  F(F>hgg)=nan;
  S=F;
%
  [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
  F(F>hgg)=0;
  U=F;

  [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
  F(F>hgg)=0;
  V=F;

  [F,n,m,nlr] = read_hycom(fina,finb,'thknss');
  F=F./rg;
  F(F>hgg)=0;
  dH=F; 
%  [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
%
% Initialize arrays:
  Ts=squeeze(T(:,IN(1:end-1)))*0;
  Ss=squeeze(S(:,IN(1:end-1)))*0;
  Vs=squeeze(V(:,IN(1:end-1)))*0;
  dHs=squeeze(dH(:,IN(1:end-1)))*0;
  
% Go over small segments along the section
  for ik=1:npt-1
    i1=IIs(ik);
    i2=IIs(ik+1);
    j1=JJs(ik);
    j2=JJs(ik+1);
    dd=DX(ik);  % segment length, m
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

    h1 = HH(j1,i1);
    h2 = HH(j2,i2);
    
    if h1>=0 ; 
      Ts(:,ik)=NaN;
      Ss(:,ik)=NaN;
      Vs(:,ik)=NaN;
      continue; 
    end;

%
    if di==0  % Y-section, V*norm=0
      V0=U(:,j1,i1);
      dh1 = dH(:,j1,i1-1);
      dh2 = dH(:,j1,i1);
      dh0 = 0.5*(dh1+dh2);
      h0=0.5*(HH(j1,i1)+HH(j1,i1-1));
      t1  = T(:,j1,i1-1);
      t2  = T(:,j1,i1);
      s1  = S(:,j1,i1-1);
      s2  = S(:,j1,i1);
      Sdh = dh1+dh2;
      T0 = (dh1.*t1 + dh2.*t2)./Sdh;
      S0 = (dh1.*s1 + dh2.*s2)./Sdh;
      V0(dh0==0)=nan;
      T0(dh0==0)=nan;
      S0(dh0==0)=nan;
    elseif dj==0
      V0=V(:,j1,i1);
      dh1 = dH(:,j1-1,i1);
      dh2 = dH(:,j1,i1);
      dh0 = 0.5*(dh1+dh2);
      h0=0.5*(HH(j1,i1)+HH(j1-1,i1));
      t1  = T(:,j1-1,i1);
      t2  = T(:,j1,i1);
      s1  = S(:,j1-1,i1);
      s2  = S(:,j1,i1);
      Sdh = dh1+dh2;
      T0 = (dh1.*t1 + dh2.*t2)./Sdh;
      S0 = (dh1.*s1 + dh2.*s2)./Sdh;
      V0(dh0==0)=nan;
      T0(dh0==0)=nan;
      S0(dh0==0)=nan;
    end
    
    if (abs(nansum(dh0)-abs(h0))/abs(h0)) > 0.01
      fprintf('Check depths: segm i1=%i, j1=%i\n',i1,j1);
      fprintf(' h0=%8.2f, sum(dh)=%8.2f\n',abs(h0),nansum(dh0));
      keyboard;
    end
    
    Ts(:,ik)=T0;
    Ss(:,ik)=S0;
    Vs(:,ik)=V0;
    dHs(:,ik)=dh0;    
    
  end
  

%  dH=abs(diff(ZZh,1));
%  AA = sub_sct_interp2z(SCT,U,V,T,S,ZZf,ZZh,HH,dH,DX,DY);    
%keyboard
  SCT.TM(cc)=dnmb;
  SCT.Temp(cc,:,:)=Ts;
  SCT.Saln(cc,:,:)=Ss;
  SCT.Unrm(cc,:,:)=Vs;
  SCT.dH_thkn(cc,:,:)=dHs;

  fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);

  if mod(cc,30)==0 && s_mat==1
    fprintf('Saving %s\n',fmat);
    save(fmat,'SCT');
  end
  
  
end

if s_mat==1
  fprintf('Saving %s\n',fmat);
  save(fmat,'SCT');
end




