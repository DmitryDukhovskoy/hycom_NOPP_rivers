% Vertical section in the Arctic Ocean
% GOFS3.5 reanalysis
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

dnmb = datenum(2017,07,01);
DV = datevec(dnmb);


s_fig = 0;
pfld  = 'temp';
%pfld  = 'salin';
fld0  = pfld;
%xname = 'BerNatl'; % Bering - N.Atl
%xname = 'BeaufNatl'; % Beauf. Shelf - South Norway
xname = 'BeaufIcel';
%xname = 'LaptSea'; % Laptev Sea

hg    = 2^100;
rg    = 9806;
% Change Tv and nlev for different ARCc fields
Tv    = 11; % bathym v11
nlev  = 41; % 41

pthmat2 = '/Net/ocean/ddmitry/HYCOM/ARCc0.08/112/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);
btx='vert_TS_GLB004_GOFS35.m';


ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LN04 = nc_varget(ftopo,'Longitude');
LT04 = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);


% Section:
switch(xname);
  case('BeaufIcel');
%
%  cs1=30;
%  cs2=35;  
  cs1=3.36;  % log scale S=28.79
  cs2=3.56;  % log scale S=35.16
  ct1=-2;
  ct2=5;
end


yr   = DV(1);
imo  = DV(2);
iday = DV(3);
jday = dnmb-datenum(yr,1,1)+1;

expt = 216;
pthycom = '/Net/kronos/ddmitry/hycom/GOFS3.5/';
fina  = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthycom,expt,yr,jday);
finb  = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthycom,expt,yr,jday);


% Get ARCc grid indices
flgrd = sprintf('%sregional.grid',pthycom);
fltopo = sprintf('%sdepth_GLBc0.04_27.a',pthycom);
GRD = read_grid_bath(flgrd,fltopo);
pln=GRD.PLON;
plt=GRD.PLAT;
hh=GRD.Topo;
I=find(hh>1e20);
hh=-1*hh;
hh(I)=100;

ln11=LN04(1,1);
ln21=LN04(end,1);
ln12=LN04(1,end);
ln22=LN04(end,end);

lt11=LT04(1,1);
lt21=LT04(end,1);
lt12=LT04(1,end);
lt22=LT04(end,end);

dd=sqrt((pln-ln11).^2+(plt-lt11).^2);
[j11,i11]=find(dd==0);
dd=sqrt((pln-ln12).^2+(plt-lt12).^2);
[j12,i12]=find(dd==0);
dd=sqrt((pln-ln21).^2+(plt-lt21).^2);
[j21,i21]=find(dd==0);
dd=sqrt((pln-ln22).^2+(plt-lt22).^2);
[j22,i22]=find(dd<1e-5);

%plot(i11,j11,'r.','Markersize',12)
%plot(i12,j12,'g.','Markersize',12)
%plot(i21,j21,'c.','Markersize',12)
%plot(i22,j22,'m.','Markersize',12)

IJ=[i11, j11; ...
    i12, j12; ...
    i21, j21; ...
    i22, j22];

LON = sub_Glb2Arc(pln,IJ);
LAT = sub_Glb2Arc(plt,IJ);
HH  = sub_Glb2Arc(hh,IJ);

%
% Specify sections:
switch(xname);
 case('BeaufIcel');
  XYs=[-144.349166870117          69.8827896118164;...
       157.7998046875             89.9287109375; ...
       1.6041259765625          76.6115264892578; ...
       -15.9317932128906          66.0183715820312];

  nx=size(XYs,1);
  for kx=1:nx
    x0=XYs(kx,1);
    y0=XYs(kx,2);
    dd=distance_spheric_coord(y0,x0,LAT,LON);
    [j,i]=find(dd==min(min(dd)));
    IJs(kx,1)=i;
    IJs(kx,2)=j;
  end

  xl1=1;
  xl2=6735;
  yl1=-4500;
  yl2=0;
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


Lmsk=HH(:,:);
I=find(isnan(Lmsk));
Lmsk=Lmsk*0+1;
Lmsk(I)=0;

f_map=0;
if f_map==1
  figure(10); clf;
  contour(Lmsk,[0.01 0.01],'k');
  hold on;
  plot(IJs(:,1),IJs(:,2),'k.','Markersize',10);
end



% ---------------------------------
% Read GLBc 0.04 GOFS3.5
% Read by layers
% subsample ARCc region
% and rotate the grid
% ---------------------------------
fprintf('Reading %s\n',datestr(dnmb));

ll=41;
SS=[];
%pfld='salin';
if strncmp(pfld,'salin',4);
  for kk=1:ll
    fprintf('Reading %s layer %i\n',pfld,kk);
    [Fr,n,m,l1] = read_hycom(fina,finb,pfld,'r_layer',kk);
    Fr=squeeze(Fr);
    Fr(Fr>hgg)=nan;
    dmm = sub_Glb2Arc(Fr,IJ);
    SS(kk,:,:)= dmm;
  end
  S = squeeze(SS(:,INDs));
  [a1,a2]=size(S);
  % Prepare for plotting - add extra bogus layer
  % at the bottom, Matlab won't show it
  S(end+1,1:a2)=S(end,:);
  AA=S;
  A0=S;
end


TT=[];
%pfld='temp';
if strncmp(pfld,'temp',4)
  for kk=1:ll
    fprintf('Reading %s layer %i\n',pfld,kk);
    [Fr,n,m,l1] = read_hycom(fina,finb,pfld,'r_layer',kk);
    Fr=squeeze(Fr);
    Fr(Fr>hgg)=nan;
    dmm = sub_Glb2Arc(Fr,IJ);
    TT(kk,:,:)= dmm;
  end
  T = squeeze(TT(:,INDs));
  [a1,a2]=size(T);
  % Prepare for plotting - add extra bogus layer
  % at the bottom, Matlab won't show it
  T(end+1,1:a2)=T(end,:);
  AA=T;
end
   
dH=[];
pf='thknss';
for kk=1:ll
  fprintf('Reading %s layer %i\n',pf,kk);
  [Fr,n,m,l1] = read_hycom(fina,finb,pf,'r_layer',kk);
  Fr=squeeze(Fr);
  Fr(Fr>hgg)=nan;
  Fr(Fr<0.1)=nan;
  Fr=Fr/rg;
  dmm = sub_Glb2Arc(Fr,IJ);
  dH(kk,:,:)= dmm;
end
l=ll;
Dsec=squeeze(dH(:,INDs));
Dsec(Dsec==0)=nan;
% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
clear ZZb
Dsec(isnan(Dsec))=0;
ZZb(1,:)=-Dsec(1,:);
for kk=2:ll
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

%
% Distances:

x1=Xl(1);
y1=Yl(1);
for ii=1:nS
  x2=Xl(ii);    
  y2=Yl(ii);    
  dXX(ii)=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
  x1=x2;        
  y1=y2;        
end
I=find(dXX==0);
dXX(I)=0.01;
XL1=cumsum(dXX);

dmm=ones(ll+1,1);
XL=dmm*XL1;


stl=sprintf('GOFS3.5-%i, %s, %4.4i-%2.2i, %s',...
            expt,fld0,yr,imo,xname);                                                                 


nf=1;
xdr = 1;
f_layer=0;
[ma,na]=size(AA);
s0=[];
if strncmp(pfld,'salin',4)
% Convert to log scale for better depiction of high-value S
%  A0=AA;
                
  s0=35.6;
  flg=1;        
  AA=salin2exp(s0,A0,flg);
%  c1=min(min(AA));
%  c2=max(max(AA));
  c1=1.1316;  % to match arc008
  c2=5.56293; 
  cntrS=[34.8:0.05:35];
  cntr=exp((35.6-cntrS).^(-1));
else            
  c1=ct1;
  c2=ct2;       
  cntr=[-1.:0.5:4];
end             


xl1=1;
xl2=(max(XL1));
yl1=-4500;
yl2=0;

sub_plot_xsection3(nf,XL,ZZ,AA,stl,pfld,f_layer,...
                     xl1,xl2,yl1,yl2,c1,c2,cntr,s0);

bottom_text(btx,'pwd',1);






