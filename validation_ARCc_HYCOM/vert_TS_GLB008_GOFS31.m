% Vertical section in the Arctic Ocean
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

dnmb = datenum(2015,07,01);
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

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat2 = '/Net/ocean/ddmitry/HYCOM/ARCc0.08/112/data_mat/';
btx='vert_TS_GLB008_GOFS31.m';


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

expt = 563;
pthycom = sprintf('/nexsan/hycom/GLBv0.08_%3.3i/data/hindcasts/%4.4i/',expt,yr); 
fnm  = sprintf('%shycom_GLBv0.08_%3.3i_%4.4i%2.2i%2.2i12_t000.nc',...
                pthycom,expt,yr,imo,iday);

fprintf('Opening %s\n',fnm);

LON = nc_varget(fnm,'lon');
LAT = nc_varget(fnm,'lat');
Z   = nc_varget(fnm,'depth');
if strncmp(pfld,'salin',4)
  F = squeeze(nc_varget(fnm,'salinity'));
else
  F = squeeze(nc_varget(fnm,'water_temp'));
end
Lmsk = squeeze(nc_varget(fnm,'surf_el'));
Lmsk = Lmsk*0+1;
Lmsk(isnan(Lmsk))=0;


%
% Get coordinates of the section line saved in vert_TS_xsect008.m
fout=sprintf('%sarc008_xsect_%s.mat',pthmat2,xname);
fprintf('Loading %s\n',fout);
XYs=load(fout);

Xv=XYs.Xl;
Yv=XYs.Yl;
nv=size(Xv,1);
cc=0;
for in=1:nv
  x0=Xv(in);
  y0=Yv(in);


  DI=abs(LON-x0);
  ii=find(DI==min(DI),1);
  DJ=abs(LAT-y0);
  jj=find(DJ==min(DJ),1);

  if cc==0
    cc=cc+1;
    IJs(cc,1)=ii;
    IJs(cc,2)=jj;
  else  % check previous point
    di=abs(IJs(cc,1)-ii);
    dj=abs(IJs(cc,2)-jj);
    if di>0 | dj>0
      cc=cc+1;
      IJs(cc,1)=ii;
      IJs(cc,2)=jj;
    end
  end

end
nS=cc;

f_map=0;
if f_map==1
  figure(10); clf;
  contour(Lmsk,[0.01 0.01],'k');
  hold on;
  plot(IJs(:,1),IJs(:,2),'k.','Markersize',10);
end

% Distances:
for ii=1:nS
  i0=IJs(ii,1); 
  j0=IJs(ii,2); 
  Xl(ii)=LON(i0);
  Yl(ii)=LAT(j0);
end

x1=Xl(1);
y1=Yl(1);
for ii=1:nS
  x2=Xl(ii);    
  y2=Yl(ii);    
  dXX(ii)=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
  x1=x2;        
  y1=y2;        
end
XL1=cumsum(dXX);

% Get section:
clear AA LTs LNs
for ii=1:nS
  i0=IJs(ii,1);
  j0=IJs(ii,2);
  AA(:,ii)=squeeze(F(:,j0,i0));
  LTs(ii)=LAT(j0);
  LNs(ii)=LON(i0);
end

stl=sprintf('GOFS3.1-%i, %s, %4.4i-%2.2i, %s',...
            expt,fld0,yr,imo,xname);                                                                 
nf=1;
xdr = 1;
f_layer=0;
[ma,na]=size(AA);
s0=[];
if strncmp(pfld,'salin',4)
% Convert to log scale for better depiction of high-value S
  A0=AA;
                
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

dmm=ones(length(XL1),1);
ZZ=-Z*dmm';
dmm=ones(length(Z),1);
XL=dmm*XL1;

sub_plot_xsct3_btm(nf,XL,ZZ,AA,stl,fld0,f_layer,...
                   xl1,xl2,yl1,yl2,c1,c2,cntr,s0);


bottom_text(btx,'pwd',1);






