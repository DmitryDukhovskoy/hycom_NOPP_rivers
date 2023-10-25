% Extract T and Layer thicknes of Atl. Water
% in specified regions (Canada Basin, BG, ...)
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps/
startup

format long g
clear all
close

f_mat=1;

regn = 'ARCc0.08';
expt = 110;

pfld  = 'temp';
%f_extr = 1;  % =0 - load in extracted depth of Atl. Water, =1 -extract
s_fig  = 0;

sfig=0;

rg = 9806;

pthfig  = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/%3.3i/fig_AtlLayer/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);

YRPLT=[];
cc=0;
for iyr=1993:2016
  for im=1:12
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=im;
    YRPLT(cc,3)=15;
  end
end
nrc=cc;

% Get topo:
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn] = size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY;

Adpth=zeros(mm,nn)*nan;
nxx=0;
%mm=2520;
%nn=1600;
%ll=33;

% Region: 
% Deep Canada Basin:
IJ=[     466        1524
         479        1628
         601        1660
         691        1516
         804        1411
         679        1247
         530        1460];

nij=length(IJ);

[X,Y]=meshgrid([1:nn],[1:mm]);
IX=inpolygon(X,Y,IJ(:,1),IJ(:,2));
IN=find(IX==1 & HH<-1000);
nI = length(IN);

f_bth=0;
if f_bth==1
figure(1); clf;
contour(HH,[0 0],'k');
hold on;
contour(HH,[-5000:1000:-100],'b');
plot(X(IN),Y(IN),'r.');
axis('equal');
end

icc=0;
for ik=1:nrc
  yr=YRPLT(ik,1);
  mo=YRPLT(ik,2);
  md=YRPLT(ik,3);
  dnmb=datenum(yr,mo,md);
  dv0=datevec(dnmb);
  iday=dnmb-datenum(yr,1,1)+1;
%if f_extr == 1
  fprintf('Reading %i/%2.2i/%2.2i\n',dv0(1:3));
  pthbin = sprintf('/nexsan/archive/ARCc0.08_110/data/%4i/',yr);
 
  dnmb=datenum(yr,1,1)+iday-1;
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  icc=icc+1;
  
  if ~exist('ZM','var');
    [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
  end
  
 % Layer thickness:
%  [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
%  F=squeeze(F(plr,:,:));
%  F=F./rg;
%  F(F>1e10)=0;
%  F(F<1e-2)=0;
%  dH=squeeze(F); 

  tic;
  [F,n,m,l] = read_hycom(fina,finb,pfld);
  F(F>1e6)=nan;
  T = F;

  Zmax = HH*nan; % depth of max T
  Zt0  = HH*nan; % depth of upper 0
  Tmax = HH*nan; % mean max T
  dHatl= HH*nan; % thickness of Atl. Layer
  Tav  = HH*nan; % mean T of Atl. Layer
  
  for i=1:nI
    i0 = IN(i);

    tt  = squeeze(T(:,i0));
    tt0 = squeeze(T(:,i0));
    z0  = ZZ(:,i0);
    iz0 = max(find(z0>-100));
    iz2 = max(find(z0>-900));
    if ~isempty(iz0)
      tt(1:iz0) = -999;
      tt(iz2:end) = -999;
      tt0(1:iz0) = -999;
    end
    
%    tt(1:10) = -10;
    tm = max(tt);  % max T within the layer:
    Tmax(i0) = tm;

    iz = find(tt==max(tt),1);
    Zmax(i0) = ZM(iz,i0);
    iz = min(find(tt>=0));
%keyboard
% Depth max T
% Interpolate to find exact z
% upper 0 interface
    iz1=iz-1;
    iz2=iz;
    t1 = tt(iz1);
    t2 = tt(iz2);
    z2 = ZZ(iz2,i0);
    z1 = ZZ(iz1,i0);
    t0 = 0;
    dtdz = (t2-t1)/(z2-z1);
    dz0 = (t0-t1)/(dtdz);
    zt0 = z1+dz0;
    Zt0(i0) = zt0;
    
% Thickness of warm Atl. Layer
% Interpolate to get exact thickness
    Ipz = find(tt0>0);
    iz1=Ipz(1)-1;
    iz2=Ipz(1);
    z1=ZM(iz1,i0);
    z2=ZM(iz2,i0);
    t2=tt0(iz2);
    t1=tt0(iz1);
    t0=0;
    dtdz = (t2-t1)/(z2-z1);
    dz0 = (t0-t1)/(dtdz);
    zt0 = z1+dz0;
    zup = zt0;
    dzup = abs(dz0);
    
    iz1=Ipz(end);
    iz2=iz1+1;
    z1=ZM(iz1,i0);
    z2=ZM(iz2,i0);
    t2=tt0(iz2);
    t1=tt0(iz1);
    dtdz = (t2-t1)/(z2-z1);
    dz0 = (t0-t1)/(dtdz);
    zbt = z1+dz0;
    dzbt= abs(dz0);
    
    dh = abs(zbt-zup); % Atl.L. thickness
    
    dHatl(i0)=dh;
    
% Mean T within the Atl. Layer:
    tt0 = squeeze(T(:,i0));
    iz1 = Ipz(1)-1;
    iz2 = Ipz(end)+1;
    dzz = diff(ZZ(:,i0));
    t1 = 0.5*(t0+tt0(Ipz(1)));
    t2 = 0.5*(t0+tt0(Ipz(end)));
    dZ = abs([dzup;dzz(Ipz);dzbt]);
    tmm = [t1;tt0(Ipz);t2]; 
    tav = (tmm'*dZ)./sum(dZ);
    Tav(i0)=tav;

% Plot T;
    chckT=0;
    if chckT==1
      figure(10); clf;
      plot(T(:,i0),ZM(:,i0),'.-');
      hold on;
      plot([0 0],[min(ZM(:,i0)) 0],'r--');
      plot(tt0(Ipz),ZM(Ipz,i0),'go');
      plot(t0,zt0,'m*'); % located 0
      plot(t0,zbt,'m*'); % located 0
      plot([tav tav],[ZM(iz2,i0) ZM(iz1,i0)],'g--');
    end
    
  end
  
  ATLW.dnmb(icc,1)   = dnmb;
  ATLW.Tmax(icc,1)   = nansum(Tmax(IN).*Acell(IN))./sum(Acell(IN)); % max T
  ATLW.Z_T0(icc,1)   = nansum(Zt0(IN).*Acell(IN))./sum(Acell(IN));  % Depth of max T
  ATLW.Z_Tmax(icc,1) = nansum(Zmax(IN).*Acell(IN))./sum(Acell(IN));  % Depth of max T
  ATLW.dHatl(icc,1)  = nansum(dHatl(IN).*Acell(IN))./sum(Acell(IN)); % AtlL thckn
  ATLW.Tatl_av(icc,1)= nansum(Tav(IN).*Acell(IN))./sum(Acell(IN)); % avrg T atl L. 
  
  tmx = ATLW.Tmax(icc,1);
  zmx = ATLW.Z_Tmax(icc,1);
  z0c = ATLW.Z_T0(icc,1);
  dha = ATLW.dHatl(icc,1);
  tav = ATLW.Tatl_av(icc,1);
  
  if f_mat==1 & mod(icc,12)==0
    fmat=sprintf('%sarc08_110_atlw_TdH.mat',pthmat);
    fprintf('saving %s\n',fmat);
    save(fmat,'ATLW');
  end
  
  fprintf('-----------------------------------\n')
  fprintf('AvT=%6.2f, Dpth_0C=%6.2f m, maxT=%6.2f\n',...
	  tav,z0c,tmx);
  fprintf('dpthTmax=%6.2f m, Atl Thck=%6.2fm\n',...
	  tmx,dha);
  fprintf('-----------------------------------\n')
	  
  fprintf('Processed 1 record %6.4f min\n\n',toc/60);
%keyboard

end




    








