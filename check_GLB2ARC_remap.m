% Read instant output archive files:
% check remapped files:
% GLBb T07 -> ARCc T07
% ARCc T07 -> ARCc T11
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear


pthbin1 = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
pthbin1 = '/Net/kronos/ddmitry/hycom/ARCb0.08/output/';
pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_73.7/';  
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

rg=9806;  % convert pressure to depth, m
yr=2019;
iday=84;
hr=0;
TV=11;

% Topo
fltopo=sprintf('%sdepth_ARCc0.08_11.nc',pthtopo);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn] = size(HH);
HH(HH>0)=nan;


%fld='thknss';
fld='temp';
%fld='salin';
%fld = 'bl_dpth';
c1=-2;
c2=18;

% remapping produced on Gordon
fina = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_%2.2i_stp1.a',pthbin1,TV,yr,iday,hr);
finb = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_%2.2i_stp1.b',pthbin1,TV,yr,iday,hr);
[F,n,m,l] = read_hycom(fina,finb,fld,'r_layer',nlr);
F(F>1e20)=nan;
A1=squeeze(F);
figure(1); clf;
pcolor(A1); shading flat
caxis([c1 c2]);
colorbar
stl=sprintf('Gordon ARCb 0.08 %s',fld);
title(stl);

% Remapping produced on COAPS machine
fina = sprintf('%sarchv_arcT%2.2iL41_stp1.%4.4i_%3.3i_12.a',pthbin,TV,yr,iday);
finb = sprintf('%sarchv_arcT%2.2iL41_stp1.%4.4i_%3.3i_12.b',pthbin,TV,yr,iday);


[F,n,m,l] = read_hycom(fina,finb,fld,'r_layer',nlr);
F(F>1e20)=nan;
A=squeeze(F);
figure(3); clf;
pcolor(A); shading flat
caxis([c1 c2]);
colorbar
stl=sprintf('COAPS ARCb 0.08 %s',fld);
title(stl);


gina = sprintf('%s737_archm.%i_%3.3i_12.a',pthglb,yr,iday);
ginb = sprintf('%s737_archm.%i_%3.3i_12.b',pthglb,yr,iday);


[F,n,m,l] = read_hycom(gina,ginb,fld,'r_layer',nlr);
F(F>1e20)=nan;
Ag=squeeze(F);
figure(2);
pcolor(Ag); shading flat
caxis([-2 18]);
set(gca
colorbar
title(fld);




%keyboard

% Section:
jj1=1300;
jj2=jj1;
ii1=200;
ii2=1400;
Tr = squeeze(F(:,jj1,ii1:ii2));
[a1,a2]=size(Tr);


fld='thknss';
[F,n,m,l] = read_hycom(fina,finb,fld);
F(F>1e10)=nan;
F(F<0.1)=nan;
F=F./rg;
Dsec=squeeze(F(:,jj1,ii1:ii2)); 
Dsec(Dsec==0)=nan;

% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
clear ZZb
Dsec(isnan(Dsec))=0;
ZZb(1,:)=-Dsec(1,:);
for kk=2:l
  ZZb(kk,:)=ZZb(kk-1,:)-Dsec(kk,:);
end

% For plotting need to have #of layers +1 interfaces
% add surface layer, otherwise all values
% will be shifted 1 layer down
[nl,npb]=size(ZZb);
ZZ=zeros(nl+1,npb);
ZZ(2:nl+1,:)=ZZb;
Tr(nl+1,:)=Tr(nl,:);

Nintf=size(ZZ,1);
for ii=1:a2
  I=ii1+ii-1;
  XL(1:Nintf,ii)=ones(Nintf,1)+I;
end

figure(2); clf;
pcolor(XL,ZZ,Tr); shading flat;
hold on;
caxis([c1 c2]);
hold on;
colorbar

zb=min(min(ZZ));
x=XL(1,:);
set(gca,'Color',[0 0 0],'tickdir','out');
set(gca,'xlim',[min(x) max(x)],'ylim',[1.02*zb 0]);
set(gca,'xtick',[0:200:max(x)],'ytick',[-4500:500:0]);
set(gca,'fontsize',16);



