% Compare target densities 
% in 32 vs 41 layers
% Use climatology relax fields interpolated into
% 32 abd 41 layers, see RELAX/relax_arc08.com
% 
% Do not use archive files from different 
% simulations with 32 and 41 layers as
% differences in topo and density fields ->
% not consistent vertical layers
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

ptharc  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
pthbin  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/110/';
pthglb  = pthbin;
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

rg=9806;  % convert pressure to depth, m
yr=1993;
iday=1;
hr=0;
TV=11;
mo=7; % month, 1 or 7

% Topo:
fltopo=sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn] = size(HH);
HH(HH>0)=nan;


% 32 layers, TOPO 11
% HYCOM output form GLBa
%fina = sprintf('%sarchv_arcT%2.2i.%4.4i_%3.3i_%2.2i.a',ptharc,TV,yr,iday,hr);
%finb = sprintf('%sarchv_arcT%2.2i.%4.4i_%3.3i_%2.2i.b',ptharc,TV,yr,iday,hr);
% Climatology 
fina = sprintf('%s32layers_T%2.2i/relax.m%2.2i.a',pthbin,TV,mo);
finb = sprintf('%s32layers_T%2.2i/relax.m%2.2i.b',pthbin,TV,mo);

% Sections "West-East", i.e. fixed j:
%sct = 'WE_CArct';
sct = 'Pacif_OB2';
%sct = 'Pacif_OB';
switch(lower(sct))
 case('we_carct');
% Central Arctic:
  jj1=1430;
  jj2=jj1;
  ii1=200;
  ii2=1450;
 case('atl_ob')
% Atl. OB:
  jj1=5;
  jj2=jj1;
  ii1=190;
  ii2=1050;
 case('pacif_ob');
% Pacific OB:
  jj1=2516;
  jj2=jj1;
  ii1=1;
  ii2=nn;
 case('pacif_ob2');
% Pacific OB:
  jj1=mm-1;
  jj2=jj1;
  ii1=1;
  ii2=nn;
end



%fld='thknss';
fld='temp';

cmp=[1 1 1];
[F,n,m,l] = read_hycom(fina,finb,fld);
F(F>1e20)=nan;
%T = squeeze(F(1,:,:));
%hold on;
%plot([ii1 ii2],[jj1 jj2],'m-');

A = squeeze(F(:,jj1,ii1:ii2));
[a1,a2]=size(A);
A(~isnan(A))=1;

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


Nintf=size(ZZb,1);
for ii=1:a2
  I=ii1+ii-1;
  XL(1:Nintf,ii)=ones(Nintf,1)+I;
end


% Get 41 layers:
%finRa=sprintf('%srestart_GLB2ARC_T11_41lr_108a.a',pthglb);
%finRb=sprintf('%srestart_GLB2ARC_T11_41lr_108a.b',pthglb);
finRa = sprintf('%s41layers_T%2.2i/relax.m%2.2i.a',pthbin,TV,mo);
finRb = sprintf('%s41layers_T%2.2i/relax.m%2.2i.b',pthbin,TV,mo);

ID=n;
JD=m;
fld='thknss';
%fld='dp';
%Fg = read_hycom_restart(finRa,finRb,fld,ID,JD);
[Fg,nr,mr,lr] = read_hycom(finRa,finRb,fld);
lr=size(Fg,1);

Fg(Fg>1e10)=nan;
Fg(Fg<0.1)=nan;

Fg=Fg./rg;
Dr=squeeze(Fg(:,jj1,ii1:ii2)); 
Dr(Dr==0)=nan;

% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
clear Zr
Dr(isnan(Dr))=0;
Zr(1,:)=-Dr(1,:);
for kk=2:lr
  Zr(kk,:)=Zr(kk-1,:)-Dr(kk,:);
end



figure(1); clf;

pcolor(XL,ZZb,A); shading flat
hold on;

set(gca,'color',[0 0 0]);

lr00=24; % layer below which isopycnals are the same in 32 and 41
for kk=1:lr
  zz=ZZb(kk,:);
  xx=XL(kk,:);
  if mod(kk,5)==0;
    plot(xx,zz,'m-');
  elseif kk==lr00
    plot(xx,zz,'g-');
  else
    plot(xx,zz,'w-');
  end    
end


% FIeld with 41 layers:
lrXX=24;
for kk=1:lr
  zz=Zr(kk,:);
  xx=XL(1,:);
  if mod(kk,5)==0;
    plot(xx,zz,'k--','Color',[0.4 0.4 0.4]);
  elseif kk==lrXX
    plot(xx,zz,'c--');
  else
    plot(xx,zz,'k--','Color',[0.7 0.7 0.7]);
  end    
end
stt=sprintf('%s, 32(blue/red)vs 41(black/grey) layers, ARCc T%2.2i, HPC mo=%2.2i',...
	   sct,TV,mo);
title(stt);

txtbt='target_dens_32v41.m';
bottom_text(txtbt,'pwd',1);


