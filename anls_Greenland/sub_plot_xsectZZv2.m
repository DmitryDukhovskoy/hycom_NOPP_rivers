function sub_plot_xsectZZv2(nf,xx,zz,A,title_str,fld0,xdr);
% plot vertical section of tracer, density, etc
% interpolated into regular vertical grid ZZ
%
ZZ=zz;
plr=20; % highlight this interface in HYCOM

% Colormap:
switch(fld0)
 case('salin');
    na=50;
    cl1=colormap_gray(na);
    cl2=colormap_purple(na);
    cl3=colormap_blue(na);
    cl4=colormap_green(na);
    cl5=colormap_yellow(na);
    cl6=colormap_red(na);

    cmp=[cl1;cl2;cl3;cl4;cl5;cl6];
    cmp=smooth_colormap(cmp,20);
    nint=length(cmp);
%    c1=32.5;
%    c2=34.9;
    c1=33.5;
    c2=35.;
 case('temp');
    na=100;
    cl1=colormap_gray(na);
    cl2=flipud(colormap_purple(na));
    cl3=flipud(colormap_blue(na));
    cl4=colormap_green(na);
    cl5=colormap_yellow(na);
    cl6=colormap_red(na);

%    cmp=[cl1;cl2;cl3;cl4;cl5;cl6];
    cmp=[cl3;cl6];
%    cmp=smooth_colormap(cmp,20);
    cmp=smooth_colormap(cmp,20);
    nint=length(cmp);
%    c1=32.5;
%    c2=34.9;
    c1=-2;
    c2=2;
  case('tracer');
%  cmp=colormap_cold(350);
%  nint=length(cmp);
%  c1=-5;
%  c2=2;
  CMP = sub_conc_clrmp01;
  cmp=CMP.colormap;
  cnt=CMP.intervals;
  nint=length(cmp);
  c1=CMP.c1;
  c2=CMP.c2;
 case('brunt'); % Brunt-V. freq., log scale
  CMP = sub_conc_clrmp02(-8,-4);
  cmp=CMP.colormap;
  cnt=CMP.intervals;
  nint=length(cmp);
  c1=CMP.c1;
  c2=CMP.c2;
  
end;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals



figure(nf);
clf
% 
axes('position',[0.1 0.1 0.8 0.8]);
%xx=XL;
%zz=ZZav;
%A=lTr;
% Find bottom:
[ll,nn]=size(A);

for ii=1:nn
  dmm = zz(:,ii);
  ddm = diff(dmm);
  k = find(ddm==0,1);
  Zbt(ii,1) = zz(k+1,ii);
end
zb=min(Zbt);
fprintf('Bottom zb=%6.1f\n',zb);
%keyboard
pcolor(xx,ZZ,A); shading interp;
%pcolor(xx,zz,A); shading flat;
caxis([c1 c2]);
hold on;
colormap(cmp);

%keyboard

x=xx(1,:);
set(gca,'Color',[0 0 0],'tickdir','out');
set(gca,'xlim',[min(min(xx)) max(max(xx))],'ylim',[1.02*zb 0]);
if xdr<0
  set(gca,'xdir','reverse');
end
set(gca,'xtick',[0:50:max(x)],'ytick',[-1000:50:0]);
set(gca,'fontsize',16);

%  lt0=mean(LAT(j1,i1:i2));
%text(xx(8),200,title_str,'Fontsize',16);
title(title_str,'Fontsize',16);

hght=[];
lngth=[];
mint=20;
mbx=20;
fsz=14;
bxc='k';
posc=[0.91 0.1 0.8 0.05];
aend=1;
[az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);



return
