function sub_plot_xsection(nf,xx,zz,A,title_str,fld0,f_layer);
% called from xsection_isopycnals_atlantic.m
% to plot vertical section of tracer conc.
% with or without v.layers
%
ZZ=zz;
plr=20; % highlight this interface in HYCOM

% Colormap:
switch(fld0)
 case('salin');
  nint=200;
  cmp = colormap(parula(nint));

  c1=30;
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
    cmp=[cl3;cl4;cl5;cl6];
%    cmp=smooth_colormap(cmp,20);
    cmp=smooth_colormap(cmp,20);
    nint=length(cmp);
%    c1=32.5;
%    c2=34.9;
    c1=-2;
    c2=6;
  case('tracer');
%  cmp=colormap_cold(350);
%  nint=length(cmp);
%  c1=-5;
%  c2=2;
  CMP = sub_conc_clrmp01;
  cmp=CMP.colormap;
  cnt=CMP.intervals;
  nint=length(cmp);
%  c1=CMP.c1;
%  c2=CMP.c2;
  c1=-3.9;
  c2=1;
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
for k=ll:-1:1
  sb=A(k,:);
  I=find(~isnan(sb));
  if ~isempty(I), break; end;
end
dmm=zz(k,I);
zb=min(dmm);
fprintf('Bottom zb=%6.1f\n',zb);

%pcolor(xx,ZZ,A); shading interp; % will truncate the bottom in z-levels
pcolor(xx,zz,A); shading flat;
caxis([c1 c2]);
hold on;
colormap(cmp);

%keyboard
if f_layer>0
  nl=size(zz,1);
  for k=1:nl
    z=ZZ(k,:);
    x=xx(k,:);
    plot(x,z,'k-','linewidth',1);
    if (k==plr+1), % bottom interface of the layer
      plot(x,z,'k-','linewidth',1.8);
    end
  end
end;

x=xx(1,:);
set(gca,'Color',[0 0 0],'tickdir','out');
set(gca,'xlim',[min(min(xx)) max(max(xx))],...
	'ylim',[1.12*zb 0],...
	'xtick',[0:250:max(max(xx))],...
	'ytick',[-4500:250:0]);
%set(gca,'fontsize',16);

%  lt0=mean(LAT(j1,i1:i2));
%text(xx(8),200,title_str,'Fontsize',16);
title(title_str,'Interpreter','none');


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