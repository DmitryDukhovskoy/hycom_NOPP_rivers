function sub_plot_xsectZZ(nf,xx,zz,A,title_str,fld0,xdr,Hs);
% plot vertical section of tracer, density, etc
% interpolated into regular vertical grid ZZ
%
ZZ=zz;
plr=20; % highlight this interface in HYCOM

SLNB =[28;30;32;32.5;33;33.5;34;34.5;35;35.5];
% Colormap:
hcb=[];
ISL=[];
switch(fld0)
 case('salin');
    na=200;
    
    i1=max(find(SLNB<=min(min(A))));
    if isempty(i1), i1=1; end;
    i2=max(find(SLNB<=max(max(A))));
    if isempty(i2), i2=length(SLNB); end;
    c1=SLNB(i1);
    c2=SLNB(i2);
    cmp = colormap(parula(na));
%    CMP = colormap_sclr2(na,c1,c2);
%    cmp = CMP.colormap;
%    cl1=colormap_gray(na);
    cl2=colormap_purple(100);
    for k=1:3;
      cl2(k,:)=[1 1 1];
    end
%    cl3=colormap_blue(na);
%    cl4=colormap_green(na);
%    cl5=colormap_yellow(na);
%    cl6=colormap_red(na);
%
%    cmp=[cl1;cl2;cl3;cl4;cl5;cl6];
    cmp=[cl2;cmp];
    cmp=smooth_colormap(cmp,20);
    nint=length(cmp);
    cr1=28;
    cr2=35;
    dcr=0.2;
    cr0=33.5;
    
%    ISL=[24:1:32,33:0.1:33.4,33.5:0.05:35];
    if min(min(A))<30
      ISL=[24:1:32,33:0.1:33.4,33.5:0.05:35];
    else
      ISL=[30:0.5:33.5,34.:0.05:35];
    end
    
    
% Use log scale for S to expand
% high S values
%    s0=28;  
%    Ssm=A;
%    Ssm(Ssm<s0)=s0;
%cf=1e-7;
% lp=8;
%lSsm=((Ssm-s0).^lp)*cf; % expand high S values
%			 %  lSsm=((Ssm-s0).^8)*1e-5; % expand high S values
%cf=1e-8;
%    cf=1e-5;
%    lp=9;
%    A=((Ssm-s0).^lp)*cf; % expand high S values
%			 %  lSsm=((Ssm-s0).^8)*1e-5; % expand high S values
%    c1=0;
%    c2=18;
%    
%    CMP=colormap_WB(200,c1,(c2-c1)/2);
%    cmp1=CMP.colormap;
%    cmp1=flipud(cmp1); 
%    CMP=colormap_WOR(200,(c2-c1)/2,c2);
%    cmp2=CMP.colormap;
%
%    cmp=[cmp1;cmp2];
%    tcks=[c1:3:c2];
%    
%    ss=(tcks*1/cf).^(1/lp)+s0;
%    for kk=1:length(tcks)
%      ltck{kk}=sprintf('%3.1f',ss(kk));
%    end
%    
%   hcb=[0.92 0.32 0.01 0.55]; 
    
 case('temp');
    na=240;
    c1=-2;
    c2=4;
    CMP = colormap_sclr2(na,c1,c2);
    cmp = CMP.colormap;
    cr1=-2;
    cr2=10;
    dcr=0.5;
    cr0=0;
    
    nint=length(cmp);
 case('dens');
    na=200;
    cmp = colormap(parula(na));
    nint=length(cmp);
    c1=25.5;
    c2=27.5;
%    CMP = colormap_sclr2(na,c1,c2);
%    cmp = CMP.colormap;
    cr1=20;
    cr2=28;
    dcr=0.2;
    cr0=26;
    nint=length(cmp);
    
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
  c1=-4;
  c2=3;
  cr1=c1;
  cr2=c2;
  cr0=0;
  dcr=0.5;
  
 case('brunt'); % Brunt-V. freq., log scale
  CMP = sub_conc_clrmp02(-8,-4);
  cmp=CMP.colormap;
  cnt=CMP.intervals;
  nint=length(cmp);
  c1=CMP.c1;
  c2=CMP.c2;
  
end;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals


fvsoff=logical(0);
if nf<0
  nf=abs(nf);
  fvsoff=logical(1);
end;

figure(nf); clf;
if fvsoff, set(gcf,'Visible','off'); end;
% 
axes('position',[0.1 0.3 0.8 0.6]);
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
dmm=zz(k);
zb=min(Hs);
fprintf('Bottom zb=%6.1f\n',zb);
%keyboard
pcolor(xx,ZZ,A); shading interp;
%pcolor(xx,zz,A); shading flat;
caxis([c1 c2]);
hold on;
colormap(cmp);

%hold on;
%plot(xx,Hs,'r');
%keyboard
if isempty(ISL)
  contour(xx,ZZ,A,[cr1:dcr:cr2],'Color',[0.3 0.3 0.3]);
  [c,h]=contour(xx,ZZ,A,[cr0:dcr:cr2],'Color',[0.3 0.3 0.3]);
else
  [c,h]=contour(xx,ZZ,A,ISL,'Color',[0.3 0.3 0.3]);
end
%keyboard
clabel(c,h,'Fontsize',8,'Labelspacing',200);
%contour(xx,ZZ,A,[cr0 cr0],'k','Linewidth',1.6);

%x=xx(1,:);
set(gca,'Color',[0.8 0.8 0.8],...
	'tickdir','out',...
        'xlim',[min(min(xx)) max(max(xx))],...
	'ylim',[1.02*zb 0]);
if xdr<0
  set(gca,'xdir','reverse');
end
if abs(zb)>400
  set(gca,'xtick',[0:20:max(xx)],...
	'ytick',[-2000:100:0]);
else
  set(gca,'xtick',[0:20:max(xx)],...
	'ytick',[-500:50:0]);
end
set(gca,'fontsize',14);

xx=xx(:);
xhp=[xx(1); xx; xx(end)];
yhp=[min(Hs)-1000; Hs; min(Hs)-1000];
patch(xhp,yhp,[0 0 0]);

%  lt0=mean(LAT(j1,i1:i2));
%text(xx(8),200,title_str,'Fontsize',16);
title(title_str,'Fontsize',12);

if ~isempty(hcb)

  hb=colorbar;
  set(hb,'Position',hcb,...
	 'Ticks',tcks,...
	 'TickLabels',ltck,...
	 'Ticklength',0.025);
else
  hght=[];
  lngth=[];
  mint=20;
  mbx=20;
  fsz=12;
  bxc='k';
  posc=[0.91 0.3 0.6 0.03];
  aend=0;
  [az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);
end

return