function sub_plot_hflxZcntr(FF,ZZ,xx,nf,title_str,Hs,c1,c2);
% plot vertical section of tracer, density, etc
% interpolated into regular vertical grid ZZ
%
% Colormap:
na = 100;
cl2 = colormap_orange(na);
cl1 = colormap_blue(na);
for ik=1:10
  cl1(ik,:) = [1 1 1];
  cl2(ik,:) = [1 1 1];
end
cl1 = flipud(cl1);
cmp = [cl1;cl2];
cmp = smooth_colormap(cmp,10);

%c1   = -0.6;
%c2   =  0.6;
nint = length(cmp);
cnt  = (c1:(c2-c1)/nint:c2);  % construct array of intervals


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
%FF=lTr;
% Find bottom:
[ll,nn]=size(FF);
for k=ll:-1:1
  sb=FF(k,:);
  I=find(~isnan(sb));
  if ~isempty(I), break; end;
end
zb=min(min(ZZ));
fprintf('Bottom zb=%6.1f\n',zb);
%keyboard
pcolor(xx,ZZ,FF); shading interp;
%pcolor(xx,zz,FF); shading flat;
caxis([c1 c2]);
hold on;
colormap(cmp);

hold on;
%plot(xx,Hs,'r');
%keyboard
%contour(xx,ZZ,FF,[cr1:dcr:cr2],'Color',[0.3 0.3 0.3]);
%contour(xx,ZZ,FF,[cr1:1:cr2],'k');

%x=xx(1,:);
set(gca,'Color',[0.8 0.8 0.8],...
	'tickdir','out',...
        'xlim',[min(min(xx)) max(max(xx))],...
	'ylim',[1.002*zb 0]);
set(gca,'xtick',[0:500:max(max(xx))],...
	'ytick',[-3000:250:0]);

set(gca,'fontsize',12);
%keyboard

x1x=xx(1,:);
x1x=x1x(:);
xhp=[x1x(1); x1x; x1x(end)];
Hs=Hs(:);
yhp=[min(Hs)-1000; Hs; min(Hs)-1000];
patch(xhp,yhp,[0 0 0]);

%  lt0=mean(LAT(j1,i1:i2));
%text(xx(8),200,title_str,'Fontsize',16);
title(title_str,'Fontsize',12);

chb = colorbar;
set(chb,'Position',[0.92 0.12 0.012 0.8],...
	'TickLength',0.03,...
	'FontSize',12);

return
