function sub_plot_TxsctZ(nfg,Hb,F,ZZ,X,dnmb,stl,cntr,cc0);
% Plot T xsection
% on Z grid
%
F=[F(1,:);F]; % plotting pcolor
[ma,na]=size(F);
if size(ZZ,1)<ma
  ZZ=[ZZ(1,:)*0;ZZ];
end


c1=-1.8;
c2=4.8;
dc = 0.6;
nint = 200;
%  CMP = create_colormap5(nint,c1,c2);
%  cmp = CMP.colormap;
%  cnt = CMP.intervals;
%
%  nav = 15;
%  cmp = smooth_colormap(cmp,nav);
cmp = colormap(jet(200));


for k=1:ma
  XX(k,:)=X;
end

%keyboard

figure(nfg); clf;
AX=axes('Position',[0.1 0.3 0.75 0.6]);
pcolor(X,ZZ,F); shading interp
hold on
%pcolor(XL,ZZ,F); shading flat
colormap(cmp);
caxis([c1 c2]);
%keyboard

if ~isempty(cntr),
  cc1=min(cntr);
  cc2=max(cntr);
  dcc=cntr(2)-cntr(1);
  contour(XX,ZZ,F,[cc1:dcc:cc2],'Color',[1 1 1]);
  contour(XX,ZZ,F,[cc0 cc0],'Linewidth',1.6,'Color',[1 1 1]);
end  

hc = colorbar;
set(hc,'Position',[0.89 0.3 0.022 0.6],...
       'Ticks',[c1:dc:c2],...
       'TickLength',0.042,...
       'Fontsize',12);

%keyboard
%stt=sprintf('%s, Volume Flux, m3/s',nm);
%title(stt);
xl1=X(1);
xl2=X(end);
if xl1>xl2
  xl1=xl2;
  xl2=XL(1);
end

set(gca,'Color',[0 0 0],...
	'tickdir','out',...
	'xlim',[xl1 xl2],...
	'xtick',[-20:5:20],...
	'ylim',[min(min(ZZ))-1 0],...
	'ytick',[-5000:500:0],...
	'Fontsize',12);

title(stl);

return