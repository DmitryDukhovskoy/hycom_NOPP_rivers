function sub_plot_Vxsct(nfg,Hb,F,ZZ,X,dnmb,stl,cntr,cc0,mbtm);
% Plot V xsection
% 
%
F=[F(1,:);F]; % plotting pcolor
[ma,na]=size(F);
if size(ZZ,1)<ma
  ZZ=[ZZ(1,:)*0;ZZ];
end



c1=-0.3;
c2=0.3;
nint = 200;

cl1=flipud(colormap_blue(100));
cl1(end,:)=[1 1 1];
cl2=colormap_red(100);
cl2(1,:)=[1 1 1];
cmp=[cl1;cl2];
cmp=smooth_colormap(cmp,5);

if size(ZZ,2)>1
  for k=1:ma
    XX(k,:)=X;
  end
else
  XX=X;
end

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

dc=0.1;
hc = colorbar;
set(hc,'Position',[0.89 0.3 0.022 0.6],...
       'Ticks',[c1:dc:c2],...
       'TickLength',0.042,...
       'Fontsize',12);


if mbtm==1
  xbv=[X(1),X,X(end)];
  zbv=[-9000,Hb,-9000];
  patch(xbv,zbv,[0 0 0]);
end


%keyboard
%stt=sprintf('%s, Volume Flux, m3/s',nm);
%title(stt);
xl1=X(1);
xl2=X(end);
if xl1>xl2
  xl1=xl2;
  xl2=XL(1);
end

yl1=round(min(Hb))-10;

set(gca,'Color',[0 0 0],...
	'tickdir','out',...
	'xlim',[xl1 xl2],...
	'xtick',[-20:5:20],...
	'ylim',[min(min(ZZ))-1 0],...
	'ytick',[-5000:500:0],...
	'Fontsize',12);

title(stl);

return