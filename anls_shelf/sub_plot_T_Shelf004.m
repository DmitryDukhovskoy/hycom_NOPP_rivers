function sub_plot_T_Shelf004(nf,TSZ,Tsm,btx,stl,Sc0,nsct,...
			     xl1,xl2,yl1,yl2,Tcntr,...
			     c1,c2,dxx);
% Plot Summer and winter T fields (months: js1:js2)
% sections
ZZ=TSZ.ZZlevels;
ZM=TSZ.ZM;
nzm=length(ZM);
nzz=length(ZZ);

% For plotting add an extra layer 
%Un(nzz,:)=Un(nzz-1,:);

dw=0.4;
dh=0.38;
%POS=[0.07 0.56 dw dh; ...
%     0.53 0.56 dw dh; ...
%     0.07 0.08 dw dh; ...
%     0.53 0.08 dw dh];

%Scntr=[30:0.1:34.8,34.9,34.95,35];

figure(nf); clf;
  
Tsm = sub_fill_bottom_nans(Tsm);
% Filter to get rid of z-level disc.
WT=[0.3, 0.3, 0.25;...
    0.3, 0.3, 0.3;...
    0.3, 0.3, 0.3];
WT=WT./(sum(sum(WT)));
[nw,mw]=size(WT);
di=floor(nw/2);
dj=di;

[mm,nn]=size(Tsm);
for it=1:6
dmm=Tsm;
for ii=2:nn-1
  for jj=2:mm-1
    i1=ii-di;
    i2=ii+di;
    j1=jj-dj;
    j2=jj+dj;
    i1=max([1,i1]);
    i2=min([nn,i2]);
    j1=max([1,j1]);
    j2=min([mm,j2]);
    A=Tsm(j1:j2,i1:i2);
    df=sum(sum(A.*WT));
    dmm(jj,ii)=df;
  end
end
Tsm=dmm;
end


LL=TSZ(nsct).Dist_origin*1e-3; % km
nm=TSZ(nsct).Name;
Hb=TSZ(nsct).Hbottom;

%keyboard  

%  hcb=[];
dw=0.8;
dh=0.5;
POS=[0.1 0.3 dw dh];
hcb=[0.92 0.35 0.006 0.4]; 

ip=1;   
%c1=25;
%c2=35;
CMP=create_colormap7(400,c1,0);
cmp=CMP.colormap;

ps=POS(ip,:);
axes('Position',ps);
hold on;
pcolor(LL,ZM,Tsm); shading interp;
caxis([c1 c2]);
colormap(cmp);

hbx=[LL(1);LL;LL(end)];
hby=[-5000;Hb;-5000];
fill(hbx,hby,[0 0 0]);

set(gca,'tickdir','out',...
	'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[0:dxx:xl2],...
	'ytick',[-800:50:0],...
	'Fontsize',14);
xlabel('Distance, km');

tcks=[c1:1:c2];
%ss=(tcks*1/cf).^(1/lp)+s0;
for kk=1:length(tcks)
  ltck{kk}=sprintf('%3.1f',tcks(kk));
end


if ~isempty(hcb);
  hb=colorbar;
  set(hb,'Position',hcb,...
	 'Ticks',tcks,...
	 'TickLabels',ltck,...
	 'Ticklength',0.025,...
	 'Fontsize',14);
  
end



contour(LL,ZM,Tsm,Tcntr,'k');
contour(LL,ZM,Tsm,[Sc0 Sc0],'k','Linewidth',1.6);

stl=sprintf('%s, %s',stl,nm);
title(stl,'Interpreter','none');




return