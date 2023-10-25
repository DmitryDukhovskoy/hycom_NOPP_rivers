function sub_plot_S_Shelf004(nf,TSZ,Ssm,btx,stl,Sc0,nsct,...
			     xl1,xl2,yl1,yl2,Scntr,...
			     c1,c2,dxx);
% Plot Summer and winter S fields (months: js1:js2)
% sections
% Sc0 - salinity defining the plume front
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

% keyboard
% 
Ssm = sub_fill_bottom_nans(Ssm);
% Filter to get rid of z-level disc.
WT=[0.3, 0.3, 0.25;...
    0.3, 0.3, 0.3;...
    0.3, 0.3, 0.3];
WT=WT./(sum(sum(WT)));
[nw,mw]=size(WT);
di=floor(nw/2);
dj=di;

[mm,nn]=size(Ssm);
for it=1:6
dmm=Ssm;
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
    A=Ssm(j1:j2,i1:i2);
    df=sum(sum(A.*WT));
    dmm(jj,ii)=df;
  end
end
Ssm=dmm;
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
CMP=colormap_WB(200,c1,(c2-c1)/2);
cmp1=CMP.colormap;
cmp1=flipud(cmp1); 
CMP=colormap_WOR(200,(c2-c1)/2,c2);
cmp2=CMP.colormap;
cmp=[cmp1;cmp2];
cmp=smooth_colormap(cmp,9);
%CMP=create_colormap_freshwater();

ps=POS(ip,:);
axes('Position',ps);
hold on;
pcolor(LL,ZM,Ssm); shading interp;
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



contour(LL,ZM,Ssm,Scntr,'k');
contour(LL,ZM,Ssm,[Sc0 Sc0],'k','Linewidth',1.6);

stl=sprintf('%s, %s',stl,nm);
title(stl,'Interpreter','none');
  


return