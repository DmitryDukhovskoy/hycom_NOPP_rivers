function sub_plot_grsctS(ip,LL,ZM,Ssm,xl1,xl2,yl1,yl2,Hb,hcb,POS);
% Plot S across Greenland shelf sections
%%
% ip - section #
%
%s0=25;  
s0=30;  
Ssm(Ssm<s0)=s0;
%cf=1e-7;
% lp=8;
%lSsm=((Ssm-s0).^lp)*cf; % expand high S values
%			 %  lSsm=((Ssm-s0).^8)*1e-5; % expand high S values
%cf=1e-8;
cf=1e-5;
lp=9;
lSsm=((Ssm-s0).^lp)*cf; % expand high S values
			 %  lSsm=((Ssm-s0).^8)*1e-5; % expand high S values
c1=0;
c2=18;
			 
%  keyboard
CMP=colormap_WB(200,c1,(c2-c1)/2);
cmp1=CMP.colormap;
cmp1=flipud(cmp1); 
CMP=colormap_WOR(200,(c2-c1)/2,c2);
cmp2=CMP.colormap;
  
cmp=[cmp1;cmp2];
cmp=smooth_colormap(cmp,9);
  
ps=POS(ip,:);
axes('Position',ps);
hold on;
pcolor(LL,ZM,lSsm); shading interp;
caxis([c1 c2]);
colormap(cmp);

hbx=[LL(1);LL;LL(end)];
hby=[-5000;Hb;-5000];
fill(hbx,hby,[0 0 0]);

set(gca,'tickdir','out',...
	'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[0:20:xl2],...
	'ytick',[-4000:200:0],...
	'Fontsize',14);
xlabel('Distance, km');

tcks=[c1:3:c2];
ss=(tcks*1/cf).^(1/lp)+s0;
for kk=1:length(tcks)
  ltck{kk}=sprintf('%3.1f',ss(kk));
end


if ~isempty(hcb);
  hb=colorbar;
  set(hb,'Position',hcb,...
	 'Ticks',tcks,...
	 'TickLabels',ltck,...
	 'Ticklength',0.025,...
	 'Fontsize',14);
  
end

return