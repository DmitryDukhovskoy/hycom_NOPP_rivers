function sub_plot_grsctT(c1,c2,ip,LL,ZM,Fld,xl1,xl2,yl1,yl2,Hb,hcb,POS);
% Plot T across Greenland shelf sections
%%
% ip - section #
%
dw=0.4;
dh=0.38;
%POS=[0.07 0.56 dw dh; ...
%     0.53 0.56 dw dh; ...
%     0.07 0.08 dw dh; ...
%     0.53 0.08 dw dh];

%keyboard

CMP=create_colormap7(400,c1,0);
cmp=CMP.colormap;

ps=POS(ip,:);
axes('Position',ps);
hold on;
pcolor(LL,ZM,Fld); shading interp;
caxis([c1 c2]);
colormap(cmp);

hbx=[LL(1);LL;LL(end)];
hby=[-5000;Hb;-5000];
fill(hbx,hby,[0 0 0]);

set(gca,'tickdir','out',...
	'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[0:20:xl2],...
	'ytick',[-4000:200:0]);
xlabel('Distance, km');
if ~isempty(hcb);
  hb=colorbar;
  set(hb,'Position',hcb,...
	 'Ticklength',0.025,...
	 'Fontsize',14);
  
end

return