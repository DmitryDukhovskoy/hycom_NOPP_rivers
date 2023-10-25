function sub_plot_flxsct1d(AA,xsct,fnmb,yl1,yl2,yt1,yt2,dy);
% Plot fluxes across section

figure(fnmb); clf;
axes('Position',[0.08 0.5 0.88 0.4]);
plot(xsct,AA,'Color',[0 0 0],'linewidth',2.5);
hold on;
if isempty(yl1);
  yl1 = 1.02*min(AA);
  yl2 = 1.02*max(AA);
end

set(gca,'tickdir','out',...
	'xlim',[0 max(xsct)],...
	'ylim',[yl1 yl2], ...
	'xtick',[0:500:10000],...
	'ytick',[yt1:dy:yt2],...
	'xgrid','on',...
	'ygrid','on',...
	'FontSize',14);
	
	


return