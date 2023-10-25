function sub_plot_Month_flxsct1d(AA,dmo,xsct,fnmb,yt1,yt2,dy,stl);
% Plot monthly fluxes across section
% all months on 1 diagram
cmp = [0,0,1;
       0,0.5,1;
       0,1,1;
       0,1,0.5;
       0,1,0;
       0.5,1,0;
       1,1,0;
       1,0.5,0;
       1,0,0;
       1,0,0.5;
       1,0,1;
       0.5,0,1];

figure(fnmb); clf;
axes('Position',[0.08 0.5 0.88 0.4]);
hold on;
for im=1:dmo:12
  a=AA(im,:);
  clr = cmp(im,:);
  plot(xsct,a,'Color',clr,'linewidth',1.6);
end

yl1 = 1.05*min(min(AA));
yl2 = 1.05*max(max(AA));


set(gca,'tickdir','out',...
	'xlim',[0 max(xsct)],...
	'ylim',[yl1 yl2], ...
	'xtick',[0:500:10000],...
	'ytick',[yt1:dy:yt2],...
	'xgrid','on',...
	'ygrid','on');
	
title(stl);

axes('Position',[0.08 0.2 0.2 0.25]);
hold on;
cc=0;
for im=1:dmo:12
  clr = cmp(im,:);
  cc=cc+1;
  plot([0 1],[-cc -cc],'Color',clr,'linewidth',1.6);
  str=sprintf('%2.2i',im);
  text(1.2, -cc, str,'Fontsize',11);
end
set(gca,'xlim',[0 4],...
	'ylim',[-cc-1 0],...
	'box','on',...
	'visible','off');




return