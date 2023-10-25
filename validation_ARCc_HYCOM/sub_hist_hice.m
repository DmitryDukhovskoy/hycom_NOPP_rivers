function sub_hist_hice(hice,hbin,nf);
% Plot histogram of ice thicknesses

I=find(hice>0);
nhh = hist(hice(I),hbin);
dh=hbin(2)-hbin(1);
nhh=nhh/(sum(nhh)*dh);

figure(nf); clf;
axes('Position',[0.08 0.6 0.8 0.3]);
hold on
hb = bar(hbin,nhh,0.95);
set(hb,'FaceColor',[0 0.7 0.8]);

set(gca,'tickdir','out',...
        'xlim',[-0.1 max(hbin)+dh],...
        'ylim',[0 1.1*max(nhh)],...
        'xtick',[0:dh:20],...
        'fontsize',14,...
        'box','on');


return

