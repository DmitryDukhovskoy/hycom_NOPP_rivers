% Plot quick figures
% Using vol_intgr_trcr_00[84].m
% To analyze concentration difference
% in the simulations
% first run vol_intgr_trcr_00[84].m
% in different matlab sessions
% IT will create figures with conc. fields 
% zoom in onto some region defined by i1,i2,j1,j2
% Will need to jump between the 0.08 and 0.04 
% matlab windows by C&P into prompts

% 0.08:
i1=630; i2=660; j1=460; j2=510;
%i1=540; i2=555; j1=460; j2=500;
regn = 'ARCc0.08';

% 0.04
i1=630*2; i2=660*2; j1=460*2; j2=510*2;
%i1=540*2; i2=555*2; j1=460*2; j2=500*2;
regn = 'ARCc0.04';
 
% 
A=lTr(j1:j2,i1:i2);
ar=Acell(j1:j2,i1:i2);
hh=HH(j1:j2,i1:i2);
M=nansum(nansum(A.*ar))*1e-9; % roughly gives the total mass in the region

A(hh<-500)=nan;
figure(3); clf;
pcolor(A); shading flat
hold on
caxis([4 7]);
colorbar
contour(hh,[-2000:200:0],'k');

stl = sprintf('%s, IntgrM=%6.1f kg',regn,M);
title(stl);
% title('ARCc0.08');
% title('ARCc0.04');

I=find(A>1e-3);
dx=0.5;
xx=[0:dx:8];
figure(4);
hb=hist(A(I),xx);
% Normalize:
ahb = nansum(hb*dx);
nhb = hb/ahb;
bar(xx,nhb);
title(stl);
set(gca,'tickdir','out',...
	'xtick',[0:8],...
	'ygrid','on',...
	'ylim',[0 0.75]);







