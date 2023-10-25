% plot T/S diagram
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup;

close all
clear

% Plot sigmas:
SGR=[32:0.05:35.9];
TGR=[-1.8:0.1:11];
ns=length(SGR);
nt=length(TGR);
clear sgm0
for ik=1:nt;  % 
  rho=sw_dens0(SGR,ones(1,ns)*TGR(ik))-1000;
  sgm0(ik,1:length(rho))=rho;
end;

figure(1); clf;
contour(SGR,TGR,sgm0,[10:0.2:30],'b','linewidth',1.5);
hold on
contour(SGR,TGR,sgm0,[27 27],'r','linewidth',1.5);
%contour(SGR,TGR,sgm0,[21 21],'r');
axis('equal');
set(gca,'xlim',[33.5 35.8],'ylim',[2.5 7.],'tickdir','out',...
	'xgrid','on','ygrid','on');
set(gca,'xtick',[32:0.25:35.8],'ytick',[-1.5:0.5:11]);

pthf = '/Net/ocean/ddmitry/vector_winds/fig_EN4/';
  fignm=sprintf('%sTSdgr',pthf);
  fprintf('Saving figure: %s \n',fignm);
  print('-dpng','-f1','-r0',fignm);
