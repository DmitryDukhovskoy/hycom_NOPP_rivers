function sub_plot_dSSH(pthfig,fmat);
% Plot Lapl. SSH in the BG
% from HYCOM simulation
% extracted in calc_dSSH.m
% Laplacian <0 => convexed surface with local maximum 
fprintf('Loading %s\n',fmat);
load(fmat);

cf=1e9;
TM = LAPLE.TM;
dE = LAPLE.d2E_BG*cf;

DV=datevec(TM);
nc=length(TM);
TT=(TM-TM(1))/365.25+DV(1,1);

% Annual means:
% Year is Jan - Dec
cp=0;
clear dEY
for iyr=DV(1,1):DV(end,1),
  cp=cp+1;
  IY=find(DV(:,1)==iyr);
  dEY(cp,1)=mean(dE(IY));
end


% Winter d2E:
IS=find(DV(:,2)>4 & DV(:,2)<10);
IW=find(DV(:,2)<=4 | DV(:,2)>9);

clear TY dEW dES
cp=0;
for iyr=DV(1,1):DV(end,1),
% Winter: Oct - March (previous/current year)
%  IYW=find(DV(:,1)==iyr & (DV(:,2)<4 | DV(:,2)>9));
  IYW=find((DV(:,1)==iyr-1 & DV(:,2)>9) |...
	   (DV(:,1)==iyr & DV(:,2)<4));
  dmm=mean(dE(IYW));
  vmm=mean(dE(IYW));
  cp=cp+1;
  dEW(cp,1)=dmm; % winter curl by years
  TY(cp,1)=iyr;
  
% Summer
  IYS=find(DV(:,1)==iyr & (DV(:,2)>4 & DV(:,2)<10));
  smm=mean(dE(IYS));
  dES(cp,1)=smm;
end



irr=find(TT<1997);
ipp=find(TT>=1997);
mn90=mean(dE(irr)); % annual mean LaplSSH prior 1997
mn00=mean(dE(ipp)); % annueal mean LaplSSH after 1997

iw1=find(TY<1997);
iw2=find(TY>=1997);
wmn90=mean(dE(iw1)); % winter mean d2E prior 1997
wmn00=mean(dE(iw2)); % winter mean d2E >1997


figure(1); clf;
axes('Position',[0.08 0.7 0.85 0.25]);
bb=bar(TY,dEY);
set(bb,'FaceColor',[0. 0.8 0.3],'EdgeColor','none');
set(gca,'tickdir','out',...
	'ylim',[-15 0],...
	'ytick',[-16:2:2],...
	'xlim',[1992.5 2015.5],...
	'xtick',[1993:2016]);
sll=sprintf('ARCc0.08_110, Annual d2E*%1.0e BG',cf);
title(sll,'Interpreter','none');

% Winter
axes('Position',[0.08 0.39 0.85 0.25]);
bb=bar(TY,dEW);
set(bb,'FaceColor',[0. 0.8 0.3],'EdgeColor','none');
set(gca,'tickdir','out',...
	'ylim',[-15 0],...
	'ytick',[-16:2:2],...
	'xlim',[1992.5 2015.5],...
	'xtick',[1993:2016]);
sll=sprintf('ARCc0.08_110, Winter d2E*%1.0e BG',cf);
title(sll,'Interpreter','none');

axes('Position',[0.08 0.07 0.85 0.25]);
bb=bar(TY,dES);
set(bb,'FaceColor',[0. 0.8 0.3],'EdgeColor','none');
set(gca,'tickdir','out',...
	'ylim',[-15 0],...
	'ytick',[-16:2:2],...
	'xlim',[1992.5 2015.5],...
	'xtick',[1993:2016]);
sll=sprintf('ARCc0.08_110, Summer d2E*%1.0e BG',cf);
title(sll,'Interpreter','none');

btx='sub_plot_dSSH.m';
bottom_text(btx,'pwd',1);

%keyboard

return