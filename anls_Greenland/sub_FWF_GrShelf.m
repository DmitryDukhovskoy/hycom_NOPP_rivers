function sub_FWF_GrShelf(FWF,UZGR,nf,stl)
% Caclulate FWF across sections on Gr Shelf

nsct=length('UZGR');

dw=0.4;
dh=0.18;
POS=[0.07 0.80 dw dh; ...
     0.07 0.55 dw dh; ...
     0.07 0.30 dw dh; ...
     0.07 0.05 dw dh];
POS2=[0.53 0.80 dw dh; ...
     0.53 0.55 dw dh; ...
     0.53 0.30 dw dh; ...
     0.53 0.05 dw dh];

CLR1=[0 0.5 0.8; ...
     0 0.7 0.3; ...
     0.8 0 0.3; ...
     0.7 0. 0.9];

CLR2=[0 0.6 1; ...
     0 1 0.6; ...
     1 0 0.4; ...
     0.9 0. 0.9];

mnth=[1.5:12.5];
figure(nf); clf;

% Plot fluxes on GrShelf bounded by the contour 
for ip=1:nsct
  sname=UZGR(ip).Name;
  fwf=FWF(ip).FWF_GrSh*1e-3; % mSv
  pL=prctile(fwf,25);
  pU=prctile(fwf,75);
  mn=mean(fwf,1);
  dmm=fwf(:);
  amn=nanmean(dmm);
  stdv=std(dmm);
  
  
  pos=POS(ip,:);
  axes('position',pos);
  clr=CLR1(ip,:);
  plot(mnth,mn,'Linewidth',2.5,'Color',clr);
  hold on;
  for ik=1:12
    plot([ik+0.5 ik+0.5],[pL(ik) pU(ik)],'k--','Color',clr);
  end
  
  stxt=sprintf('%4.1f+/-%4.1f',amn,stdv);
  text(1,0.8*min(pL),stxt);
  set(gca,'tickdir','out',...
	  'xlim',[1 12.8],...
	  'xtick',[1:12],...
	  'ylim',[1.05*min(pL) max([0, max(pU)])]);
  stlL=sprintf('%s GrShCntr: %s',stl,sname);
  title(stlL,'Interpreter','none');
  
end

% Whole Section: Use this for
% Fluxes in Davis, Denmark Straits
for ip=1:nsct
  sname=UZGR(ip).Name;
  fwf=FWF(ip).FWF_sect*1e-3; % mSv
  pL=prctile(fwf,25);
  pU=prctile(fwf,75);
  mn=mean(fwf,1);
  dmm=fwf(:);
  amn=nanmean(dmm);
  stdv=std(dmm);
  
  pos=POS2(ip,:);
  axes('position',pos);
  clr=CLR1(ip,:);
  plot(mnth,mn,'Linewidth',2.5,'Color',clr);
  hold on;
  for ik=1:12
    plot([ik+0.5 ik+0.5],[pL(ik) pU(ik)],'k--','Color',clr);
  end
  
  stxt=sprintf('%4.1f+/-%4.1f',amn,stdv);
  text(1,0.8*min(pL),stxt);
  set(gca,'tickdir','out',...
	  'xlim',[1 12.8],...
	  'xtick',[1:12],...
	  'ylim',[1.05*min(pL) max([0, max(pU)])]);
  stlL=sprintf('%s WhSect: %s',stl,sname);
  title(stlL,'Interpreter','none');
end




return