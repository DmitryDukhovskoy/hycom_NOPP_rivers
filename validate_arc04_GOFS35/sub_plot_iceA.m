function sub_plot_iceA(fld,A,TM,sttl,fgn);

fprintf('Plotting %s\n',fld);
stl = sprintf('amed = A.%s_mn;',fld);
eval(stl);
stl = sprintf('ap1 = A.%s_iqr(:,1);',fld);
eval(stl);
stl = sprintf('ap2 = A.%s_iqr(:,2);',fld);
eval(stl);

yl2 = max(ap2);
yl1 = min([0,min(ap1)]);

DV = datevec(TM);
dm = diff(DV(:,2));
tmm = TM-TM(1);
iM = find(dm==1);

clear dtm d15
icc=0;
for imo=1:12
  icc=icc+1;
  dtm(icc,1) = datenum(DV(1,1),imo,1);
end
dtm=dtm-TM(1);

icc=0;
for imo=1:12
  icc=icc+1;
  d15(icc,1) = datenum(DV(1,1),imo,15);
  smo{imo}=sprintf('%2.2i',imo);
end;
d15=d15-TM(1);

dtm(dtm<0)=nan;
d15(d15<0)=nan;

%amed
%keyboard

figure(fgn); clf;
set(gcf,'Position',[1627         844         872         474]);
axes('Position',[0.09 0.5 0.83 0.42]);
plot(tmm,amed,'-','Linewidth',2.5);
hold on;
plot(tmm,ap1,'--','Linewidth',1.6,'Color',[0.4 0.4 0.4]);
plot(tmm,ap2,'--','Linewidth',1.6,'Color',[0.4 0.4 0.4]);

for ii=1:length(dtm)
  x1=dtm(ii);
  plot([x1 x1],[yl1 yl2],'--','Color',[0.6 0.6 0.6]);
end

ytxt = yl1-0.1*yl2;
  
for ii=1:length(dtm)
  x1=d15(ii);
  plot([x1 x1],[yl1 yl2],':','Color',[0.8 0.8 0.8]);
  stl = smo{ii};
  text(x1, ytxt,stl,'Fontsize',14);
end

set(gca,'tickdir','out',...
        'xlim',[1 max(tmm)],...
        'xtick',[ ],...
        'ylim',[yl1 yl2],...
        'ygrid','on',...
        'Fontsize',14);

%sttl = sprintf('%s %s',texpt,fld);
title(sttl);
 
%keyboard

return

