function [data_min,data_max,time_min,time_max] = sub_plot_iceA_all(fld,A,TM,tm0,clr);

fprintf('Plotting %s\n',fld);
stl = sprintf('amed = A.%s_mn;',fld);
eval(stl);
%stl = sprintf('ap1 = A.%s_iqr(:,1);',fld);
%eval(stl);
%stl = sprintf('ap2 = A.%s_iqr(:,2);',fld);
%eval(stl);

%yl2 = max(ap2);
%yl1 = min([0,min(ap1)]);

DV = datevec(TM);
dm = diff(DV(:,2));
tmm = TM-tm0;


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

%if fstart == 1
%	figure(fgn); clf;
%	set(gcf,'Position',[1627         844         872         474]);
%	axes('Position',[0.09 0.5 0.83 0.42]);
%end
plot(tmm,amed,'-','Linewidth',2.5,'Color',clr);
hold on;
%plot(tmm,ap1,'--','Linewidth',1.6,'Color',[0.4 0.4 0.4]);
%plot(tmm,ap2,'--','Linewidth',1.6,'Color',[0.4 0.4 0.4]);

time_min = min(tmm);
time_max = max(tmm);
data_min = min(amed);
data_max = max(amed);


return

