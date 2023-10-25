function sub_setgca_iceA_all(fgn,EXPT,CLR,TM,tm0,ymin,ymax,lpos);


DV = datevec(TM);
dm = diff(DV(:,2));
tmm = TM-tm0;

yr1 = min(DV(:,1));
yr2 = max(DV(:,1));

%keyboard

clear dtm d15
icc=0;
for iyr=yr1:yr2
	for imo=1:12
		icc=icc+1;
		dtm(icc,1) = datenum(iyr,imo,1);
	end
end
dtm=dtm-tm0;

icc=0;
for iyr=yr1:yr2
	for imo=1:12
		icc=icc+1;
		d15(icc,1) = datenum(iyr,imo,15);
		smo{imo}=sprintf('%2.2i',imo);
	end;
end
d15=d15-tm0;

%dtm(dtm<0)=nan;
%d15(d15<0)=nan;

for ii=1:length(dtm)
	x1=dtm(ii);
  if isnan(x1); continue; end;
	plot([x1 x1],[ymin ymax],'--','Color',[0.6 0.6 0.6]);
end

yl1=min([ymin,0]);
if yl1<0, yl1=1.05*yl1; end
yl2=1.05*ymax;
ytxt = yl1-0.1*yl2;
	
for ii=1:length(dtm)
	x1=d15(ii);
	plot([x1 x1],[yl1 yl2],':','Color',[0.8 0.8 0.8]);
	stl = smo{ii};
	text(x1, ytxt,stl,'Fontsize',14);
end

set(gca,'tickdir','out',...
				'xlim',[tmm(1) max(tmm)],...
				'xtick',[ ],...
				'ylim',[yl1 yl2],...
				'ygrid','on',...
				'Fontsize',14);

%keyboard
axes('Position',lpos);
hold on;
yy0=length(EXPT)+1;
xx1=0.05;
xx2=xx1+0.2;
xx3=xx2+0.1;

for ixx=1:length(EXPT);
  nm = sprintf('%s %3.2f_%s',EXPT(ixx).cice,EXPT(ixx).res,EXPT(ixx).cice_opt);
  clr = CLR(ixx,:);
  iy = yy0-ixx+1;
  plot([xx1 xx2],[iy iy],'-','Color',clr,'Linewidth',2);
  text(xx3,iy,nm,'Fontsize',12,'Interpreter','none');
end

set(gca,'xlim',[0 5],...
        'ylim',[1 yy0+1],...
        'visible','off')




return

