function sub_check_sect(nf,HH,SCT)
% Plot boxes, segments etc.
%
%
btx='sub_check_sect';
CLR=[0 0.4 0.8; ...
     0.8 0.4 0; ...
     1 0.2 0; ...
     0 1 0; ...
     0.8 0 0.6; ...
     0.8 1 0];

figure(nf); clf;
contour(HH,[0 0],'k');
hold on;
%contour(HH,[-5000:500:-100],'Color',[0.9 0.9 0.9]);

nsct=length(SCT);
for isct=1:nsct
		iGTs=SCT(isct).sgmX;
		jGTs=SCT(isct).sgmY;
		IJp=SCT(isct).IJ_indx;
		IJadj=SCT(isct).IJ_adj_indx;
		Nrm=SCT(isct).Norm;
		nn=length(iGTs);

		clr=CLR(isct,:);
		for j=1:nn
				x1=iGTs(j,1);
				x2=iGTs(j,2);
				y1=jGTs(j,1);
				y2=jGTs(j,2);
				plot([x1 x2],[y1 y2],'.-','Color',clr);
				xn0=0.5*(x1+x2);
				yn0=0.5*(y1+y2);
				xn1=xn0+Nrm(j,1);
				yn1=yn0+Nrm(j,2);

				ppx=IJp(j,1);
				ppy=IJp(j,2);
				apx=IJadj(j,1);
				apy=IJadj(j,2);
				plot([xn0 xn1],[yn0 yn1],'-','Color',[0.8 0.3 0]);
				plot(ppx,ppy,'.','Markersize',18,'Color',[0.1 0.8 0.9]);
				plot(apx,apy,'*','Color',[0.8 0.1 0.9]);
  end
end

axis('equal');
set(gca,'xlim',[300 1150],...
       	'ylim',[50 800]);

bottom_text(btx,'pwd',1);


return
