function sub_check_boxes(nf,HH,GC,BOX)
% Plot boxes, segments etc.
%
%
btx='sub_check_boxes';
CLR=[0 0.4 0.8; ...
     0.8 0.4 0; ...
     1 0.2 0; ...
     0 1 0; ...
     0.8 0 0.6; ...
     0.8 1 0];

figure(nf); clf;
contour(HH,[0 0],'k');
hold on;
contour(HH,[-5000:500:-100],'Color',[0.9 0.9 0.9]);

Igs=GC.sgmX;
Jgs=GC.sgmY;
NrmG=GC.Norm;
ijPg=GC.IJ_indx;
ijADJg=GC.IJ_adj_indx;

ngr=length(GC.Ictr);
nbx=length(BOX);
for ibx=1:nbx
		for isct=1:2
				iGTs=BOX(ibx).S(isct).sgmX;
				jGTs=BOX(ibx).S(isct).sgmY;
				Igtx=BOX(ibx).S(isct).Intrcp_Gate;
				Igt0=BOX(ibx).S(isct).Start_Gate;
				Igrx(isct)=BOX(ibx).S(isct).Intrcp_Greenl;
				IJp=BOX(ibx).S(isct).IJ_indx;
				IJadj=BOX(ibx).S(isct).adjIndx_I1J1;
				Nrm=BOX(ibx).S(isct).Norm;
				plot(iGTs(:,1),jGTs(:,1),'-','Color',[0.5 0.5 0.5]);

				if Igt0>Igtx
						ii1=Igtx;
						ii2=Igt0;
				else
						ii1=Igt0;
						ii2=Igtx;
				end

				clr=CLR(ibx,:);
				for j=ii1:ii2
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
						if isct==1
								plot([xn0 xn1],[yn0 yn1],'-','Color',[0.8 0.3 0]);
								plot(ppx,ppy,'.','Markersize',18,'Color',[0.1 0.8 0.9]);
								plot(apx,apy,'*','Color',[0.1 0.8 0.9]);
						else
								plot([xn0 xn1],[yn0 yn1],'-','Color',[0.2 0.9 0]);
								plot(ppx,ppy,'.','Markersize',12,'Color',[0.2 0.3 0.]);
								plot(apx,apy,'o','Color',[0.2 0.3 0.]);
						end
				end
		end

%
% For 1st segment special care 
% of the cut of the contour
  ix1=[];
  ix2=[];
  if Igrx(1)>Igrx(2) & ibx==1
    di=ngr-Igrx(1);
    ix1=-di;
    ix2=Igrx(2);
  end

  if isempty(ix1);
				if Igrx(1)<Igrx(2)
						ix1=Igrx(1);
						ix2=Igrx(2);
				else
						ix1=Igrx(2);
						ix2=Igrx(1);
				end
  end
    
		for j=ix1:ix2
    j1=j;
    if j<=0;
      j1=ngr+j;
    end
				x1=Igs(j1,1);
				x2=Igs(j1,2);
				y1=Jgs(j1,1);
				y2=Jgs(j1,2);

				xn0=0.5*(x1+x2);
				yn0=0.5*(y1+y2);
				xn1=xn0+NrmG(j1,1);
				yn1=yn0+NrmG(j1,2);

				plot([x1 x2],[y1 y2],'-','Linewidth',2,'Color',clr);
    plot([xn0 xn1],[yn0 yn1],'-','Color',[0 0.8 1]);

		end

end

axis('equal');
set(gca,'xlim',[450 1100],...
	'ylim',[350 1200]);

bottom_text(btx,'pwd',1);


return
