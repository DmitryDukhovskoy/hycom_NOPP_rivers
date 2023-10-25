function sub_anls_diff(S1,S2,HH,fav,mm,nn,nlr,LRS);
% compare S at different levels
% from 2 experiments
% fav = 0 - need to divide by # of records
%     = 1 - division has been already done

if fav==0
  for ilr=1:nlr
    cnc = S1(ilr).nrec;
    dmm = S1(ilr).Fld;
    S1(ilr).Fld = dmm./cnc;

    cnc = S2(ilr).nrec;
    dmm = S2(ilr).Fld;
    S2(ilr).Fld = dmm./cnc;
  end
end

i1=375;
i2=1075;
j1=207;
j2=634;

Dmm=HH*0;
Dmm(j1:j2,i1:i2)=HH(j1:j2,i1:i2);

cl1=colormap_blue(100);
cl2=colormap_red(100);
for ik=1:12;
  cl1(ik,:)=[1 1 1];
  cl2(ik,:)=[1 1 1];
end
cl1=flipud(cl1);
cmp = [cl1;cl2];
cmp = smooth_colormap(cmp,10);

c1 = -0.1;
c2 = 0.1;


xl1=300;
xl2=1200;
yl1=200;
yl2=1100;

Lmsk = HH*0;
Lmsk(HH<0)=1;
lmp=[0 0 0; 1 1 1];

Idp=find(Dmm<-1000);

for ilr=1:nlr
  F1 = S1(ilr).Fld;
  F2 = S2(ilr).Fld;
  F1(F1==0)=nan;
  F2(F2==0)=nan;
  dF = F2-F1;
  zb2 = LRS(ilr+1);
  dF(HH>zb2)=nan;
  
  j0=313;
  i0=833;
  fprintf('ilr=%i, Max dF=%6.4f, min dF=%6.4f\n',...
	  ilr,max(dF(Idp)),min(dF(Idp)));
%  if ilr==4 & max(abs(dF(Idp)))>0.2,
%    fprintf('Max dF=%6.4f, min dF=%6.4f\n',max(dF(Idp)),min(dF(Idp)));
%    keyboard;
%  end
  
  
  iplt=0;
  if iplt==1
    figure(ilr); clf;
    pcolor(Lmsk); shading flat;
    colormap(lmp);
    freezeColors;
    hold on;

    pcolor(dF); shading flat;
    caxis([c1 c2]);
  %  hold on;
  %  contour(HH,[0 0],'w','linewidth',1.2);
    axis('equal');
    set(gca,'xlim',[xl1 xl2],...
	    'ylim',[yl1 yl2],...
	    'Color',[0.8 0.8 0.8],...
	    'xtick',[],...
	    'ytick',[]);
    colormap(cmp);
    chb=colorbar('location','eastoutside');
    set(chb,'Fontsize',14,'ticklength',0.02);

    ctl=sprintf('Fld: %s, dF=Gr-NoGr, lr=%i %i, %i',pfld,S1(ilr).depth_av(1),...
		S1(ilr).depth_av(2),yr);
    title(ctl);

    btx='dltS_GreenldExp.m';
    bottom_text(btx,'pwd',1);
  end
  
end



return