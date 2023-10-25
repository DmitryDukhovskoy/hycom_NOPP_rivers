function sub_plot_S_GrShelf(nf,SZ,SZGR,js1,js2,jw1,jw2,btx,stlS,stlW,POS);
% Plot Summer and winter S fields (months: js1:js2)
% sections
nsct=length(SZGR);
ZZ=SZGR(1).ZZlevels;
ZM=SZGR(1).ZM;
nzm=length(ZM);
nzz=length(ZZ);

% For plotting add an extra layer 
%Un(nzz,:)=Un(nzz-1,:);

dw=0.4;
dh=0.38;
%POS=[0.07 0.56 dw dh; ...
%     0.53 0.56 dw dh; ...
%     0.07 0.08 dw dh; ...
%     0.53 0.08 dw dh];

Scntr=[30:0.2:34.8,34.9,34.95,35];

figure('Position',[912 529 1544 804]);
clf;
for ip=1:nsct
  xl1=[];
  xl2=[];
  yl1=[];
  yl2=[];
  Ssm=SZ(ip).ssm(js1:js2,:,:);
  Ssm=squeeze(nanmean(Ssm,1));
  
  Ssm = sub_fill_bottom_nans(Ssm);

  LL=SZGR(ip).Dist_origin*1e-3; % km
  nm=SZGR(ip).Name;
  Hb=SZGR(ip).Hbottom;
%keyboard  
  if ip==2
    xl1=0;
    xl2=90;
%    xl2=150;
    yl1=-1000;
%    yl1=-2500;
    yl2=0;
  end
  
  if isempty(xl1)
    xl1=LL(1);
    xl2=LL(end);
  end
  if isempty(yl1);
    yl1=min(Hb)-10;
    yl2=0;
  end
  
  nmfld='S';
  hcb=[];
  if ip==2, hcb=[0.94 0.57 0.006 0.38]; end;
  if ip==4, hcb=[0.94 0.08 0.006 0.38]; end;
  
  sub_plot_grsctS(ip,LL,ZM,Ssm,xl1,xl2,yl1,yl2,Hb,hcb,POS);
  contour(LL,ZM,Ssm,Scntr,'k');
  contour(LL,ZM,Ssm,[34 34],'k','Linewidth',1.6);

  stl=sprintf('%s, %s',stlS,nm);
  title(stl,'Interpreter','none');
  
end

%set(gcf,'Position',[912 529 1544 804]);
bottom_text(btx,'pwd',1);


% Winter
%jj1=11;
%jj2=3;
figure('Position',[912 529 1544 804]);
clf;
for ip=1:nsct
  xl1=[];
  xl2=[];
  yl1=[];
  yl2=[];
  Ssm=SZ(ip).ssm([1:jw2,jw1:12],:,:);
  Ssm=squeeze(nanmean(Ssm,1));

  Ssm = sub_fill_bottom_nans(Ssm);
  
  LL=SZGR(ip).Dist_origin*1e-3; % km
  nm=SZGR(ip).Name;
  Hb=SZGR(ip).Hbottom;
  
  if ip==2
    xl1=0;
    xl2=90;
%    xl2=150;
    yl1=-1000;
%    yl1=-2500;
    yl2=0;
  end
  
  if isempty(xl1)
    xl1=LL(1);
    xl2=LL(end);
  end
  if isempty(yl1);
    yl1=min(Hb)-10;
    yl2=0;
  end
  
%  stl=sprintf('008-%3.3i, %s, S, %i-%i,mo:%i-%i',expt,nm,YR1,YR2,jj1,jj2);
  nmfld='S';
  hcb=[];
  if ip==2, hcb=[0.94 0.57 0.006 0.38]; end;
  if ip==4, hcb=[0.94 0.08 0.006 0.38]; end;
  
  sub_plot_grsctS(ip,LL,ZM,Ssm,xl1,xl2,yl1,yl2,Hb,hcb,POS);
  contour(LL,ZM,Ssm,Scntr,'k');
  contour(LL,ZM,Ssm,[34 34],'k','Linewidth',1.6);

  stl=sprintf('%s, %s',stlW,nm);
  title(stl,'Interpreter','none');
  
end

%set(gcf,'Position',[912 529 1544 804]);
bottom_text(btx,'pwd',1);


return
