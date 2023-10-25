function sub_plot_UV_GrShelf(nf,UV,UZGR,js1,js2,jw1,jw2,btx,stlS,stlW,POS);
% =================
% Plotting UV
% =================
nsct=length('UZGR');
ZZ=UZGR(1).ZZlevels;
ZM=UZGR(1).ZM;
nzm=length(ZM);
nzz=length(ZZ);

% For plotting add an extra layer 
%Un(nzz,:)=Un(nzz-1,:);

%dw=0.4;
%dh=0.38;
%POS=[0.07 0.56 dw dh; ...
%     0.53 0.56 dw dh; ...
%     0.07 0.08 dw dh; ...
%     0.53 0.08 dw dh];


c1=-0.4;
c2=0.4;
    
figure(nf); clf;
for ip=1:nsct
  xl1=[];
  xl2=[];
  yl1=[];
  yl2=[];
  switch(ip)
   case(1)
    Usm=UV(ip).vsm(js1:js2,:,:);
   case(2)
    Usm=UV(ip).vsm(js1:js2,:,:);
    xl1=0;
    xl2=90;
%    xl2=150;
    yl1=-1000;
%    yl1=-2500;
    yl2=0;
   case(3)
    Usm=UV(ip).vsm(js1:js2,:,:);
   case(4)
    Usm=UV(ip).usm(js1:js2,:,:);
  end
  
  Usm=squeeze(nanmean(Usm,1));
  Usm = sub_fill_bottom_nans(Usm);
  
  LL=UZGR(ip).Dist_origin*1e-3; % km
  nm=UZGR(ip).Name;
  Hb=UZGR(ip).Hbottom;
  
  if isempty(xl1)
    xl1=LL(1);
    xl2=LL(end);
  end
  if isempty(yl1);
    yl1=min(Hb)-10;
    yl2=0;
  end
  
  hcb=[];
  if ip==4, hcb=[0.95 0.08 0.005 0.38]; end;
  sub_plot_grsctU(c1,c2,ip,LL,ZM,Usm,xl1,xl2,yl1,yl2,Hb,hcb,POS);
  contour(LL,ZM,Usm,[-0.8:0.05:-0.01],'k');
  contour(LL,ZM,Usm,[-0.2 -0.2],'k','Linewidth',1.6);
  contour(LL,ZM,Usm,[0.1:0.05:0.8],'k--');
  contour(LL,ZM,Usm,[0.2 0.2],'k--','Linewidth',1.6);
  
  stl=sprintf('%s, %s',stlS,nm);
  title(stl,'Interpreter','none');
  
end

%set(gcf,'Position',[701 362 1840 973]);
set(gcf,'Position',[912 529 1544 804]);
bottom_text(btx,'pwd',1);
  
% Winter along-flow U/V 
figure(nf+1); clf;
for ip=1:nsct
  if ip==1 | ip==2 | ip==3
    Usm=UV(ip).vsm([1:jw2,jw1:12],:,:);
  else
    Usm=UV(ip).usm([1:jw2,jw1:12],:,:);
  end
  
  xl1=[];
  xl2=[];
  yl1=[];
  yl2=[];
  switch(ip)
   case(1)
    Usm=UV(ip).vsm([1:jw2,jw1:12],:,:);
   case(2)
    Usm=UV(ip).vsm([1:jw2,jw1:12],:,:);
    xl1=0;
    xl2=90;
%    xl2=150;
    yl1=-1000;
%    yl1=-2500;
    yl2=0;
   case(3)
    Usm=UV(ip).vsm([1:jw2,jw1:12],:,:);
   case(4)
    Usm=UV(ip).usm([1:jw2,jw1:12],:,:);
  end
  
  Usm=squeeze(nanmean(Usm,1));
  Usm = sub_fill_bottom_nans(Usm);
  
  LL=UZGR(ip).Dist_origin*1e-3; % km
  nm=UZGR(ip).Name;
  Hb=UZGR(ip).Hbottom;
  
  if isempty(xl1)
    xl1=LL(1);
    xl2=LL(end);
  end
  if isempty(yl1);
    yl1=min(Hb)-10;
    yl2=0;
  end
  
%  stl=sprintf('arc08-%3.3i, %s, U m/s, ONDJFM %i-%i',expt,nm,YR1,YR2);

  nmfld='UV';
  hcb=[];
  if ip==4, hcb=[0.95 0.08 0.005 0.38]; end;
  sub_plot_grsctU(c1,c2,ip,LL,ZM,Usm,xl1,xl2,yl1,yl2,Hb,hcb,POS);
  contour(LL,ZM,Usm,[-0.8:0.1:-0.01],'k');
  contour(LL,ZM,Usm,[-0.2 -0.2],'k','Linewidth',1.6);
  contour(LL,ZM,Usm,[0.1:0.1:0.8],'k--');
  contour(LL,ZM,Usm,[0.2 0.2],'k--','Linewidth',1.6);
  stl=sprintf('%s, %s',stlW,nm);
  title(stl,'Interpreter','none');
  
end

%set(gcf,'Position',[701 362 1840 973]);
set(gcf,'Position',[912 529 1544 804]);
bottom_text(btx,'pwd',1);



return