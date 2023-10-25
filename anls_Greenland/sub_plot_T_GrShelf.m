function sub_plot_T_GrShelf(nf,SZ,TZGR,js1,js2,jw1,jw2,btx,stlS,stlW,POS);
% Plot Summer and winter T fields (months: js1:js2)
% sections
nsct=length(TZGR);
ZZ=TZGR(1).ZZlevels;
ZM=TZGR(1).ZM;
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

Tcntr=[-2:0.5:8];

c1=-2;
c2=8;

figure('Position',[912 529 1544 804]); 
clf;
for ip=1:nsct
  xl1=[];
  xl2=[];
  yl1=[];
  yl2=[];
  Tsm=SZ(ip).tsm(js1:js2,:,:);
  Tsm=squeeze(nanmean(Tsm,1));
  
  Tsm = sub_fill_bottom_nans(Tsm);

  LL=TZGR(ip).Dist_origin*1e-3; % km
  nm=TZGR(ip).Name;
  Hb=TZGR(ip).Hbottom;
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
  
  nmfld='T';
  hcb=[];
  if ip==2, hcb=[0.94 0.57 0.006 0.38]; end;
  if ip==4, hcb=[0.94 0.08 0.006 0.38]; end;

  sub_plot_grsctT(c1,c2,ip,LL,ZM,Tsm,xl1,xl2,yl1,yl2,Hb,hcb,POS);
  contour(LL,ZM,Tsm,Tcntr,'k');
  contour(LL,ZM,Tsm,[0 0],'k','Linewidth',1.6);

  stl=sprintf('%s, %s',stlS,nm);
  title(stl,'Interpreter','none');
  
end

%set(gcf,'Position',[912 529 1544 804]);
bottom_text(btx,'pwd',1);


% Winter
%jj1=11;
%jj2=3;
%figure(nf+1); clf;
figure('Position',[912 529 1544 804]); 
clf;
for ip=1:nsct
  xl1=[];
  xl2=[];
  yl1=[];
  yl2=[];
  Tsm=SZ(ip).tsm([1:jw2,jw1:12],:,:);
  Tsm=squeeze(nanmean(Tsm,1));

  Tsm = sub_fill_bottom_nans(Tsm);
  
  LL=TZGR(ip).Dist_origin*1e-3; % km
  nm=TZGR(ip).Name;
  Hb=TZGR(ip).Hbottom;
  
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
  
  sub_plot_grsctT(c1,c2,ip,LL,ZM,Tsm,xl1,xl2,yl1,yl2,Hb,hcb,POS);
  contour(LL,ZM,Tsm,Tcntr,'k');
  contour(LL,ZM,Tsm,[0 0],'k','Linewidth',1.6);

  stl=sprintf('%s, %s',stlW,nm);
  title(stl,'Interpreter','none');
  
end

%set(gcf,'Position',[912 529 1544 804]);
bottom_text(btx,'pwd',1);



return
