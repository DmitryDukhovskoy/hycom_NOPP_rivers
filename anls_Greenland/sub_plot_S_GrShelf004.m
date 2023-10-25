function sub_plot_S_GrShelf004(nf,UTSZ,Ssm,btx,stl,Sc0);
% Plot Summer and winter S fields (months: js1:js2)
% sections
% Sc0 - salinity defining the plume front
ZZ=UTSZ.ZZlevels;
ZM=UTSZ.ZM;
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

Scntr=[30:0.1:34.8,34.9,34.95,35];

figure('Position',[912 529 1544 804]);
clf;

xl1=[];
xl2=[];
yl1=[];
yl2=[];
  
Ssm = sub_fill_bottom_nans(Ssm);
% Filter to get rid of z-level disc.
WT=[0.3, 0.3, 0.25;...
    0.3, 0.3, 0.3;...
    0.3, 0.3, 0.3];
WT=WT./(sum(sum(WT)));
[nw,mw]=size(WT);
di=floor(nw/2);
dj=di;

[mm,nn]=size(Ssm);
for it=1:6
dmm=Ssm;
for ii=2:nn-1
  for jj=2:mm-1
    i1=ii-di;
    i2=ii+di;
    j1=jj-dj;
    j2=jj+dj;
    i1=max([1,i1]);
    i2=min([nn,i2]);
    j1=max([1,j1]);
    j2=min([mm,j2]);
    A=Ssm(j1:j2,i1:i2);
    df=sum(sum(A.*WT));
    dmm(jj,ii)=df;
  end
end
Ssm=dmm;
end



LL=UTSZ.Dist_origin*1e-3; % km
nm=UTSZ.Name;
Hb=UTSZ.Hbottom;

%keyboard  
xl1=0;
xl2=90;
yl1=-1000;
yl2=0;
  
nmfld='S';

%  hcb=[];
dw=0.6;
dh=0.38;
POS=[0.1 0.56 dw dh];
hcb=[0.74 0.57 0.006 0.38]; 

ip=1;   
sub_plot_grsctS(ip,LL,ZM,Ssm,xl1,xl2,yl1,yl2,Hb,hcb,POS);
contour(LL,ZM,Ssm,Scntr,'k');
contour(LL,ZM,Ssm,[Sc0 Sc0],'k','Linewidth',1.6);

stl=sprintf('%s, %s',stl,nm);
title(stl,'Interpreter','none');
  


return
