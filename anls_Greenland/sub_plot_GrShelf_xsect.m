function sub_plot_GrShelf_xsect(dnmb,nfg,fmat,xl1,xl2,yl1,yl2,nm,Sc0);
% Plot daily vertical U,T,S
% along Gr Shelf transect
% for analysis of coastal currents

fprintf('Loading %s\n',fmat);
load(fmat);

TM=UTSZ.Time;
ZZ=UTSZ.ZZlevels;
ZM=UTSZ.ZM;
nzm=length(ZM);
nzz=length(ZZ);
LL=UTSZ.Dist_origin*1e-3; % km
nm=UTSZ.Name;
Hb=UTSZ.Hbottom;

% Find closest day to specified date
dfT=max(diff(TM));
dTM=abs(TM-dnmb);
indx=find(dTM==min(dTM),1);
dTd=abs(dnmb-TM(indx));
if dTd>dfT,
  error('Specifed date %s is not in %s\n',datestr(dnmb),fmat);
else
  fprintf('Requested date %s, plotting %s\n',...
	  datestr(dnmb),datestr(TM(indx)));
end


S=squeeze(UTSZ.Saln(indx,:,:));
T=squeeze(UTSZ.Temp(indx,:,:));
U=squeeze(UTSZ.Unrm(indx,:,:));
[mm,nn]=size(S);


% Plot S
Ssm = sub_fill_bottom_nans(S);
% Filter to get rid of z-level disc.
WT=[0.3, 0.3, 0.25;...
    0.3, 0.3, 0.3;...
    0.3, 0.3, 0.3];
WT=WT./(sum(sum(WT)));
[nw,mw]=size(WT);
di=floor(nw/2);
dj=di;

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


Scntr=[30:0.1:34.8,34.9,34.95,35];
%Sc0=34.0;

figure(nfg); clf;
POS=[0.08 0.25 0.82 0.6];
hcb=[0.92 0.35 0.008 0.4];
ip=1;
sub_plot_grsctS(ip,LL,ZM,Ssm,xl1,xl2,yl1,yl2,Hb,hcb,POS);
contour(LL,ZM,Ssm,Scntr,'k');
contour(LL,ZM,Ssm,[Sc0 Sc0],'k','Linewidth',1.6);

dv=datevec(dnmb);
stl=sprintf('%s, %4.4i/%2.2i/%2.2i, S0=%4.2f',nm,dv(1:3),Sc0);
title(stl,'Fontsize',12,'Interpreter','none');

axes('Position',[0.08 0.1 0.4 0.05]);
text(0,0,fmat,'Interpreter','none');
set(gca,'visible','off');


return