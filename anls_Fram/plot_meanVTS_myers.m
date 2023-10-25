addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=2002;
YR2=2016;

pthdat='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_Myers/FRAM/';

btx = 'plot_meanVTS_myers.m';


fout=sprintf('%smonthly_VTS.mat',pthdat);
fprintf('Loading %s\n',fout);
load(fout);

nyr=length(VTS);
aa=squeeze(VTS(1).V(1,:,:)*0);
ss=aa;
tt=aa;
for iy=1:nyr
  for im=1:12
    aa=aa+squeeze(VTS(iy).V(im,:,:));
    ss=ss+squeeze(VTS(iy).S(im,:,:));
    tt=tt+squeeze(VTS(iy).T(im,:,:));
  end
end
VM=aa/(nyr*12);
SM=ss/(nyr*12);
TM=tt/(nyr*12);

% StDv
aa=VM*0;
ss=aa;
tt=aa;
for iy=1:nyr
  for im=1:12
    dmm=squeeze(VTS(iy).V(im,:,:));
    aa=aa+(dmm-VM).^2;
    dmm=squeeze(VTS(iy).S(im,:,:));
    ss=ss+(dmm-SM).^2;
    dmm=squeeze(VTS(iy).T(im,:,:));
    tt=tt+(dmm-TM).^2;
%    fprintf('Max t=%8.4f\n',max(max(tt)));
%    keyboard
  end
end
sgmV=sqrt(aa/(nyr*12));
sgmS=sqrt(ss/(nyr*12));
sgmT=sqrt(tt/(nyr*12));

X=VTS(1).X;
nxx=length(X);
Z=VTS(1).Z;
nzz=length(Z);
ZM=Z;
dZ=abs(diff(Z));
dZ(nzz)=dZ(nzz-1);
ZZ(1)=0;
for iz=1:nzz
  ZZ(iz+1)=ZM(iz)-0.5*dZ(iz);
end

% Get bottom depths:
dmm=squeeze(VTS(1).V(1,:,:));
for ii=1:nxx
  ib=max(find(~isnan(dmm(:,ii))));
  if isempty(ib);
    Zb(ii)=10;
  else
    Zb(ii)=ZZ(ib+1);
  end
end


VM=sub_fill_bottom_nans(VM);
SM=sub_fill_bottom_nans(SM);
TM=sub_fill_bottom_nans(TM);


% Interpolate bottom into higher res grid
%Xi=[0:nxx-1];
%dn=3;
%Xnew=[0:0.1:Xi(end)];
%Pns = piecewise_interp(Xi,Zb,dn,Xnew);
Xnew=[X(1):0.01:X(end)];
Pns = spline_cub(X,Zb,Xnew,'natural');

Xbtm=[Xnew(1),Xnew,Xnew(end)];
Zbtm=[-5000,Pns',-5000];


% 
cl1=flipud(colormap_blue(200));
cl2=colormap_red(200);
for ik=1:5
  cl1(end-(ik-1),:)=[1 1 1];
  cl2(ik,:)=[1 1 1];
end
cmp1=[cl1;cl2];
cmp1=smooth_colormap(cmp1,5);


% AWI moorings:
SBE = Fram_moorAWI;
nf=length(SBE);


% ================
% Plot mean U and std dv
% ================
figure(1); clf;
pcolor(X,Z,VM); shading flat;
colormap(cmp1);
caxis([-0.2 0.2]);
hold on;
contour(X,Z,sgmV,[0:0.01:0.5],'Color',[0 0 0]);
contour(X,Z,sgmV,[0.01 0.01],'Color',[0 0 0],'Linewidth',1.6);
%plot(X,Hb,'Color',[0 0 0],'Linewidth',4);
fill(Xbtm,Zbtm,[0 0 0]);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  dmm=abs(x-X);
  im=find(dmm==min(dmm),1);
  zb=Zb(im);
%  plot(IJ(1),IJ(2),'r.','Markersize',16);
  plot(x,zb-100,'r^','Markersize',9,'MarkerFaceColor',[1 0 0]);
  nm=SBE(ik).Name;
  text(x,zb-200,nm,'Fontsize',12,'Color',[1 0 0]);
end


hb=colorbar;
set(hb,'Fontsize',14,...
       'Position',[0.92 0.1 0.02 0.8]);

xl1=-19.2;
xl2=15.7;
yl1=-3550;
yl2=0;
set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[-3500:500:0],...
	'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'Fontsize',14);

title('0.4 NEMO mean V, StDev (0.01), 2002-2016');
bottom_text(btx,'pwd',1);


% ================
% Plot mean S
% ================
figure(2); clf;
cs1=32;
cs2=35.;
%CMP = colormap_sclr2(200,cs1,cs2);
%cmpS1=CMP.colormap;
cmpS=colormap(parula(400));

pcolor(X,Z,SM); shading flat;
colormap(cmpS);
caxis([cs1 cs2]);
hold on;
%contour(X,Z,sgmS,[0:0.1:1],'Color',[0 0 0]);
%contour(X,Z,sgmS,[0.01:0.02:0.099],'--','Color',[0 0 0]);
%contour(X,Z,sgmS,[0.01 0.01],'--','Color',[0 0 0],'Linewidth',1.6);
%contour(X,Z,SM,[34.9:0.01:34.99],'k-','Color',[0.4 0.4 0.4]);
%contour(X,Z,SM,[34.95 34.95],'k-','Color',[0.4 0.4 0.4],'Linewidth',1.6);
contour(X,Z,SM,[34.8:0.02:36.0],'k-','Color',[0.3 0.3 0.3]);
contour(X,Z,SM,[34.94 34.94],'k-','Color',[0.3 0.3 0.3],'Linewidth',1.6);
contour(X,Z,SM,[34.5 34.5],'k-','Color',[0.6 0.6 0.6],'Linewidth',1);
contour(X,Z,SM,[32.5:0.5:34.],'k-','Color',[0.6 0.6 0.6],'Linewidth',1);

%plot(X,Hb,'Color',[0 0 0],'Linewidth',4);
fill(Xbtm,Zbtm,[0 0 0]);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  dmm=abs(x-X);
  im=find(dmm==min(dmm),1);
  zb=Zb(im);
%  plot(IJ(1),IJ(2),'r.','Markersize',16);
  plot(x,zb-100,'r^','Markersize',9,'MarkerFaceColor',[1 0 0]);
  nm=SBE(ik).Name;
  text(x,zb-200,nm,'Fontsize',12,'Color',[1 0 0]);
end



hb1=colorbar;
set(hb1,'Fontsize',14,...
       'Position',[0.92 0.1 0.02 0.8]);

set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[-3500:500:0],...
	'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'Fontsize',14);

title('NEMO mean S, StDev, 2002-2016');
bottom_text(btx,'pwd',1);




% ==========================  
%   Plot Mean T 
% ==========================  
figure(5); clf;
ct1=-2;
ct2=5;
nint=400;
CMP = create_colormap5(nint,ct1,ct2);
cmpT = CMP.colormap;

pcolor(X,Z,TM); shading flat;
colormap(cmpT);
caxis([ct1 ct2]);
hold on;
%contour(X,Z,sgmT,[0:0.5:2],'Color',[0 0 0]);
%contour(X,Z,sgmT,[0.0:0.05:0.49],'--','Color',[0 0 0]);
%contour(X,Z,sgmT,[0.05 0.05],'--','Color',[0 0 0],'Linewidth',1.6);
%contour(X,Z,TM,[-1:8],'k-','Color',[0.6 0.6 0.6]);
%contour(X,Z,TM,[0 0],'k-','Color',[0.6 0.6 0.6],'Linewidth',1.6);
contour(X,Z,TM,[-2:1:8],'k-','Color',[0.6 0.6 0.6]);
contour(X,Z,TM,[0 0],'k-','Color',[0.6 0.6 0.6],'Linewidth',1.6);

%plot(X,Hb,'Color',[0 0 0],'Linewidth',4);
fill(Xbtm,Zbtm,[0 0 0]);

for ik=1:nf
  y=SBE(ik).lat_lon(1);
  x=SBE(ik).lat_lon(2);
  dmm=abs(x-X);
  im=find(dmm==min(dmm),1);
  zb=Zb(im);
%  plot(IJ(1),IJ(2),'r.','Markersize',16);
  plot(x,zb-100,'r^','Markersize',9,'MarkerFaceColor',[1 0 0]);
  nm=SBE(ik).Name;
  text(x,zb-200,nm,'Fontsize',12,'Color',[1 0 0]);
end


hb1=colorbar;
set(hb1,'Fontsize',14,...
       'Position',[0.92 0.1 0.02 0.8]);

set(gca,'tickdir','out',...
	'xtick',[-20:4:20],...
	'ytick',[-3500:500:0],...
	'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'Fontsize',14);

title('NEMO mean T, StDev, 2002-2016');
bottom_text(btx,'pwd',1);

