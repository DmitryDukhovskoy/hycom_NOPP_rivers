function GC = sub_greenl_isobath(HG,LON,LAT);
% IJ indices of Greenland coast
% for clacluating upwelling index
HH = HG;
[mm,nn]=size(HG);
[II,JJ] = meshgrid([1:nn],[1:mm]);


% Define polygon with Greenland
PGR = [  573         938
         598        1061
         616        1069
         643        1075
         689        1079
         755        1085
         759        1081
	 840        1073
	 899        1060
         949        1054
         990         991
        1000         902
        1003         862
         994         817
         963         774
         918         682
         912         627
         888         611
         848         593
         787         523
         754         518
         692         486
         652         443
         618         388
         573         396
         536         442
         510         537
         514         600
	 516         622
	 513         641
         517         658
         508         725
         523         750
         521         790
         553         821
         557         854
         557         894]; 
%
IN = inpolygon(II,JJ,PGR(:,1),PGR(:,2));
HG(~IN )=-10000;
figure(10); clf;
contour(HG,[0 0],'k');
hold on
[ca,ha] = contour(HG,[-800 -800],'b');
close(10);

%keyboard

% select around-Greenland contour:
GC = struct;
lca = length(ca);
i1=0;
i2=0;
cc = 0;
while i2<lca
  i1=i2+1;
  npp = ca(2,i1);
  i2=i1+npp;
  i1=i1+1;
  x1=ca(1,i1);
  y1=ca(2,i1);
  x2=ca(1,i2);
  y2=ca(2,i2);
  daa=sqrt((x1-x2).^2+(y1-y2).^2);
  if npp>1000 & daa<1
    cc = cc+1;
    GC(cc).X=ca(1,i1:i2);
    GC(cc).Y=ca(2,i1:i2);
  end
end

if cc>1, error('More than 1 Greenland contours found ...'); end
%
% Find HYCOM index
%Iocn = find(HH<0 & (HG>-1000));
Iocn = find(HG>-1000);
[Jo,Io]=ind2sub(size(HH),Iocn);
X=GC.X;
Y=GC.Y;
ncr = length(X);
fprintf('Finding Greenland contour along isobath ...\n');
dik = round(ncr*0.005);
cca=0;
for ik=1:dik:ncr
  x1=X(ik);
  y1=Y(ik);
  dd=sqrt((x1-Io).^2+(y1-Jo).^2);
  idx=find(dd==min(dd),1);
  
  cca=cca+1;
  GC.Iocn(cca)=Io(idx);
  GC.Jocn(cca)=Jo(idx);
end
cca=cca+1;
GC.Iocn(cca)=GC.Iocn(1);
GC.Jocn(cca)=GC.Jocn(1);

% Reorder points starting from xS,yS
% and go in c/clckwise direction
X  = GC.Iocn;
Y  = GC.Jocn;
X(end) = nan; % closed contour - repeated 1st pnt
Y(end) = nan;
nx = length(X);
xS = 691;
yS = 1080;
x0 = 724;
y0 = 781;
d  = sqrt((X-xS).^2+(Y-yS).^2);
iS = find(d==min(d),1); 
%
% Assuming that the points go 
% in order, otherwise need 
% to write a search algorithm
% for unordered group of points
clear Xn Yn
Xn=X(iS:end-1);
Yn=Y(iS:end-1);
Xn=[Xn,X(1:iS-1)];
Yn=[Yn,Y(1:iS-1)];
Xn(end+1)=Xn(1);
Yn(end+1)=Yn(1);

GC.Iocn = Xn;
GC.Jocn = Yn;

% Find grid points between the segment points
IIs = [];
JJs = [];
for ii=1:cca-1
  i1=GC.Iocn(ii);
  i2=GC.Iocn(ii+1);
  j1=GC.Jocn(ii);
  j2=GC.Jocn(ii+1);
  [I,J]=sub_xsct_indx(i1,j1,i2,j2);
  I=I(:);
  J=J(:);
  
  IIs=[IIs;I];
  JJs=[JJs;J];
end

ni = length(IIs);
% Check for repeated points
dii=diff(IIs);
djj=diff(JJs);
sdd=sqrt(dii.^2+djj.^2);
I=find(sdd>0);
IIs=IIs(I);
JJs=JJs(I);


ni = length(IIs);

GC.cntr_Iindx=IIs;
GC.cntr_Jindx=JJs;

GC = rmfield(GC,'Iocn');
GC = rmfield(GC,'Jocn');
GC = rmfield(GC,'X');
GC = rmfield(GC,'Y');

% Find unit normal pointed inside the contour
% norm is in the middle of a segment, i.e. either nrm_y=0 or nrm_x=0
clear Nrm
for ip=1:ni-1
  i0 = GC.cntr_Iindx(ip);
  j0 = GC.cntr_Jindx(ip);
  i1 = GC.cntr_Iindx(ip+1);
  j1 = GC.cntr_Jindx(ip+1);

  if i1==i0 % Y section, norm is along X axis
    ing=i1-1;
    ips=i1+1;
    jng=0.5*(j1+j0);
    jps=jng;
    nrm_ps=inpolygon(ips,jps,IIs,JJs);
    if nrm_ps
      Nrm(ip,1)=1;
      Nrm(ip,2)=0;
    else
      Nrm(ip,1)=-1;
      Nrm(ip,2)=0;
    end
% Check: plot norms:
%    plot([i1 i1+Nrm(ip,1)],[jps jps],'c-')
  else  % X section
    ing=0.5*(i1+i0);
    ips=ing;
    jng=j1-1;
    jps=j1+1;
    nrm_ps=inpolygon(ips,jps,IIs,JJs);
    if nrm_ps
      Nrm(ip,1)=0;
      Nrm(ip,2)=1;
    else
      Nrm(ip,1)=0;
      Nrm(ip,2)=-1;
    end
% Check: plot norms:
%    plot([ips ips],[j1 j1+Nrm(ip,2)],'r-')
    
  end
  
end
GC.Norm_in = Nrm;
    
  
% Compute distances along the contour
% Get bottom depth
clear dx Hb
dx(1) = 0;
for ip=1:ni-1
  i0 = GC.cntr_Iindx(ip);
  j0 = GC.cntr_Jindx(ip);
  i1 = GC.cntr_Iindx(ip+1);
  j1 = GC.cntr_Jindx(ip+1);
  x0 = LON(j0,i0);
  y0 = LAT(j0,i0);
  x1 = LON(j1,i1);
  y1 = LAT(j1,i1);
  dx(ip+1) = distance_spheric_coord(y1,x1,y0,x0);
  Hb(ip) = HH(j0,i0);
end
Hb(ni)=Hb(1);
dst = cumsum(dx);

GC.Distance_m = dst;
GC.Hbottom = Hb;

I = GC.cntr_Iindx;
J = GC.cntr_Jindx;
II = sub2ind(size(HH),J,I);
%[X,Y]=meshgrid([1:nn],[1:mm]);
GC.cntr_lon = LON(II);
GC.cntr_lat = LAT(II);

return