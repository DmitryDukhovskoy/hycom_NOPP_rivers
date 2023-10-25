function GC = sub_greenl_coast08(HG);
% IJ indices of Greenland coast
% for clacluating upwelling index
HH = HG;
[mm,nn]=size(HG);
[II,JJ] = meshgrid([1:nn],[1:mm]);


% Defin polygon with Greenland
PGR = [  749        1067
         814        1063
         885        1062
         924        1004
         956         976
         925         918
         918         840
         909         781
         899         732
         883         683
         877         645
         836         606
         702         538
         676         526
         617         417
         584         417
	 554         454
         542         489
         536         607
         553         690
         584         776
         617         869
         623         925
         590         954
         585         978
         597        1014
         619        1038
	 629        1044
         662        1037
         680        1036
         688        1057
         738        1065];
%
IN = inpolygon(II,JJ,PGR(:,1),PGR(:,2));
HG(~IN)=nan;
figure(10); clf;
[ca,ha] = contour(HG,[0 0],'k');
close(10);

% select Greenland contour:
GC = struct;
lca = length(ca);
i1=0;
i2=0;
cc = 0;
while i1<lca
  i1=i1+1;
  npp = ca(2,i1);
  i2=i1+npp;
  i1=i1+1;
  if npp>3000
    cc = cc+1;
    GC.X=ca(1,i1:i2);
    GC.Y=ca(2,i1:i2);
  end
end

if cc>1, error('More than 1 Greenland contour found ...'); end
% Find local normalunit vectors
% and closest ocean point
Iocn = find(HH<0 & ~isnan(HG));
[Jo,Io]=ind2sub(size(HH),Iocn);
X=GC.X;
Y=GC.Y;
ncr = length(X);
fprintf('Finding ocean pts and offshore normal ...\n');
for ik=1:ncr
  x1=X(ik);
  y1=Y(ik);
  dd=sqrt((x1-Io).^2+(y1-Jo).^2);
  idx=find(dd==min(dd));
  
  GC.Iocn(ik)=Io(idx);
  GC.Jocn(ik)=Jo(idx);

% Find normal to local tangent line
% approximated as a straight
% line between 2 midpoints
% x0-x1 and x1-x2
  if ik==1
    ikb = ncr+(ik-1);
    x0=GC.X(ikb);
    y0=GC.Y(ikb);
    x2=GC.X(ik+1);
    y2=GC.Y(ik+1);
  elseif ik>=ncr
    ikb = ik-ncr+1;
    x0=GC.X(ik-1);
    y0=GC.Y(ik-1);
    x2=GC.X(ikb);
    y2=GC.Y(ikb);
  else
    x0=GC.X(ik-1);
    y0=GC.Y(ik-1);
    x2=GC.X(ik+1);
    y2=GC.Y(ik+1);
  end

  
  a0 = y0;
  d10  = (y1-y0)/(x1-x0);   % [y1,y0]
  d21  = (y2-y1)/(x2-x1);   % [y2,y1]

  xp01 = 0.5*(x1+x0);
  yp01 = 0.5*(y1+y0);
  xp12 = 0.5*(x1+x2);
  yp12 = 0.5*(y1+y2);
  aa   = (yp12-yp01)/(xp12-xp01);  
  c  = y1-aa*x1;
  
% Normal line:
% y = anm*x+cnm
  anm = -1/aa;
  cnm = y1-anm*x1;
% Find unit normal directed off land
% Make 1 step in positive direction along the normal line
% and in negatve dir - see where ocean point is
% then find unit normal vector directed offshore
% normal in positive direction
  alf = atan2(anm,1); % angle of the normal line
  xps = x1+cos(alf); % x coord of unit vector in pos dir
  yps = y1+sin(alf);
%  hps = HH(round(yps), round(xps)); % 
  inps= inpolygon(xps,yps,X,Y);
  xng = x1+cos(alf+pi); % x coord of unit vector in neg dir
  yng = y1+sin(alf+pi);
  inng= inpolygon(xng,yng,X,Y); % pnt is inside Greenland coastline
%  hng = HH(round(yng), round(xng)); % 
  
  if inps*inng 
    fprintf('Could not determine off shore direction ...\n');
    keyboard
  end
  
  if inng
    xnrm = cos(alf);
    ynrm = sin(alf);
  else
    xnrm = cos(alf+pi);
    ynrm = sin(alf+pi);
  end    
    
  GC.xnrm_offshore(ik) = xnrm;
  GC.ynrm_offshore(ik) = ynrm;
  
% draw approximated tangent line 
  f_chck=0;
  if f_chck>0
   xl1=min([x0,x1,x2]);
   xl2=max([x0,x1,x2]);
   dx=(xl2-xl1)/10;
   xx=[xl1:dx:xl2+0.2];
   yy = aa*xx+c;
   ynm= anm*xx+cnm;
   
   figure(11); clf;
   hold on;
   plot(GC.X,GC.Y,'r.-')
   plot(x0,y0,'k*');
   plot(x1,y1,'b*');
   plot(x2,y2,'g*');
   axis('equal');
   plot(xx,yy,'.-');
   plot(xx,ynm,'c.-');
   plot(xps,yps,'g+');
   plot(xng,yng,'ro');
   keyboard
  end
  
end




return