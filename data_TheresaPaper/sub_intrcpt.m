function [IxGr,IxGt,IxGt0]  = sub_intrcpt(IIs,JJs,IS,JS,Igs,Jgs,IG,JG,icst,jcst,xP,yP);
% Find intercept of 2 contours and 
% connect at the corner/intercept
% IIs,JJs - n x 2 - start-end points of segments going through U/V points along the gate section
% IS,JS - corresponding grid indices (center of gird cells or  p-points)
% 
% Igs,Jgs - n x 2 U/Vpoints segment start/end points in index space, Gr contour
% IG, JG - corresponding grid incides at p-points, Gr Contour
%

% Intersection pnt:
% Find the first intercept pont going from the Gr coast along the u/v points
% in case there are several (when Gr contour and gate segments overlap over several points
%  at the intersection, gowing parallel or perpendicular to the coast)
%
% Section end closest to the coast:
ddC=min(sqrt((IIs-icst).^2+(JJs-jcst).^2)');

ni=length(IIs);
im=[];
if ddC(1)>ddC(end)
		for iss=ni:-1:1
				dd=sqrt((Igs(:,1)-IIs(iss,1)).^2+(Jgs(:,1)-JJs(iss,1)).^2);
				mdd=min(dd);
				if mdd==0, break; end
		end
else
		for iss=1:ni
				dd=sqrt((Igs(:,1)-IIs(iss,2)).^2+(Jgs(:,1)-JJs(iss,2)).^2);
				mdd=min(dd);
				if mdd==0, break; end;
		end
end
if abs(mdd)>1e-20
		fprintf('Could not locate GrContour - Gate intercept, section %i: %s\n');
		keyboard;
end
%
% igr,jgr - intercept of 2 contours, grid index
IGr=find(dd==mdd);  % Gr contour intercept
IGt=iss;            % gate intercept
igr=IG(iss);
jgr=JG(iss);
clear iss



if ddC(1)<ddC(end)
		st1=1; % Greenl coast at index 1 of the segment
  igt2=IGt+1;
  IxGt0=1;
  iSpnt=IIs(1,1);
  jSpnt=JJs(1,1);
else
		st1=0;
  igt2=IGt-1;
  IxGt0=length(IIs);
  iSpnt=IIs(end,1);
  jSpnt=JJs(end,1);
end
% Segments Gate section, 1 before the intercept, 1 after the intercept
% both segments should have similar pnt at the intercept pnt
xsS1=IIs(IGt,1);
xsE1=IIs(IGt,2);
ysS1=JJs(IGt,1);
ysE1=JJs(IGt,2);

xsS2=IIs(igt2,1);
xsE2=IIs(igt2,2);
ysS2=JJs(igt2,1);
ysE2=JJs(igt2,2);

%keyboard

%
% Gate section: Start from the coast to the interception point
% find intercepting Gr contour segments
ns=length(IIs);
ng=length(Igs);

% One of the endpoints of the intercept segment
% should intercept with the Gr contour 
isgm=1;
d1=sqrt((Igs(:,1)-IIs(IGt,isgm)).^2+(Jgs(:,1)-JJs(IGt,isgm)).^2);
d2=sqrt((Igs(:,2)-IIs(IGt,isgm)).^2+(Jgs(:,2)-JJs(IGt,isgm)).^2);
dmn1=min(d1);
dmn2=min(d2);

if dmn1>0 | dmn2>0
  isgm=2;
  d1=sqrt((Igs(:,1)-IIs(IGt,isgm)).^2+(Jgs(:,1)-JJs(IGt,isgm)).^2);
  d2=sqrt((Igs(:,2)-IIs(IGt,isgm)).^2+(Jgs(:,2)-JJs(IGt,isgm)).^2);
  dmn1=min(d1);
  dmn2=min(d2);
end

%
% There should be 2 segments on Gr UV-point contour connected
% at the interception point, the task is to pick 1 of these 
% to close off the box
if dmn1>0 | dmn2>0
  fprintf('sub_intercpt: \n');
  fprintf('2 segments interceping Gate section not found\n');
  keyboard
end

i1=find(d1==dmn1);
i2=find(d2==dmn2);

% First segment:
xS1=Igs(i1,1);
xE1=Igs(i1,2);
yS1=Jgs(i1,1);
yE1=Jgs(i1,2);

xS2=Igs(i2,1);
xE2=Igs(i2,2);
yS2=Jgs(i2,1);
yE2=Jgs(i2,2);

%  This should not be 
% Select segment that is inside the box wrt to Gate section
% and coast point
D  = sign(orientation(xsS1,ysS1,xsE1,ysE1,xP,yP));
R1 = sign(orientation(xsS1,ysS1,xsE1,ysE1,xS1,yS1));
if R1==0  % intrcept point
  R1 = sign(orientation(xsS1,ysS1,xsE1,ysE1,xE1,yE1));
end

R2 = sign(orientation(xsS1,ysS1,xsE1,ysE1,xS2,yS2));
if R2==0  % intrcept point
  R2 = sign(orientation(xsS1,ysS1,xsE1,ysE1,xE2,yE2));
end
%keyboard
%
% Greenl Segment can coincide with the gate section line
% in this case R1=0 or R2=0 and none =D
if (R1~=D & R2~=D) & R1==0
  R1=D;
end
if (R1~=D & R2~=D) & R2==0
  R2=D;
end;

if R1~=D & R2~=D 
  fprintf('sub_intrcp: could not locate correct segment inside box R1=%i R2=%i D=%i\n',R1,R2,D);
  keyboard
end
if R1==D & R2==D
  fprintf('sub_intrcp: both segments inside box: R1=%i R2=%i D=%i\n',R1,R2,D);
  keyboard
end

if R1==D
  IxGr=i1;
else
  IxGr=i2;
end
xS=Igs(IxGr,1);
xE=Igs(IxGr,2);
yS=Jgs(IxGr,1);
yE=Jgs(IxGr,2);

%
% Now check that the Gate section segment is the right one:
% wrt to the gr segment
d1s=sqrt((xsS1-iSpnt).^2+(ysS1-jSpnt).^2);
d1e=sqrt((xsE1-iSpnt).^2+(ysE1-jSpnt).^2);
d2s=sqrt((xsS2-iSpnt).^2+(ysS2-jSpnt).^2);
d2e=sqrt((xsE2-iSpnt).^2+(ysE2-jSpnt).^2);

if min([d1s,d1e])<min([d2s,d2e])
  IxGt=IGt;
else
  IxGt=igt2;
end;
xsS=IIs(IxGt,1);
xsE=IIs(IxGt,2);
ysS=JJs(IxGt,1);
ysE=JJs(IxGt,2);
iss=IS(IxGt);
jss=JS(IxGt);



f_chck=0;
if f_chck==1
  figure(12); clf;
  hold on;
  plot(Igs(:,1),Jgs(:,1),'.-')
  plot(IIs(:,1),JJs(:,1),'.-');
  axis('equal')
%  plot(igr,jgr,'g.','Markersize',18); % intercept p-point
  plot(iss,jss,'.','Markersize',18); % p-point Gate section
  text(iss,jss+0.1,sprintf('GateSct p-pnt (%i %i)',iss,jss),'Fontsize',12);
  plot(icst,jcst,'b.','Markersize',14);
  plot(xP,yP,'r*');

  plot([xS1 xE1],[yS1 yE1],'-','Color',[0.5 0.5 0.5],'Linewidth',4); % possible segment Gr cntr
  plot([xS2 xE2],[yS2 yE2],'-','Color',[0.4 0.4 0.4],'Linewidth',4); % possible segment Gr cntr
  plot([xsS1 xsE1],[ysS1 ysE1],'-','Color',[0.8 0.3 0],'Linewidth',4); % selected segm gate
  plot([xsS2 xsE2],[ysS2 ysE2],'-','Color',[0.8 0.3 0],'Linewidth',4); % selected segm gate
  plot([xS xE],[yS yE],'-','Color',[0 0.8 1],'Linewidth',1.8);  % selected segment Gr cntr
  plot([xsS xsE],[ysS ysE],'-','Color',[0 0.8 1],'Linewidth',1.8); % selected segm gate
  ig=IG(IxGr); 
  jg=JG(IxGr);
  plot(ig,jg,'.','Markersize',16); % P-point grid cell for Gr Contour flux
  text(ig,jg-0.1,sprintf('GrCntr p-pnt (%i %i)',ig,jg),'Fontsize',12);

  btx='sub_intrcp.m';
  bottom_text(btx,'pwd',1);
%keyboard
end



return
