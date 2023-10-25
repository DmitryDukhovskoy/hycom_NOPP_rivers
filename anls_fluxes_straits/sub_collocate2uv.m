function CLC=sub_collocate2uv(U,T,S,dH,HH,i0,j0,i1,j1);
% collocate p-point variables with U or V pnt
% only 2D and 3D array 
%
% Use hycom approach to derive dp or dh at u/v points
% Otherwise very high transports may exist
% at the near bottom grid cells with dh--> 0
% 
% see dpudpv.F in hycom note that depthu=uvdep(pbot(i,j),pbot(i-1,j),
% where uvdep=min(a,b)
%
% Land ==0, depths >0 following hycom source code
HH(HH>=0)=0;
HH=abs(HH);

% Fir very thin layers, make U=0
dHeps=0.01; 

% Land:
%    UTS.Hb(is)=HH(j0,i0);
%    UTS.Unrm(:,cc)=V(:,j0,i0);
%    UTS.Tnrm(:,cc)=T(:,j0,i0);
%    UTS.Snrm(:,cc)=S(:,j0,i0);
%    UTS.dH(:,cc)=dH(:,j0,i0);
if abs(HH(j0,i0))<=1e-10
 	CLC.Un=U(:,j0,i0);
 	CLC.Tn=T(:,j0,i0);
 	CLC.Sn=S(:,j0,i0);
 	CLC.dHn=dH(:,j0,i0);
 	CLC.Hn=HH(j0,i0);
	 return
end


%
dm=size(U);
ndm=max(size(dm));
if ndm>3, error('sub_collocate2uv.m: Dimension > 3'); end
%keyboard
if ndm==3
  un=U(:,j0,i0);
  t1=T(:,j0,i0);
  t2=T(:,j1,i1);
  s1=S(:,j0,i0);
  s2=S(:,j1,i1);
  p1=cumsum(dH(:,j0,i0));
  p2=cumsum(dH(:,j1,i1));
  dh1=dH(:,j0,i0);
  dh2=dH(:,j1,i1);
else
  un=U(j0,i0);
  t1=T(j0,i0);
  t2=T(j1,i1);
  s1=S(j0,i0);
  s2=S(j1,i1);
  p1=cumsum(dH(j0,i0));
  p2=cumsum(dH(j1,i1));
  dh1=dH(j0,i0);
  dh2=dH(j1,i1);
end
p1=[0;p1];
p2=[0;p2];
h1=HH(j0,i0);
h2=HH(j1,i1);

% Implement HYCOM fluxes:
depthu=min(h1,h2);
clear dHn
for k=1:length(t1)
  dpk1=min([depthu, 0.5*(p1(k+1)+p2(k+1))]);
  dpk =min([depthu, 0.5*(p1(k)+p2(k))]);
  dp0 =dpk1-dpk;
  dHn(k,1)=dp0;
end

% Near land points:
if h2>=0
  t2=t1;
  s2=s1;
end


I=find(dHn<dHeps);
if ~isempty(I),
  un(I)=0;
end

% Interpolate into spurious high velocities
% These are typically near bottom, steep topography
% very thin layer
umx0=5.;
if max(abs(un))>=umx0
  imx=find(abs(un)==max(abs(un)));
  umx=un(imx);
  iB=min(find(dHn==0));
  unB=un(1:iB);
  I=find(abs(unB)==max(abs(unB)));
  JJ=find(abs(unB)<umx0);
  X=unB(JJ);
% To avoid errors in interpolation, grid should have distinct points
% i.e. different values
  dmm = dHn;
  dmm(dmm==0.)=1.e-6;
  ZZ=cumsum(dmm(1:iB));
  ZX=ZZ(JJ);
%keyboard
  Xi=interp1(ZX,X,ZZ,'pchip');
  un(1:iB)=Xi;
  umxF=un(imx);
  umxA=max(abs(un));
  

  fprintf('sub_collocate2uv: Unrealistically high U=%6.2f level=%i dh=%6.4f\n',...
          umx,I,dHn(I));
  fprintf('sub_collocate2uv: Fixed U=%8.6d, overall max|U|=%6.2f \n',  umxF,umxA);
%  keyboard
end

CLC.Un=un;
CLC.Tn=0.5*(t1+t2);
CLC.Sn=0.5*(s1+s2);
%dHn=0.5*(dh1+dh2); % these layer thickness not correct near the bottom
sH=abs(nansum(dHn));
%Hn=abs(0.5*(h1+h2)); % this is not consistent with HYCOM
Hn=abs(depthu);
CLC.dHn=dHn;
CLC.Hn=Hn;

if abs(sH>1e-10) & abs(Hn-sH)/sH>0.01
fprintf('** ERR sub_collocate2uv:  depths i=%i, j=%i, H=%6.2f dh=%6.2f\n',...
	    i0,j0,Hn,sH);
%    keyboard
end

%keyboard

return
