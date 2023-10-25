function CLC=sub_collocate2uv(U,dH,HH,i0,j0,i1,j1);
% collocate deoth with U or V pnt
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
if abs(HH(j0,i0))<=1e-10
 	CLC.Un=U(:,j0,i0);
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
  p1=cumsum(dH(:,j0,i0));
  p2=cumsum(dH(:,j1,i1));
  dh1=dH(:,j0,i0);
  dh2=dH(:,j1,i1);
else
  un=U(j0,i0);
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
for k=1:length(un)
  dpk1=min([depthu, 0.5*(p1(k+1)+p2(k+1))]);
  dpk =min([depthu, 0.5*(p1(k)+p2(k))]);
  dp0 =dpk1-dpk;
  dHn(k,1)=dp0;
end


I=find(dHn<dHeps);
if ~isempty(I),
  un(I)=0;
end

% Interpolate into spurious high velocities
% These are typically near bottom, steep topography
% very thin layer
umx0=3;
if max(abs(un))>=umx0
  imx=find(abs(un)==max(abs(un)));
  umx=un(imx);
  iB=min(find(dHn==0));
  unB=un(1:iB);
  I=find(abs(unB)==max(abs(unB)));
  JJ=find(abs(unB)<umx0);
  X=unB(JJ);
  ZZ=cumsum(dHn(1:iB));
  ZX=ZZ(JJ);
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
sH=abs(nansum(dHn));
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
