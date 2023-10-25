    function Psi = stream_fn(Tx,Ty,DX,DY,i0,j0,mthd,LM);
% Calculate stream function for
% mass fluxes Tx, Ty that are
% u*layer thickness 
% Tx, Ty can be depth integrated over
% several layers
% LM - land mask: 0 - land, 1 - ocean
%      make Psi = 0 over the land
% DX, DY - grid spacing, not coordinates!
% mthd = 'simps' - use simpson integration
%                  otherwise use trapezoid quadrature
%
% Define a streamfunction, Psi(x,y), such that 
%	u = -d(Psi)/dy and v = d(Psi)/dx. 		(1)
% By definition:	 du/dx+dv/dy = 0			(2)
% so the field is non-divergent. 
% Given (1), u=-d(Psi)/dy -> Psi = -Int{y0:y}u dy + a(x)	(3)
%and therefore:	d(Psi)/dx = -Int{y0:y}du/dx dy + da/dx
%Given (2), d(Psi)/dx = Int{y0:y}(dv/dy) dy + da/dx
%		     = v(x,y) - v(x,y_0) + da/dx	(4)
%But also, by definition:   d(Psi)/dx = v(x,y)
%Therefore: 	da/dx = v(x,y_0)
%		a(x) = Int{x_0:x}v(x,y_0) dx		(5)
%From (3) and (5):  
%      Psi(x,y) = -Int{y0:y}u(x,y) dy + Int{x_0:x}v(x,y_0) dx  (6)
%
%Algorithm: Choose an (x_0,y_0). Integrate the second (indefinite)
%integral of (6) along y_0. The result is a function of x. Then at
%each x, indefinitely integrate the first integral of (6) along y.
%The streamfunction is the sum of those two.
% 
% FOr better accuracy: I estimate psi by intgrating u over dy 
% and then psi by integrating v over dx
% and psi is the average of two
%
% Integration - quadrature: trapeziod gives inaccurate integration
% use simpson cumulative summation for regularly spaced data
%
% D. Dukhovskoy, 12/12/2014
%
% To check stream function - use check_streamfn.m
%
fprintf('Calculating Psi stream function ...\n');

[mm,nn]=size(Tx);
clear Psi*
%keyboard

Tx(isnan(Tx)) = 0;
Ty(isnan(Ty)) = 0;

if ~isempty(mthd) & strncmp(mthd,'sim',3)
  fprintf('stream_fn:  Simpson integration \n');
  
% Calculate Psi1:
  dx=DX(j0,:);
  X=cumsum(dx);
% Interpolate into regular grid:
  ty=Ty(j0,:);
  dxx=(X(end)-X(1))/(2*(nn-1));
  Xi=[X(1):dxx:X(end)];
  tyi=interp1(X,ty,Xi,'spline');
  D1i = cumsimp(tyi)*dxx; % integr(x0:x)V(x,y0)*dx =  a(x)
  D1 = interp1(Xi,D1i,X);
% Integration of U over Y at every x pnt:
  for i=1:nn
    dy=DY(:,i);
    Y=cumsum(dy);
    tx=Tx(:,i);
    dyy=(Y(end)-Y(1))/(2*(mm-1));
    Yi=[Y(1):dyy:Y(end)];
    txi=interp1(Y,tx,Yi,'spline');
    D2 = -cumsimp(txi)*dyy;
    dmm=D1(i)+D2;
    dmm2=interp1(Yi,dmm,Y,'spline');
%    landM=LM(:,i);
%    dmm2 = check_land(dmm2,landM);
    Psi1(:,i)=dmm2;
  end

% Calculate Psi2:
  dy=DY(:,i0);
  Y=cumsum(dy);
  tx=Tx(:,i0);
  dyy=(Y(end)-Y(1))/(2*(mm-1));
  Yi=[Y(1):dyy:Y(end)];
  txi=interp1(Y,tx,Yi,'spline');
  D1i = -cumsimp(txi)*dyy;
  D1 = interp1(Yi,D1i,Y);
  
  for j=1:mm
    dx-DX(j,:);
    X=cumsum(dx);
    ty=Ty(j,:);
    dxx=(X(end)-X(1))/(2*(nn-1));
    Xi=[X(1):dxx:X(end)];
    tyi=interp1(X,ty,Xi,'spline');
    D2 = cumsimp(tyi)*dxx;
    dmm=D1(j)+D2;
    dmm2=interp1(Xi,dmm,X,'spline');
%    landM=LM(j,:);
%    dmm2 = check_land(dmm2,landM);
     Psi2(j,:)=dmm2;
  end
%  keyboard
  Psi=0.5*(Psi1+Psi2);
  
else
  fprintf('Trapezoid integration ...\n');
%
% Calculate Psi1
  for i=1:nn
    dx=DX(j0,1:i);
    X=cumsum(dx);
    ty=Ty(j0,1:i);
    if i==1,
      D1=nansum(ty.*dx);
    else
      D1 = trapz(X,Ty(j0,1:i));
    end
    for j=1:mm
      dy=DY(1:j,i);
      Y=cumsum(dy);
      tx=Tx(1:j,i);
      if j==1
	D2=-nansum(Tx(1:j,i).*dy);
      else
	D2 = -trapz(Y,Tx(1:j,i));
      end

      Psi1(j,i)=D1+D2;
    end
  end

  % Calculate Psi2
  fprintf('Calculating Psi2 ...\n');
  for j=1:mm
    dy=DY(1:j,i0);
    Y=cumsum(dy);
    tx=Tx(1:j,i0);
    if j==1
      D1=-nansum(Tx(1:j,i0).*dy);
    else
      D1 = -trapz(Y,Tx(1:j,i0));
    end

    for i=1:nn
      dx=DX(j,1:i);
      X=cumsum(dx);
      ty=Ty(j,1:i);
      if i==1
	D2=nansum(Ty(j,1:i).*dx);
      else
	D2=trapz(X,Ty(j,1:i));
      end
    end
    Psi2(j,i)=D1+D2;
  end

  %keyboard
  Psi=0.5*(Psi1+Psi2);
end

return

function psiL = check_land(psi,lm);
% Stream function over the land should be 
% brought to 0 by subtracting/adding dPhsi
% same dPsi should be subtracted/added 
% on the opposite end of the landmass
psiL=psi;
J=find(lm==1);
if isempty(J); return; end; % no ocean points
I=find(lm==0);
if isempty(I); return; end; % no land points
di=diff(I);
i1=0; % on the 1st land point
i2=0; % on the last land point

while i2<length(lm);
  la=lm;
  if (i2>0);
    la(1:i2)=1;
  end
  i1=min(find(la==0)); % land starts
  if isempty(i1); break; end;  % no land
  lb=lm;
  lb(1:i1)=0;
  i2=min(find(lb==1))-1; % land ends
  if isempty(i2); i2=1e6; end;

  if i1==1, continue; end;
  dpsi=-psiL(i1);
  psiL(i1:end)=psiL(i1:end)+dpsi;
  
end

return

