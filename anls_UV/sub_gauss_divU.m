function divU = sub_gauss_divU(HH,Utr,Vtr,dij);
% Compute area-mean divergence
% Using Gauss theorem, intgr(over boundary) U*n*dA
% Utr, Vtr are already Transports, i.e.
% Utr=u*dy*dz, Vtr=v*dx*dh
%
% HH is mask: =0 - skip point, =1 - valid point
%
% Sum transports:
% over Y faces:
% Note orientation of the normal unit vector
% directed outside the contour

[mm,nn]=size(HH);
divU=HH*nan;
Utr(isnan(Utr)) = 0;
Vtr(isnan(Vtr)) = 0;
for ii=dij+1:nn-dij-1
  if mod(ii,100)==0
    fprintf('----->  Calculating divU %4.1f...%%\n',ii/nn*100);
  end

  for jj=dij+1:mm-dij-1
    if (HH(jj,ii)==0), continue; end
    dm1 = -1*sum(Utr(jj-dij:jj+dij,ii-dij));    % west bndry
    dm2 = sum(Utr(jj-dij:jj+dij,ii+dij+1));     % east bndry
    dm3 = -1*sum(Vtr(jj-dij,ii-dij:ii+dij));    % south
    dm4 = sum(Vtr(jj+dij+1,ii-dij:ii+dij));     % north

    divU(jj,ii) = dm1+dm2+dm3+dm4;
  end
end

divU(HH==0)=nan;

return