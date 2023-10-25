function dvuf = sub_fltr(dvu,pgrd,Hmsk);
% Spatial filtering 
% pgrd - # of grid points
% filter is pgrd x pgrd
% exponential
[mm,nn] = size(Hmsk);
IN      = find(Hmsk>0);
nin     = length(IN);
dvuf    = dvu*nan;
di      = (pgrd-1)/2;

[X,Y] = meshgrid([-di:di],[-di:di]);

sx = di/3;
sy = di/3;
rr = 0.009;
Fw = 1/(2*pi*sx*sy*sqrt(1-rr^2))*...
     exp(-1/(2*(1-rr^2))*(X.^2/sx^2+Y.^2/sy^2-...
			  2*rr*X.*Y/(sx*sy)));
% Make sure that total weight = 1
nsm = sum(sum(Fw));
Fw  = Fw/nsm;
fprintf('Filtering, %i pnts\n',pgrd);

for ip=1:nin
  if mod(ip,2000)==0
    fprintf(' ::: Filtering, %4.1f%% ...\n',ip/nin*100);
  end
  
  I0=IN(ip);
  [j0,i0]=ind2sub(size(Hmsk),I0);
  if isnan(dvu(j0,i0)); continue; end;

%  j1 = max([j0-di,1]);
%  j2 = min([j0+di,mm]);
%  i1 = max([i0-di,1]);
%  i2 = min([i0+di,nn]);
  j1 = j0-di;
  j2 = j0+di;
  i1 = i0-di;
  i2 = i0+di;
  if j1>0 & j2<=mm &...
     i1>0 & i2<nn
    A = dvu(j1:j2,i1:i2);
  else
    for jj=j1:j2
      jc = jj-j1+1;
      for ii=i1:i2
	ic = ii-i1+1;
	if jj<1 | ii<1 |...
	   jj>mm | ii>nn
	  A(jc,ic) = nan;
	else
	  A(jc,ic) = dvu(jj,ii);
	end
      end
    end
  end
  
  inan = find(isnan(A));
  A(inan)=nanmean(nanmean(A));
  
  af = sum(sum(A.*Fw));
  dvuf(j0,i0)=af;
end






return