function dvuf = sub_fltr(dvu,pgrd,Hmsk);
% Spatial filtering 
% pgrd - # of grid points
% filter is pgrd x pgrd
% Gaussian 2D filter

IN = find(Hmsk>0);
nin = length(IN);
dvuf = dvu*nan;
di = (pgrd-1)/2;

if floor(pgrd/2)*2 == di
  error('pgrid has to be odd number: %i',pgrd);
end

[X,Y]=meshgrid([-di:di],[-di:di]);

sx=di/3;
sy=di/3;
rr=0.009;
Fw = 1/(2*pi*sx*sy*sqrt(1-rr^2))*...
     exp(-1/(2*(1-rr^2))*(X.^2/sx^2+Y.^2/sy^2-...
			  2*rr*X.*Y/(sx*sy)));
% Make sure that total weight = 1
nsm = sum(sum(Fw));
Fw = Fw/nsm;
fprintf('Filtering, %i pnts\n',pgrd);

%keyboard

[jdm,idm] = size(dvu);
for ip=1:nin
  if mod(ip,100000)==0
    fprintf(' ::: Filtering, %4.1f%% ...\n',ip/nin*100);
  end
  
  I0=IN(ip);
  [j0,i0]=ind2sub(size(Hmsk),I0);
  if isnan(dvu(j0,i0)); continue; end;
 
%  fprintf('ip=%i i0=%i, j0=%i\n',ip,i0,j0); 
%  if ip==23235;
%    keyboard;
%  end

  if j0-di<=0 | i0-di<=0 | j0+di>jdm | i0+di>idm
    continue; % near the bndry
  end

  A=dvu(j0-di:j0+di,i0-di:i0+di);
  inan = find(isnan(A));
  A(inan)=nanmean(nanmean(A));
  
  af = sum(sum(A.*Fw));
  dvuf(j0,i0)=af;
end






return
