function ufx = sub_check_uv(u,dh,flg);
% Quick check and fix of spurious U,V in the
% mean U,V fields HYCOM
% u is a 1D (depth) array
% better to have bottom values = nan if not
% the bottom values may be filled with the last 
% nonspurious U

rmx = 15; % assume rmx IQR as maximum difference 

uold = u;
nl = length(u);
ufx = u;
%u(u==0) = nan;
ib = max(find(~isnan(u)));
if ib==1,   return;  end
ib = find(abs(u)>1e-4);
if isempty(ib), return; end;

au  = abs(u);
mau = nanmedian(au);
iqau= iqr(au); 
if max(iqau)<1e-3, return; end;
dau = abs(au-mau)/iqau; % measure in IQR distances

%fprintf('max au=%8.6f\n',max(dau))
%keyboard
if max(dau)<rmx; return; end;

dz=dh;
dz(dz==0)=0.1;
zz=-cumsum(dz);
zz=[0;zz];
clear zm
for k=1:nl
  zm(k,1)=zz(k)-0.5*dh(k);
end



clear u0 zm0
cc=0;
dau(isnan(dau))=0;
for ik=1:nl
  if dau(ik)<rmx
    cc=cc+1;
    zm0(cc)=zm(ik);
    u0(cc)=u(ik);
  end
end
u0=u0(:);
u0(isnan(u0))=0;
Inan=find(isnan(u));
u(Inan)=0;
zm0=zm0(:);

ui = interp1(zm0,u0,zm,'pchip');
ufx = ui;

if max(diff(Inan)>1);
  fprintf('Check Inan - nans in the middle water column\n');
%  keyboard
end
ufx(Inan)=nan;

Utr=nansum(u.*dh);
Utrf=nansum(ufx.*dh);

if flg>0
  fprintf('::: ==> Spurious original U, maxU=%6.4f minU=%6.4f\n',...
	max(u), min(u));
  fprintf('::: <== Spurious fixed U, maxU=%6.4f minU=%6.4f\n',...
	max(ufx), min(ufx));
  fprintf(':::  Transport original: %6.2f m2/s, fixed: %6.2f m2/s\n',...
	  Utr,Utrf);
%  keyboard
end

%ufx = u;

return