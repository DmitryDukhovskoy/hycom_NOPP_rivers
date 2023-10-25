function ufx = sub_check_uv(u,flg);
% Quick check and fix of spurious U,V in the
% mean U,V fields HYCOM
% u is a 1D (depth) array
% better to have bottom values = nan if not
% the bottom values may be filled with the last 
% nonspurious U
uold = u;
nl = length(u);
ufx = u;
%u(u==0) = nan;
ib = max(find(~isnan(u)));
if ib==1,   return;  end
if max(abs(u))<2.; return; end;

rmx = 20; % assume rmx IQR as maximum difference 
rmv = 30;

I=0;
cx=0;
while ~isempty(I),
  
  du  = abs(diff(u));
  du(du==0) = nan;
  %du  = abs(u);
  mdu = nanmedian(du);
  iqdu= iqr(du); 
  dbu = abs(du-mdu)/iqdu; % measure in IQR distances
  I=find(dbu>rmx); 
  %keyboard
  
  if isempty(I); 
    break; 
  end;
% Check magnitude:
%  ua = abs(u);
  mdUa = nanmedian(u);
  iqUa = iqr(u);
  ul1  = mdUa - rmv*iqUa;
  ul2  = mdUa + rmv*iqUa;
  if min(u)>=ul1 & max(u)<=ul2 & max(abs(u))<2.0, 
    break; 
  end
%keyboard  
  cx=cx+1;
  if cx>nl 
    error('sub_chck_uv: endless loop ...');
  end

  if cx==1 & flg>0,
    fprintf('::: ==> Spurious U, maxU=%6.4f medU=%6.4f minU=%6.4f\n',...
	max(u), nanmedian(u), min(u));
  end
  
  ni=length(I);
  if I(1)>1
    for ik=1:ni
      k=I(ik);
      u(k+1)=u(k);
    end
  else
    u(1) = median(u);
    for ik=2:ni
      k=I(ik);
      u(k+1)=u(k);
    end
  end  
end
if cx>0 & flg>0
  fprintf('::: <== Spurious fixed U, maxU=%6.4f minU=%6.4f\n',...
	max(u), min(u));
%  keyboard
end

ufx = u;

return