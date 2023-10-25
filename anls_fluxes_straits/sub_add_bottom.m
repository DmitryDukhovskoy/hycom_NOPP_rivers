function Tz = sub_add_bottom(Tz,Hb,ZZf);
% Add bottom - NaNs to the 2D vertical section
% Usually Needed for plotting
[nlr,np]=size(Tz);
%keyboard
for kk=1:np
  h0=Hb(kk);
  if h0>=0
    Tz(1:end,kk)=NaN;
  else
    izb=max(find(ZZf>=h0));
    if izb<nlr
      Tz(izb:end,kk)=nan;
    end
  end
end

return