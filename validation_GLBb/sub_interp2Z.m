function FI=sub_interp2Z(F,DP,ZM,ZZi,IPP);
% Interpolate fields 3D array F
% from HYCOM layers into
% fixed depths
% for all points IPP
% 
[ll,mm,nn]=size(F);
npb=length(IPP);
nli=length(ZZi);
FI=zeros(nli,mm,nn)*nan;
for ii=1:npb
  pmd=ii/npb*100;
  if mod(ii,10000)
    fprintf('Interp2Z: done %5.2f %% ...\n',pmd);
  end
  
  i0=IPP(ii);
  Ib=min(find(Dsec(:,i0)==0));
  if Ib==1, continue; end;
  ff=F(:,i0);
  zz=ZM(:,i0);
  hb=-sum(DP(:,i0));
  ibz = max(find(ZZi>=hb));
  for kl=Ib:nl
    zz(kl)=zz(kl-1)-0.1;
  end;
  zz=[0;zz];
  ff=[ff(1);ff];
  if zz(nl)>ZZi(end)
    zz(nl)=ZZi(end);
  end
  ffi = interp1(zz,ff,ZZi,'cubic');
  ffi(ibz+1:end)=nan;
  FI(:,i0)=ffi;
end


return