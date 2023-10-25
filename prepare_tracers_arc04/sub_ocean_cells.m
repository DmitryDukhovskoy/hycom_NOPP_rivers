      function IJ = sub_ocean_cells(Ntot,HH,i0,j0);
% 
% Find Ntot ocean grid cells around and including i0,j0
%

ncl=1;
di=0;
clear i;
IJ(ncl,1)=i0;
IJ(ncl,2)=j0;
while ncl<Ntot
  di=di+1;  % #gr pnts from i,j - search radius
  LL=2*di;
  ic0=i0+di;
  jc0=j0+di;
  ivct=-1+0*i;  % dir vector
  for kk=1:4*LL % grid pnts around i,j
  %	kk
    if HH(jc0,ic0)<0
      ncl=ncl+1;
      IJ(ncl,1)=ic0;
      IJ(ncl,2)=jc0;
    end
    if ncl==Ntot, break; end;
    if mod(kk-1,LL)==0 & kk>1,
      ivct=ivct*i;
    end
    ic0=ic0+real(ivct);
    jc0=jc0+imag(ivct);
  end
  if di>Ntot+1, 
    error('Could not find %i ocean pnts Runoff',Ntot);
  end
end
  




return