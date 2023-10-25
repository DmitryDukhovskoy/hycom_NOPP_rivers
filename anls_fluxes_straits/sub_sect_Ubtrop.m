function UTS=collocate_UTS_section(I,J,U,V,HH);
% Return, velocity component normal to the segment
% 1 layer only or barotropic field
%
%  HYCOM grid:
%
%          ------------   | V(i,j+1) ------------------
%          |                             |
%          |                             |
%        --- U(i,j)      *P(i,j)        ---  U(i+1,j)
%
%          |                             |
%          |                             |
%          ------------   | V(i,j) ------------------
%
%  Section follows U-V points
%  not centroids (P-points)
%   
HH(isnan(HH))=1e-10;
HH(HH>0)=1e-10;

nsc=length(I);
for is=1:nsc
  i0=I(is);
  j0=J(is);
  if is<nsc
    i1=I(is+1);
    j1=J(is+1);
  else
    i1=I(is-1);
    js=J(is-1);
  end
  
  if HH(j0,i0)>=0, 
    UTS.Hb(is)=HH(j0,i0);
    UTS.Unrm(:,is)=V(j0,i0);
    continue; 
  end;
  
  if j0==j1  % horizontal segment
    un=V(j0,i0);
    h1=HH(j0,i0);
    h2=HH(j0-1,i0);
  else
    un=U(j0,i0);
    h1=HH(j0,i0);
    h2=HH(j0,i0-1);    
  end
  

  Hn=abs(0.5*(h1+h2));
  
  UTS.Hb(is)=-Hn;
  UTS.Unrm(:,is)=un;
  
end



return