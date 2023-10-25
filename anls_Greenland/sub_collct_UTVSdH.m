function CL = sub_collct_UTVSdH(U,V,T,S,dH,i1,j1,DX,DY,pnt);
%
% Also check spurious U,V near the "wall" cells
% i.e. where next grid is bottom
%
% Collocate data at pnt = 'vpnt' or 'pnt'
% for grid pnt i1,j1
%
%     ---------------------------
%     |            |             |    
%     |   T,S,dH   |    T,S,dH   |           
%   U -     i-1,  U-     i,j     |    
%     |      j     |             |    
%     |            |             |    
%     -------|------------|------
%            V            V

vpnt = strncmp(pnt,'vpnt',4);
upnt = strncmp(pnt,'upnt',4);

if ~vpnt & ~upnt
  error('Possible options for pnt vpnt or upnt');
end

u1    = squeeze(U(:,j1,i1));
u2    = squeeze(U(:,j1,i1+1));
hu_m1 = squeeze(dH(:,j1,i1-1));
hu0   = squeeze(dH(:,j1,i1));
hu_p1 = squeeze(dH(:,j1,i1+1));
dx1   = DX(j1,i1);
u1(hu_m1<1e-1)= nan;
u1(hu0<1e-1)  = nan;
u2(hu_p1<1e-1)= nan;
u2(hu0<1e-1)  = nan;
hu1   = 0.5*(hu_m1+hu0);
hu2   = 0.5*(hu0+hu_p1);
hu1(hu1<1e-1) = 0;
hu2(hu2<1e-1) = 0;

v1    = squeeze(V(:,j1,i1));
v2    = squeeze(V(:,j1+1,i1));
hv_m1 = squeeze(dH(:,j1-1,i1));
hv0   = squeeze(dH(:,j1,i1));
hv_p1 = squeeze(dH(:,j1+1,i1));
dy1   = squeeze(DY(j1,i1));
v1(hv_m1<1e-1)= nan;
v1(hv0<1e-1)  = nan;
v2(hv_p1<1e-1)= nan;
v2(hv0<1e-1)  = nan;
hv1   = 0.5*(hv_m1+hv0);
hv2   = 0.5*(hv0+hv_p1);
hv1(hu1<1e-1) = 0;
hv2(hu2<1e-1) = 0;

f_pr=0;
if f_pr==1
  fprintf('%s\n',fina);
  fprintf('grid pnt: i=%i, j=%i\n',i1,j1);
  fprintf('Depth: %8.4fm, sum(dP): %8.4fm\n',HH(j1,i1), sum(hv0));
  fprintf('Layer     v(i,j)     dP(i,j-1)     dP(i,j)\n');
  for k=1:41
    fprintf(' %2.2i:   %11.6f %11.6f %11.6f\n',k,v1(k),hv_m1(k),hv0(k));
  end
end;

u1 = sub_chck_uv(u1);
u2 = sub_chck_uv(u2);
v1 = sub_chck_uv(v1);
v2 = sub_chck_uv(v2);

fu1 = nansum(u1.*hu1*dy1);
fu2 = nansum(u2.*hu2*dy1);
fv1 = nansum(v1.*hv1*dx1);
fv2 = nansum(v2.*hv2*dx1);

divU = -fu1+fu2-fv1+fv2;
if abs(divU)>1.0
  dltU = -divU

if vpnt % interpolate into V(i,j)
  vv  = squeeze(V(:,j1,i1));
  dh1 = dH(:,j1-1,i1);
  dh2 = dH(:,j1,  i1);
  vv(dh1<1e-1)=nan;
  vv(dh2<1e-1)=nan;
  
  
else
  
  
end




return