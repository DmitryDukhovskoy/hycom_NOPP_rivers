function SCT = sub_set_Davis_xsct(HH,LON,LAT);
% Davis Str.
% 
%

IJ = [458     658
      556     658];

cca = length(IJ);

SCT.End_Pnts = IJ;

IIs = [];
JJs = [];
for ii=1:cca-1
  i1 = IJ(ii,1);
  i2 = IJ(ii+1,1);
  j1 = IJ(ii,2);
  j2 = IJ(ii+1,2);
  [I,J] = sub_xsct_indx(i1,j1,i2,j2);
  I = I(:);
  J = J(:);
  
  IIs = [IIs;I];
  JJs = [JJs;J];
end

ni = length(IIs);

SCT.II=IIs;
SCT.JJ=JJs;

% Find unit normal pointed in the Amerasian Basin
% norm is in the middle of a segment, i.e. either nrm_y=0 or nrm_x=0
clear Nrm
for ip=1:ni-1
  i0 = SCT.II(ip);
  j0 = SCT.JJ(ip);
  i1 = SCT.II(ip+1);
  j1 = SCT.JJ(ip+1);

  if i1==i0 % Y section, norm is along X axis
    Nrm(ip,1) = -1;
    Nrm(ip,2) = 0;
  else
    Nrm(ip,1) = 0;
    Nrm(ip,2) = 1;
  end
  
end
SCT.Norm = Nrm;
    
  
% Compute distances along the contour
% Get bottom depth
clear dx Hb
dx(1) = 0;
for ip=1:ni-1
  i0 = SCT.II(ip);
  j0 = SCT.JJ(ip);
  i1 = SCT.II(ip+1);
  j1 = SCT.JJ(ip+1);
  x0 = LON(j0,i0);
  y0 = LAT(j0,i0);
  x1 = LON(j1,i1);
  y1 = LAT(j1,i1);
  dx(ip+1) = distance_spheric_coord(y1,x1,y0,x0);
  Hb(ip) = HH(j0,i0);
end
Hb(ni)=Hb(1);
dst = cumsum(dx);

SCT.Distance_m = dst;
SCT.Hbottom    = Hb;


return