function dL = sub_segm_dL(IIs,JJs,LON,LAT);
% calculate lengths of the segments [st end] pnts
% 
ni=length(IIs(:,1));
for isg=1:ni
  i0=floor(IIs(isg,1));
  i1=floor(IIs(isg,2));
  j0=floor(JJs(isg,1));
  j1=floor(JJs(isg,2));
  dl=distance_spheric_coord(LAT(j0,i0),LON(j0,i0),LAT(j1,i1),LON(j1,i1));
  dL(isg,1)=dl;
end

return
