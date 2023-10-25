% Calculate mean fluxes across sections in array SCT
% for all sections
% COnvert to GFWA km3/s
function sub_GFWA_Flx(SCT);

nsct = length(SCT);
for isc=1:nsct
  nm=SCT(isc).Name;
  Nrm=SCT(isct).Norm;
  Xn=Nrm(:,1);
  Yn=Nrm(:,2);



return
