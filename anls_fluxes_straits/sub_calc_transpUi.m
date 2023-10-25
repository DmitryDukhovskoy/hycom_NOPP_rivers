% Calculate Vol Flux from depth-interpoalted
% U fields
function vf1 = sub_calc_transpUi(Unrm,dL,ZZ);

dZ = abs(diff(ZZ));
dz1 = 0.5*abs(ZZ(2)-ZZ(1));
dZ = [dz1;dZ];
dZ(end) = dZ(end)-dz1;
[DL,DZ]=meshgrid(dL,dZ);

vf1  = nansum(Unrm.*DL.*DZ);

return


