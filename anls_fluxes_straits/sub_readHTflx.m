% Read Heat fluxes for specified exprt
%
function FLX = sub_readHTflx(fmatout);

fprintf('Loading %s\n',fmatout);
load(fmatout);
F1 = SCT;


nsct = length(SCT);
clear SCT

for ik=1:nsct
  TM1 = F1(ik).Time;
  nl = length(TM1);

  FLX(ik).Name  = F1(ik).Name;
  FLX(ik).Hflx1 = F1(ik).Hflx1_W;
  FLX(ik).Hflx2 = F1(ik).Hflx2_W;
  FLX(ik).TM    = TM1;
end

return
