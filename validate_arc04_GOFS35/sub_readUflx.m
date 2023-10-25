% Read Vol flux for specified exprt
%
function FLX = sub_readUflx(fmatout);

fprintf('Loading %s\n',fmatout);
load(fmatout);
F1 = SCT;


nsct = length(SCT);
clear SCT

for ik=1:nsct
  TM1 = F1(ik).Time;
  nl = length(TM1);

  FLX(ik).Name = F1(ik).Name;
  FLX(ik).Vol  = F1(ik).VolFlx_m3s;
  FLX(ik).TM   = TM1;
end

return
