% Combine 2 time series for fluxes via straits/segments
% expt1 - main time series, uses expt 2 to fill gaps (start or end segments)
% Length of time series = min(expt1, expt2) - max(expt1,expt2)
%
function FLX = sub_combine2flx(ix1,ix2,EXPT,YR);

expt1   = EXPT(ix1).Nmb;
texpt1  = EXPT(ix1).cice_opt; % CICE options for sens. experiments
expt2   = EXPT(ix2).Nmb;
texpt2  = EXPT(ix2).cice_opt; % CICE options for sens. experiments

pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/strait_fluxes/',expt1);
fmatout=sprintf('%shycom004_%3.3i_%s_StraitFluxesDayMean_%4.4i.mat',...
									pthmat,expt1,texpt1,YR);
fprintf('Loading %s\n',fmatout);
load(fmatout);
F1 = SCT;

pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/strait_fluxes/',expt2);
fmatout=sprintf('%shycom004_%3.3i_%s_StraitFluxesDayMean_%4.4i.mat',...
									pthmat,expt2,texpt2,YR);
fprintf('Loading %s\n',fmatout);
load(fmatout);
F2 = SCT;

nsct = length(SCT);
clear SCT

for ik=1:nsct
  TM1 = F1(ik).Time;
  TM2 = F2(ik).Time;

  I2 = [];
  if TM1(1)>TM2(1)
    i2e = max(find(TM2<TM1(1)));
    TM = TM2(1:i2e);
    TM = [TM; TM1];
    I2 = [1:i2e]';
  else
    TM = TM1;
  end

  if TM2(end)>TM1(end)
    i2s = min(find(TM2>TM1(end)));
    TM = [TM;[TM2(i2s:end)]'];
    I2 = [I2; [i2s:length(TM2)]'];
  end

  nl = length(TM1);
  c1=0;
  clear I1
  for ill=1:nl
    i1 = find(TM==TM1(ill));
    I1(ill) = i1;
  end

  clear J2
  for imm=1:length(I2)
    i2 = find(TM==TM2(I2(imm)));
    J2(imm,1) = i2;
  end
%keyboard 

  FLX(ik).Name = F1(ik).Name;
  FLX(ik).Vol(I1) = F1(ik).VolFlx_m3s;
  FLX(ik).Vol(J2) = F2(ik).VolFlx_m3s(I2);
  FLX(ik).Vol1 = F1(ik).VolFlx_m3s;
  FLX(ik).Vol2 = F2(ik).VolFlx_m3s;
  FLX(ik).TM = TM;
  FLX(ik).TM1 = TM1;
  FLX(ik).TM2 = TM2;
end




return
