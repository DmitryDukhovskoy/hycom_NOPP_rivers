function rr = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ,IntSrf);
%
% Calculate fraction of tracer (nTr=1,...,5) at every grid point (rr)
% to the total tracer mass in the whole region
%
% Note hZ is the layer thickness, scalar
% it works for cases when tracer is integrated exactly 
% within the depth limits
%
% IntSrf=1 - integrate from surface to level=ilv (its bottom interface)
% from 
% monthly-mean, depth-integrated 
% tracer concentrations 
% The Data are extracted in extr_trcr_mnth.m  <---- do not use, saves mean tr
% use: extr_MassTrcr_month.m <-- saves depth integrated mass by specified layers
%
% Tracer content calculated in vol_intgr_regn_trcr008.m
%
% Note: updated regions, now all regions
% are adjacent to each other
% So that Subpolar Gyre
% combines reg #1 (Labr)+Reg#2(IcelSea)+Reg#7 (Centr NAtl)

if expt==110
  pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_mat/';
else
  pthmat = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/%3.3i/data_matTr/',expt);
end
DV = datevec(dnmb);

hmsk=HH;
hmsk(HH<0)=nan;

iyr = DV(1,1);
imo = DV(1,2);
IN  = find(HH<0);

cc=0;
rr = [];
clear sumT
%fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
%fmat = sprintf('%sMassTr_lrs_%i%2.2i.mat',pthmat,iyr,imo);
fmat = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,iyr,imo);
fprintf('Loading %s\n',fmat);
if exist(fmat,'file')
  load(fmat);
else
  fprintf(' =========  MISSING %s\n',fmat);
  return
end
fprintf('Reading: %i/%2.2i, Tracer %i, Lev %i\n',...
	iyr,imo,nTr,ilv);

% find whole-depth layer
nlr = length(TRCR);
%if isfield('Layer_Z2',TRCR);
ibtm=0;
for ilr=1:nlr
  z2=TRCR(ilr).Layer_Z2;
  if z2<-9000, ibtm=ilr; end;
end

%keyboard
  
Tr = squeeze(TRCR(ilv).MassTr_kg);
if IntSrf==1 & ilv~=ibtm % integrate from surface to ilv
  Tr=Tr*0;
  for ikk=1:ilv
    pmm = squeeze(TRCR(ikk).MassTr_kg);
    Tr=Tr+pmm;
  end
end

Tr_dom=squeeze(TRCR(ibtm).MassTr_kg); % tracer integrated ove whole water depth
%Tr(Tr<=1e-18)=nan;
Tr_dom(Tr_dom<=1e7)=nan;
%keyboard
%vol = Acell*hZ; % grid cell vol
%intgrTr_dom = nansum(nansum(Tr.*Acell.*hZ));
MTr_dom = nansum(nansum(Tr_dom));
MTr_lr = Tr; % mass, tracer in 1 layer
rr=MTr_lr/MTr_dom; 
%keyboard

return

