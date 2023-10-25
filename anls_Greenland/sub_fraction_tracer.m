function rr = sub_fraction_tracer(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ);
%
% Calculate fraction of tracer at every grid point (rr)
% to the total tracer mass in the whole region
%
% from 
% monthly-mean, depth-integrated 
% tracer concentrations 
% The Data are extracted in extr_trcr_mnth.m  <---- do not use, saves mean tr
% use: extr_MassTrcr_month.m <-- saves depth integrated mass by specified layers

pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
DV = datevec(dnmb);

hmsk=HH;
hmsk(HH<0)=nan;

iyr = DV(1,1);
imo = DV(1,2);
IN  = find(HH<0);

cc=0;
clear sumT
fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
%fmat = sprintf('%sMassTr_lrs_%i%2.2i.mat',pthmat,iyr,imo);
fprintf('Loading %s\n',fmat);
if exist(fmat,'file')
  load(fmat);
else
  fprintf(' =========  MISSING %s\n',fmat);
  return
end
rr = [];
fprintf('Reading: %i/%2.2i, Tracer %i, Lev %i\n',...
	iyr,imo,nTr,ilv);

Tr = squeeze(TRCR(ilv).TR(nTr,:,:));
Tr(Tr<=1e-23)=nan;
Tr_dom=squeeze(TRCR(3).TR(nTr,:,:)); % tracer integrated ove whole water depth
% Very rough estimate of dz, use hZ(x,y) instead
%dz = abs(TRCR(ilv).depth_av(2)-TRCR(ilv).depth_av(1));
%dz = hZ;
vol = Acell.*hZ; % grid cell vol
%intgrTr_dom = nansum(nansum(Tr.*Acell.*hZ));
intgrTr_dom = nansum(nansum(Tr.*Acell.*hZ));
Mtr = Tr.*vol; % mass, tracer in 1 layer
rr=Mtr/intgrTr_dom; 
keyboard



return

