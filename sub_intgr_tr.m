function TRI = sub_intgr_tr(dH,Tr,HH,Acell);
% Several characteristics: 
% Integrated tracer over depths/area 
% HH - can be used as mask
%      tracer integrated where HH<0
%     make HH>0 where integration is not needed
%
%if isempty(IN),
IN=find(HH<0);
%end
%keyboard
[l,m,n]=size(dH);

% Depth-integrate:
Mtot = 0;
Vtot = 0;
Mtr = zeros(m,n); % mass per 1 m2, kg/m2
Tvrt = zeros(l,1); % mean mass in each layer
dZmn = zeros(1,1); % mean layer thkns
for k=1:l
  fprintf('Integrating over vert. layers: k=%i\n',k);
  tr = squeeze(Tr(k,:,:)); % kg/m3
  dz = squeeze(dH(k,:,:));
  t1 = tr(IN).*dz(IN).*Acell(IN); % kg
  v1 = nansum(dz(IN).*Acell(IN));
  Mtot = Mtot+nansum(t1);
  Vtot = Vtot+nansum(v1);
  Mtr(IN) = Mtr(IN)+tr(IN).*dz(IN); % kg/m2 - mass per 1 m2
  Tvrt(k,1) = nansum(t1); % mass in layers, kg
  dZmn(k,1) = nanmean(nanmean(dz(IN))); 
end

TRI.DpthIntgrMass_kg_m2 = Mtr; 
TRI.OverallMass_kg = Mtot;
TRI.OverallVol_m3  = Vtot;
TRI.Vert_LrMass_kg = Tvrt;
TRI.Mean_LThkn = dZmn;

return