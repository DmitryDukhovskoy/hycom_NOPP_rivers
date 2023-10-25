function TRI = sub_layer_avrg_trConc(dH,Tr,HH,Acell);
% Vaerage tracer conc by layers
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
Tvrt = zeros(l,1); % mean mass in each layer
dZmn = zeros(1,1); % mean layer thkns
for k=1:l
  fprintf('Integrating over vert. layers: k=%i\n',k);
  tr = squeeze(Tr(k,:,:)); % kg/m3
  dz = squeeze(dH(k,:,:)); % m
  t1 = tr(IN).*dz(IN).*Acell(IN); % kg
  v1 = dz(IN).*Acell(IN); % m3
%  Mtr(IN) = Mtr(IN)+tr(IN).*dz(IN); % kg/m2 - mass per 1 m2
  Tvrt(k,1) = nansum(t1)./nansum(v1); % mean conc in layer, kg/m3
  dZmn(k,1) = nanmean(nanmean(dz(IN))); 
end

TRI.Vert_ConcMean_kg_m3 = Tvrt;
TRI.Mean_LThkn_m = dZmn;

return