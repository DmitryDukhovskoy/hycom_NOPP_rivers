function AAf = sub_fix_enthalpie(AA,ierr,ic,ilr,PTH,YY,MM);
% In the old version of remap_piomas2arc.m
% in the very thin upper layers, when Tsurf > Tice melt
% enthalpie got >0 (should be <0)
% was bug, rerun the whole code is too long
% try to fix several points
fprintf('Fixing enthalpie errors ...\n');

%fmat = sprintf('%srest_ice_aicen_%4.4i%2.2i.mat',...
%		PTH.mat,YY,MM);
%fprintf('Loading %s\n',fmat);
%load(fmat);

fmat = sprintf('%srest_ice_vicen_%4.4i%2.2i.mat',...
		PTH.mat,YY,MM);
fprintf('Loading %s\n',fmat);
load(fmat);

%fmat = sprintf('%srest_ice_vsnon_%4.4i%2.2i.mat',...
%		PTH.mat,YY,MM);
%fprintf('Loading %s\n',fmat);
%load(fmat);

fmat = sprintf('%srest_ice_trcrn_%4.4i%2.2i.mat',...
		PTH.mat,YY,MM);
fprintf('Loading %s\n',fmat);
load(fmat);

nfx = length(ierr);
nlr = 4;

Tocn = -1.8;
Tm   = -0.4; % ice melt T
AAf = AA;
for ik=1:nfx
  i0 = ierr(ik);
%  aicen = AICEN(i0,ic);
  vicen = VICEN(i0,ic);
  Tsrf = min([TSFC(i0,ic),Tm]);
  dT = Tsrf-Tocn;
  dTl= dT/nlr; % T change in 1 layer

  Tlr_btm = (ilr-1)*dTl+Tocn;
  Tlr_top = ilr*dTl+Tocn;
  Tlr = 0.5*(Tlr_btm+Tlr_top);
  Si = 3;
  qi = sub_enthalpie('ice',Tlr,Tm); % J/m3
  enth = qi*aicen*hicen/nlr; % J/m2 - enth. per unit area
  AAf(i0) = enth;
%keyboard
end

fprintf('All fixed ...\n');


return