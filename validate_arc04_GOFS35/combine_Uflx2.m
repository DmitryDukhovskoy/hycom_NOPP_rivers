% Fluxes derived in calc_vol_transp04.m
% Transports saved in separate mat files by years
function FLX = combine_Uflx2(pthmat,YR1,YR2,expt,segm);

Vflx = [];
DN  = [];
for YR=YR1:YR2
  fmat = sprintf('%s%3.3i_Vflux_%s_%4.4i_%4.4i.mat',...
         pthmat,expt,segm,YR,YR);
  fprintf('Loading %s\n',fmat);
  load(fmat);

  DN = [DN;TM];
  Vflx = [Vflx;FLXV];
end

FLX.TM = DN;
FLX.Vol = Vflx;
  

return
