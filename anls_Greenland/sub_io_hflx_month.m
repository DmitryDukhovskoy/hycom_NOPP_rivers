function HFLX = sub_io_hflx_month(HFLX, s_mat, mold, imo, fmat, nlr)
% Process monthly fields
% and save 
%
%
if s_mat>0 & mold>0
  nrec = HFLX(mold).nrec;
  if nrec == 0 | isempty(nrec),
    error('# of saved records = 0');
  else
    fprintf('Calculating monthly mean, # of av. rcrds=%i\n',nrec);
  end

  dmm                     = HFLX(mold).Vol_flux_m3s;
  HFLX(mold).Vol_flux_m3s = dmm/nrec;
  dmm                     = HFLX(mold).Hflux_W;
  HFLX(mold).Hflux_W       = dmm/nrec;
  dmm                     = HFLX(mold).ZZ;
  HFLX(mold).ZZ            = dmm/nrec;
  
  dmm = HFLX(mold).Hflux_W;
  fprintf('   Mean Monthly HFlux = %6.4d W\n',nansum(nansum(dmm)));
%keyboard
  fprintf('=======    Saving %s\n\n',fmat);
  save(fmat,'HFLX');
else
  fprintf(' Initializing HFLX struct array...\n');
end

npp = length(HFLX(1).GrCntr_II);
if isempty(npp), 
  fprintf('Empty # of points npp=%i\n',npp);
  keyboard;
end

HFLX(imo).nrec = 0;
HFLX(imo).Vol_flux_m3s = zeros(1,npp);
HFLX(imo).Hflux_W = zeros(nlr,npp);
HFLX(imo).ZZ      = zeros(nlr,npp);

return