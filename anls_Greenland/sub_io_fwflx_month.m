function FWFLX = sub_io_fwflx_month(FWFLX, s_mat, mold, imo, fmat, nlr)
% Process monthly fields
% and save 
%
%
if s_mat>0 & mold>0
  nrec = FWFLX(mold).nrec;
  if nrec == 0 | isempty(nrec),
    error('# of saved records = 0');
  else
    fprintf('Calculating monthly mean, # of av. rcrds=%i\n',nrec);
  end

  dmm                      = FWFLX(mold).TrFlux_kg_s;
  FWFLX(mold).TrFlux_kg_s  = dmm/nrec;
  dmm                      = FWFLX(mold).ZZ;
  FWFLX(mold).ZZ           = dmm/nrec;
%keyboard

  dmm = FWFLX(mold).TrFlux_kg_s*1e-3;
  fprintf('   Mean Monthly Tracer Flux = %6.4d ton/s\n',nansum(nansum(dmm)));
%keyboard
  fprintf('=======    Saving %s\n\n',fmat);
  save(fmat,'FWFLX');
else
  fprintf('Initializing FWFlux struct array... \n');
end

npp = length(FWFLX(1).GrCntr_II);
if isempty(npp), 
  fprintf('Empty # of points npp=%i\n',npp);
  keyboard;
end

FWFLX(imo).nrec        = 0;
FWFLX(imo).TrFlux_kg_s = zeros(nlr,npp);
FWFLX(imo).ZZ          = zeros(nlr,npp);

return