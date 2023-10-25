function TFLX = sub_io_tflx_month(TFLX, s_mat, mold, imo, fmat, nlr)
% Process monthly fields
% and save 
%
%
if s_mat>0 & mold>0
  nrec = TFLX(mold).nrec;
  if nrec == 0 | isempty(nrec),
    error('# of saved records = 0');
  else
    fprintf('Calculating monthly mean, # of av. rcrds=%i\n',nrec);
  end

  dmm = TFLX(mold).VFlx_m3s;
  TFLX(mold).VFlx_m3s    = dmm/nrec;
  dmm = TFLX(mold).VFlxZ_m3s;
  TFLX(mold).VFlxZ_m3s   = dmm/nrec;
  dmm = TFLX(mold).FWflx_m3s;
  TFLX(mold).FWflx_m3s   = dmm/nrec;
  dmm = TFLX(mold).FWflxZ_m3s;
  TFLX(mold).FWflxZ_m3s  = dmm/nrec;
  for nTr=1:4
    dmm  = TFLX(mold).TrFlx_kgs(nTr,:);
    TFLX(mold).TrFlx_kgs(nTr,:)  = dmm/nrec;
    dmm = TFLX(mold).TrFlxZ_kgs(nTr,:);
    TFLX(mold).TrFlxZ_kgs(nTr,:)  = dmm/nrec;
  end
  
%keyboard
  fprintf('======= Saving: Month %i,   %s\n\n',mold,fmat);
  save(fmat,'TFLX');
end
if s_mat==0
  fprintf('!!!  Not saved mat file !!!\n\n ');
end

II                    = TFLX(1).Indx;
npp                   = length(II);
TFLX(imo).nrec        = 0;
TFLX(imo).VFlx_m3s    = zeros(1,npp); % VOl Flux
TFLX(imo).VFlxZ_m3s   = zeros(1,npp); % VOl Flux intgr to z0
TFLX(imo).FWflx_m3s   = zeros(1,npp); % Total FWFlx
TFLX(imo).FWflxZ_m3s  = zeros(1,npp); % FWFlx to depthz0
TFLX(imo).TrFlx_kgs   = zeros(4,npp); % Total Tracer Flx
TFLX(imo).TrFlxZ_kgs  = zeros(4,npp); % Tr Flx to depth z0

return