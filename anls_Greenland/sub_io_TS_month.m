function TSGR = sub_io_hflx_month(TSGR, s_mat, mold, imo, fmat, nlr, Hs)
% Process monthly fields
% and save 
%
%
if mold>0
  nrec = TSGR(mold).nrec;
  if nrec == 0 | isempty(nrec),
    error('# of saved records = 0');
  else
    fprintf('Calculating monthly mean, # of av. rcrds=%i\n',nrec);
  end

  DZ              = TSGR(mold).DZ;
  i0              = find(DZ<1e-6);
  DZ(i0)          = nan;
  dmm             = TSGR(mold).T;
  TSGR(mold).T    = dmm./DZ;
  dmm             = TSGR(mold).S;
  TSGR(mold).S    = dmm./DZ;
  dmm             = TSGR(mold).ZZ;
  TSGR(mold).ZZ   = dmm/nrec;
  dmm             = TSGR(mold).DZ;
  TSGR(mold).DZ   = dmm./nrec;
  dzm             = dmm/nrec;
%
% Check for mid-depth nan
  [m,n]=size(DZ);
  for ik=1:n
    DZ=TSGR(mold).DZ;
    S=TSGR(mold).S; 
    T=TSGR(mold).T; 
    inn=max(find(S(:,ik)>0));
    if inn==m, continue; end;
    S(inn+1:end,ik)=S(inn,ik);
    
    I=find(S(:,ik)==0);
    if ~isempty(I) & abs(DZ(1,ik))>0;
      fprintf('Middle-column nan, ik=%i\n',ik);
      keyboard;
    end
  end

  dmm = TSGR(mold).T;
  fprintf('   Mean Monthly T C= %6.4f, min = %6.4f, max = %6.4f\n',...
	  nanmean(nanmean(dmm)),min(min(dmm)),max(max(dmm)));
  dmm = TSGR(mold).S;
  fprintf('   Mean Monthly S = %6.4f, min = %6.4f, max = %6.4f \n',...
	  nanmean(nanmean(dmm)),min(min(dmm)),max(max(dmm)));
  
  I=find(Hs<0);
  hbz=sum(dzm);
  dbtm=abs(abs(Hs(I))-abs(hbz(I)));
  fprintf('   Max Bottom Mismatch with dZ = %6.4d \n',max(dbtm));
  
  if dbtm>1.0
    fprintf('Check DZ - large mismatch with bottom depth !!!\n');
    error(' Quitting ....');
  end
  
  
%keyboard
  if s_mat>0
  fprintf('=======    Saving %s\n\n',fmat);
  save(fmat,'TSGR');
  end
else
  fprintf(' Initializing TSGR struct array...\n');
end

npp = length(TSGR(1).GrCntr_II);
if isempty(npp), 
  fprintf('Empty # of points npp=%i\n',npp);
  keyboard;
end

TSGR(imo).nrec = 0;
TSGR(imo).T    = zeros(nlr,npp);
TSGR(imo).S    = zeros(nlr,npp);
TSGR(imo).ZZ   = zeros(nlr,npp);
TSGR(imo).DZ   = zeros(nlr,npp);

return