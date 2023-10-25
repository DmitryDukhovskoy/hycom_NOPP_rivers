% Fix counter bug in mat files:
close all
clear

%pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_mat/';

for iyr=2006:2006
  for imo=1:12
    fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
    fprintf('Loading %s\n',fmat);
    load(fmat);
    
    for il=1:2
      nrec = TRCR(il).nrec;
      if nrec<30
	fprintf('nrec in mat file %i, already corrected?\n',nrec);
	error('Double check nrec ...');
      end
      
      nrec0 = nrec/5;
      if mod(nrec,5) >0
	error('nrec is not divisible by # of tracers ');
      else
	fprintf('Old nrec=%i, corrected nrec %i\n',nrec,nrec0);
      end
      
      for nTr=1:5
	fprintf('Correcting, Year %i, Mo %2.2i, Lev: %i, Tracer #: %i\n',...
		iyr,imo,il,nTr);
        dmm  = squeeze(TRCR(il).TR(nTr,:,:));
	dmm = dmm*nrec/nrec0;
	TRCR(il).TR(nTr,:,:)=dmm;
	TRCR(il).nrec = nrec0;
      end
    end  % v. levels
% keyboard   
    fprintf('Saving %s\n',fmat);
    save(fmat,'TRCR');
   
  end % month
end   % year

