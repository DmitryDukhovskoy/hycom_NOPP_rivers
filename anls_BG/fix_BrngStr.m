% Fix  bug in Bering Str. FW flux
% mat files:
close all
clear

%pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_mat/';

crct = 40.9;

for iyr=2006:2006
  for imo=1:12
    fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
    fprintf('Loading %s\n',fmat);
    load(fmat);
    
    for il=1:2
%      for nTr=5:5
      nTr = 5;
      fprintf('Correcting, Year %i, Mo %2.2i, Lev: %i, Tracer #: %i\n',...
	      iyr,imo,il,nTr);
      dmm  = squeeze(TRCR(il).TR(nTr,:,:));
      cmax1 = max(max(dmm));
      dmm = dmm/crct;
      cmax2 = max(max(dmm));
      TRCR(il).TR(nTr,:,:)=dmm;
      fprintf('Max before crct %8.2f after %8.2f\n',cmax1, cmax2);
%      end
    end  % v. levels
% keyboard   
    fprintf('Saving %s\n',fmat);
    save(fmat,'TRCR');
   
  end % month
end   % year

