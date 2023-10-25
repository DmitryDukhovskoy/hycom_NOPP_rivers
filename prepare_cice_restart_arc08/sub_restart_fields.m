function A = sub_restart_fields(fld,IMSK);
%
% Prepare CICE restart fields for writing
% state variables (ice/snow thicknesses, 
% enthalpie, etc. by catogires)
% are prepared from PIOMAS reanalysis
% for 1993 Jan. 1: 
% remap_piomas2arc.m
% IMSK is not needed for cases 'iceumask'
%
% Other fields are derived from HYCOM CICE 
% test simulation - form a template CICE restart file
% read_cice_restart.m
YR = 1993; 
MM = 1;
PTH.mat = '/nexsan/people/ddmitry/hycom/ARCc0.08/103/ice_restart/';
fmat0   = 'cice_restart093j'; % template for saved rest. fields 
                            % from existing CICE restart

% sea ice mask:			   
%fmat = sprintf('%shycom_ice_mask_%4.4i%2.2i.mat',...
%		 PTH.mat,YR,MM);
%load(fmat);
			    

switch(fld),
 case ('ivel'); % ice velocity
  fmat = sprintf('%s%s_ivel.mat',PTH.mat,fmat0);
  fprintf('Loading %s\n',fmat);
  load(fmat);
%  um=nanmean(nanmean(U));
%  vm=nanmean(nanmean(V));
  U(U==0)=1e-5; % fill missing values btw HYCOM and PIOMAS ice
  V(V==0)=1e-5;
  A.Title = 'Ice Velocity Components, m/s';
  A.U=U.*IMSK; % correct to match PIOMAS ice mask
  A.V=V.*IMSK;
%keyboard
 case ('rad'); % radiation
  fmat = sprintf('%s%s_rad.mat',PTH.mat,fmat0);
  fprintf('\nLoading %s\n',fmat);
  load(fmat);
  A.Title1 = 'Radiation: (1) scale factor, ';
  A.Title2 = '(2) sh/wv down vis.dir, (3) sh/wv down vis diff,';
  A.Title3 = '(4) sh/wv down near IR direct,';
  A.Title4 = '(5) sh/wv down near IR diff, W/m2';
  A.scale_factor = scale_fct; % no ice mask is needed
  A.swvdr        = swvdr;
  A.swvdf        = swvdf;
  A.swidr        = swidr;
  A.swidf        = swidf;

 case ('ostress');  % ocean stress
  fmat = sprintf('%s%s_ostr.mat',PTH.mat,fmat0);
  fprintf('\nLoading %s\n',fmat);
  load(fmat);
  A.Title    = 'Ocean stress, N/m2';
  A.strocnxT = strocnxT.*IMSK;
  A.strocnyT = strocnyT.*IMSK;

 case ('instress');  % internal stress
  fmat = sprintf('%s%s_intstrp.mat',PTH.mat,fmat0);
  fprintf('\nLoading %s\n',fmat);
  load(fmat);
  A.Title    = 'Internal stress tensor, kg/s2';
  A.stressp_1 = stressp_1.*IMSK;
  A.stressp_3 = stressp_3.*IMSK;
  A.stressp_2 = stressp_2.*IMSK;
  A.stressp_4 = stressp_4.*IMSK;
  
  fmat = sprintf('%s%s_intstrm.mat',PTH.mat,fmat0);
  fprintf('\nLoading %s\n',fmat);
  load(fmat);
  A.stressm_1 = stressm_1.*IMSK;
  A.stressm_3 = stressm_3.*IMSK;
  A.stressm_2 = stressm_2.*IMSK;
  A.stressm_4 = stressm_4.*IMSK;
  
  fmat = sprintf('%s%s_intstr12.mat',PTH.mat,fmat0);
  fprintf('\nLoading %s\n',fmat);
  load(fmat);
  A.stress12_1 = stress12_1.*IMSK;
  A.stress12_3 = stress12_3.*IMSK;
  A.stress12_2 = stress12_2.*IMSK;
  A.stress12_4 = stress12_4.*IMSK;
  
 case ('iceumask');
% icemask created from PIOMAS and HYCOM topo
% see: remap_piomas2arc.m
%  fmat = sprintf('%s%s_iceumask.mat',PTH.mat,fmat0);
%  fprintf('Loading %s\n',fmat);
%  load(fmat);
  fmat = sprintf('%shycom_ice_mask_%4.4i%2.2i.mat',...
		 PTH.mat,YR,MM);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  A.Title = 'ice mask derived from PIOMAS';
  A.IMSK = IMSK;
end



return
