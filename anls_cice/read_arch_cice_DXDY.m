% Read CICE output
% calculate DX, DY for grid centers and save

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt='110';
yr=1995;
mo=7;
dm=30;

%pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
pthbin  = '/nexsan/people/ddmitry/hycom/ARCc0.08/110/cice/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat = '/Net/ocean/ddmitry/HYCOM/ARCc/data_mat/';

fin = sprintf('%s%s_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',...
	      pthbin,expt,yr,mo,dm);

TLON = nc_varget(fin,'TLON');
TLAT = nc_varget(fin,'TLAT');

[DX,DY] = sub_dx_dy(TLON,TLAT);
fmt = sprintf('%sarcc0.08_CICEv4_dx_dy_tgrid.mat',pthmat);
fprintf('saving %s\n',fmt);
save(fmt,'DX','DY');



