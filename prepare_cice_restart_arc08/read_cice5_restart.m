% Reads restart for CICE v5 
% netCDF file
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

YH = 117;
AB = 'f';

ncat = 5; 

pthin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/incoming/';
%restart created map_arc04_to_arc08.m:
finp = sprintf('%scice.restart08_%3.3i%s_DD.nc',pthin,YH,AB);
finp = sprintf('%scice.restart04_%3.3i%s.nc',pthin,YH,AB); % original 004 restart


% INPUT fields:
f_getnm=0;
if f_getnm==1
		fprintf('Initial file: %s\n',finp);
		finid = netcdf.open(finp,'NC_NOWRITE');
		[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(finid);
		for iv=1:nvars
				[varname vartype vardimIDs varatts] = netcdf.inqVar(finid,iv-1);
				fprintf('%i %s\n',iv,varname);
				BI(iv).Name = varname;
				BI(iv).Ndim = length(vardimIDs);
		end
		netcdf.close(finid);
end

%
% 2D fields:
vnm = 'aicen';
AA = nc_varget(finp,vnm);
ndim = length(size(AA));
if ndim==3
  AA = squeeze(sum(AA,1));
end

figure(1); clf;
pcolor(AA); shading flat;
axis('equal');
colorbar
title(finp,'Interpreter','none');

btx = 'read_cice5_restart.m';
bottom_text(btx,'pwd',1);

