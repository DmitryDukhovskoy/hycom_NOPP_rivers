% example how to read HYCOM 0.08 expt 11.0 data
% from OPeNDAP threads
% Dmitry Dukhovskoy, COAPS FSU

tic;
url='https://tds.hycom.org/thredds/dodsC/datasets/ARCc0.08/expt_11.0/data/2005_cice/';
flnm=sprintf('%s110_cice_inst.2005-01-09-00000.nc',url);
ncid = netcdf.open(flnm);
[numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);

% Lets look at the 6th variable (counting starts from 0!) that is sea ice thickness:
[name,xtype,dimids,natts] = netcdf.inqVar(ncid,5);

% Get data:
data = netcdf.getVar(ncid,5);

% Subset region you need (note indexing starts from 0!!!):
sbst = netcdf.getVar(ncid,5,[99, 89, 0],[1000, 1300, 1]); 
fprintf('Getting and subsetting data %8.6d sec\n',toc);



