function [YRm,VFWm]=sub_read_myers;
% FWC in the domain from G> Myers
pthin='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';
flnm=sprintf('%sspg_gfwa_Myers.mat',pthin);

%fid=fopen(flnm,'r');
A=load(flnm);
na=size(A,1);

YRm  = A.yr2;
VFWm = A.fwa;

return
