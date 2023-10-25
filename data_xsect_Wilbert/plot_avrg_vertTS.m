% Check saved mean T/S xsections

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

% Specify temp/salin and 0.04/0.08
pfld  = 'temp';
%pfld  = 'salin';  
ires = 0.08;  % 0.04 or 0.08

fld0  = pfld;
YR1 = 2017;
YR2 = 2020;
regn = sprintf('ARCc%3.2f',ires);
if ires == 0.08
  expt = 123;
else
  expt = 23;
end
expt_nm = 'dltEddTmltEAPJRA';

pthout = '/nexsan/people/ddmitry/hycom/data_extracted/';
fmat_out = sprintf('%shycom%3.3i_%3.3i_vsect_%s_%i-%i.mat',...
                    pthout,ires*100,expt,pfld,YR1,YR2);

fprintf('Loading %s\n',fmat_out);
load(fmat_out);

TM  = XSCT.TM;
XL  = XSCT.Sect_dist;
ZZi = XSCT.Depths;
AAi = XSCT.Mean_Field;
Hb  = XSCT.Bott_depth;
IIs = XSCT.Sect_I;
JJs = XSCT.Sect_J;
Xl  = XSCT.Sect_lon;
Yl  = XSCT.Sect_lat;


dayS = TM(1);
dayE = TM(end);
fprintf('Start %s, End %s\n',datestr(dayS),datestr(dayE));



stl=sprintf('%4.2f, %s %s, %i-%i, 170W-10E ',...
     ires,expt_nm,pfld,YR1,YR2);
fprintf('Plotting %s\n',stl);
nf=1;
xdr = 1;
sub_plot_xsectZ(nf,XL,ZZi,AAi,stl,fld0,xdr,Hb);


