% For Greenland tracer analysis
% select boxes
% Local Norm: indicates positive orientation of normal
%  [X, Y]

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.08';
expt = '110';
TV   = 11; % topo
pthfig  = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/%s/fig_domain/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
%pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/tmp_data/';

s_fig=0;
s_mat=0;
vrs = 4;  % Phase 2: added new sections to close all N. Atl., subp.gyre
txtb='ddmitry@mars: hycom_NOPP_rivers/anls_Greenland/set_sections.m';

ftopo = sprintf('%s/depth_%s_%i.nc',pthtopo,regn,TV); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

domname='whole';
%IND = smaller_domain_indices('Green');
%IND = smaller_domain_indices('NorthAtl');
%IND = smaller_domain_indices('subpolgyre');
IND = [];

if isempty(IND) % whole domain
  IND.i1 = 1;
  IND.i2 = nn;
  IND.j1 = 1;
  IND.j2 = mm;
  IND.dj = mm;
  IND.di = nn;
end

  
inc1=IND.i1;
inc2=IND.i2;
jnc1=IND.j1;
jnc2=IND.j2;
djnc=IND.dj;
dinc=IND.di;



HH=HH(jnc1:jnc2,inc1:inc2);
LON=LON(jnc1:jnc2,inc1:inc2);
LAT=LAT(jnc1:jnc2,inc1:inc2);
[mm,nn]=size(LON);
HH(1500:end,:)=nan;


ii=0;
fn=1;
sub_plot_bath(HH,LON,LAT,fn,domname);

% Denmark Strait:
is=812;
ie=854;
js=612;
je=552;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='Denmark';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,1]; % local normal, determines + flux [X,Y]

% Iceland-Faroe-Scotland 
is=951;
ie=1027;
js=512;
je=449;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='IcelandFaroe1';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,1]; % local normal, determines + flux

% Second hand:
is=1032;
ie=1164;
js=446;
je=446;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='IcelandFaroe2';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,1]; % local normal, determines + flux



% Barents Sea Opening
is=1142;
ie=1242;
js=887;
je=720;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='BarentsSO';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,-1]; % local normal, determines + flux

% Fram Strait
is=928;
ie=1071;
js=954;
je=954;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='FramStr';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,-1]; % local normal, determines + flux

% Baffin Bay:
% North
% Nares Strait
is=602;
ie=602;
js=998;
je=1052;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='Nares Strait';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,0]; % local normal, determines + flux

% Jones Sound
is=530;
ie=547;
js=1045;
je=1045;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='Jones Sound';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,-1]; % local normal, determines + flux


% Lancaster Sound
is=475;
ie=507;
js=1016;
je=1016;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='Lancaster Sound';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,-1]; % local normal, determines + flux


% Davis Str.
is=455;
ie=561;
js=661;
je=661;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='DavisStr';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,1]; % local normal, determines + flux

% Hudson Strait
is=376;
ie=376;
js=520;
je=575;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='HudsonStr';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,0]; % local normal, determines + flux


% Labrador Sea
is=434;
ie=605;
js=255;
je=420;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='Labrador Opening';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,1]; % local normal, determines + flux

% N. Atlantic s. bndry: N. Scotia - B. Biscay
is=362;
ie=1114;
js=90;
je=90;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'b.-');
I=sub2ind(size(HH),JJ,II);

ii=ii+1;
SEGM(ii).Name='SouthBndry Atl';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,1]; % local normal, determines + flux



% ----------------------
% Interior regions:
% Greenland Sea:
% ----------------------
% segm1
is=1014;
ie=1080;
js=811;
je=811;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='GreenlandSea_s1';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,-1]; % local normal, determines + flux

% Greenland Sea:
% segm2
is=ie;
ie=is;
js=je;
je=734;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='GreenlandSea_s2';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,0]; % local normal, determines + flux


% Greenland Sea:
% segm3
is=ie;
ie=1014;
js=je;
je=js;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='GreenlandSea_s3';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,1]; % local normal, determines + flux

% Greenland Sea:
% segm4
is=ie;
ie=is;
js=je;
je=811;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='GreenlandSea_s4';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,0]; % local normal, determines + flux


% ----------------------
% Icealand Sea
% ----------------------
% sgm1
is=922;
ie=922;
js=595;
je=669;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='IcelandSea_s1';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,0]; % local normal, determines + flux

% Iceland Sea
% sgm2
is=ie;
ie=995;
js=je;
je=js;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='IcelandSea_s2';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,-1]; % local normal, determines + flux

% Iceland Sea
% sgm3
is=ie;
ie=is;
js=je;
je=595;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='IcelandSea_s3';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,0]; % local normal, determines + flux

% Iceland Sea
% sgm4
is=ie;
ie=922;
js=je;
je=js;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='IcelandSea_s4';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,1]; % local normal, determines + flux

% ----------------------
% Labrador Sea
% ----------------------
% sgm1
is=460;
ie=532;
js=451;
je=js;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='LabradorSea_s1';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,-1]; % local normal, determines + flux

% sgm2
is=ie;
ie=is;
js=je;
je=363;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='LabradorSea_s2';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,0]; % local normal, determines + flux

% sgm3
is=ie;
ie=460;
js=je;
je=js;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='LabradorSea_s3';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,1]; % local normal, determines + flux

% sgm4
is=ie;
ie=is;
js=je;
je=451;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='LabradorSea_s4';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,0]; % local normal, determines + flux



% ----------------------
% Irminger Sea
% convection site: southern part, 
% Pickart et al., 2003, Deep-Sea Res.I
% ----------------------
% sgm1
is=647;
ie=719;
js=415;
je=js;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='IrmingerSea_s1';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,-1]; % local normal, determines + flux

% sgm2
is=ie;
ie=is;
js=je;
je=348;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='IrmingerSea_s2';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,0]; % local normal, determines + flux

% sgm3
is=ie;
ie=647;
js=je;
je=js;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='IrmingerSea_s3';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,1]; % local normal, determines + flux

% sgm4
is=ie;
ie=is;
js=je;
je=415;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'r.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='IrmingerSea_s4';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,0]; % local normal, determines + flux

ieddy=ii+1;
% ----------------------------------
%
% Isobaths - for eddy flux calculations:
%
% ----------------------------------
% Nordic Seas:
is=872;
ie=996;
js=598;
je=847;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Nordic_s1';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,-1]; % local normal, determines + flux

is=ie;
ie=1126;
js=je;
je=846;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Nordic_s2';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,-1]; % local normal, determines + flux

is=ie;
ie=1208;
js=je;
je=698;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Nordic_s3';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,-1]; % local normal, determines + flux

is=ie;
ie=1130;
js=je;
je=531;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Nordic_s4';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,1]; % local normal, determines + flux

is=ie;
ie=1002;
js=je;
je=js;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Nordic_s5';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,1]; % local normal, determines + flux

% Icelandic Current
is=ie;
ie=872;
js=je;
je=598;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Nordic_s6';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,1]; % local normal, determines + flux

% ---------------------------
% Eddy Segments - Irminger Sea
is=617;
ie=707;
js=398;
je=492;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Irminger_s1';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,-1]; % local normal, determines + flux

% Irminger current
is=ie;
ie=759;
js=je;
je=423;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Irminger_s2';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,-1]; % local normal, determines + flux

% E. Irminger Sea current
is=ie;
ie=722;
js=je;
je=299;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Irminger_s3';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,1]; % local normal, determines + flux

% S. Irminger Sea current
is=ie;
ie=598;
js=je;
je=306;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Irminger_s4';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,1]; % local normal, determines + flux

% SW. Irminger Sea 
is=ie;
ie=617;
js=je;
je=398;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Irminger_s5';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,0]; % local normal, determines + flux


% ----------------
% Labardor Sea eddy:
is=494;
ie=419;
js=305;
je=396;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Labr_s1';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,1]; % local normal, determines + flux

is=ie;
ie=454;
js=je;
je=487;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Labr_s2';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[1,0]; % local normal, determines + flux

is=ie;
ie=535;
js=je;
je=js;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Labr_s3';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[0,-1]; % local normal, determines + flux

is=ie;
ie=581;
js=je;
je=396;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Labr_s4';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,-1]; % local normal, determines + flux

is=ie;
ie=494;
js=je;
je=305;
[II,JJ]=sub_xsct_indx(is,js,ie,je);
plot(II,JJ,'m.-');
I=sub2ind(size(HH),JJ,II);
ii=ii+1;
SEGM(ii).Name='Eddy_Labr_s4';
SEGM(ii).I=II;
SEGM(ii).J=JJ;
SEGM(ii).LON=LON(I);
SEGM(ii).LAT=LAT(I);
SEGM(ii).Norm=[-1,1]; % local normal, determines + flux

bottom_text(txtb);
if s_fig==1
  fgnm=sprintf('%ssections_all_v%i',pthfig,vrs);
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r200',fgnm);
end



% -----------------------------
% Plot sections Lon/Lat:
% -----------------------------
ii=0;
fn=2;
sub_plot_bath(HH,LON,LAT,fn,domname);

ns=length(SEGM);
for ir=1:ieddy-1;
  nm=SEGM(ir).Name;
  i1=SEGM(ir).I(1);
  j1=SEGM(ir).J(1);
  i2=SEGM(ir).I(end);
  j2=SEGM(ir).J(end);
  ln1=SEGM(ir).LON(1);
  lt1=SEGM(ir).LAT(1);
  ln2=SEGM(ir).LON(end);
  lt2=SEGM(ir).LAT(end);

  SPP1{1,1}=sprintf('%7.3fE',ln1);
  SPP1{2,1}=sprintf('%7.3fN',lt1);
  SPP2{1,1}=sprintf('%7.3fE',ln2);
  SPP2{2,1}=sprintf('%7.3fN',lt2);

  II = SEGM(ir).I;
  JJ = SEGM(ir).J;

  is1=strfind(nm,'GreenlandSea');
  is2=strfind(nm,'IcelandSea');
  is3=strfind(nm,'LabradorSea');
  is4=strfind(nm,'IrmingerSea');
  
  text(i1-60,j1-10,SPP1,'FontSize',11);
  if isempty(is1) & isempty(is2) & isempty(is3) & isempty(is4)
    plot(II,JJ,'b.-');
    text(i2+5,j2-10,SPP2,'FontSize',11);
  else
    plot(II,JJ,'r.-');
  end
 
% Draw norms:
  nx=SEGM(ir).Norm(1);
  ny=SEGM(ir).Norm(2);
  mx=0.5*(i1+i2);
  my=0.5*(j1+j2);
  plot([mx mx+20*nx],[my my+20*ny],'g-','linewidth',2);
  
  stl=sprintf('Sections for budget calculations, vrs=%i',vrs);
  title(stl);
end;

bottom_text(txtb);
if s_fig==1
  fgnm=sprintf('%ssections_budget_v%i',pthfig,vrs);
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r200',fgnm);
end


% Plot segments for eddy fluxes:
fn=3;
sub_plot_bath(HH,LON,LAT,fn,domname);

nmo='mmmmmmmmmmmm';
for ir=ieddy:ns;
  nm=SEGM(ir).Name;
  i1=SEGM(ir).I(1);
  j1=SEGM(ir).J(1);
  i2=SEGM(ir).I(end);
  j2=SEGM(ir).J(end);
  ln1=SEGM(ir).LON(1);
  lt1=SEGM(ir).LAT(1);
  ln2=SEGM(ir).LON(end);
  lt2=SEGM(ir).LAT(end);

  II = SEGM(ir).I;
  JJ = SEGM(ir).J;

  plot(II,JJ,'m.');
  
  SPP1{1,1}=sprintf('%7.3fE',ln1);
  SPP1{2,1}=sprintf('%7.3fN',lt1);
  SPP2{1,1}=sprintf('%7.3fE',ln2);
  SPP2{2,1}=sprintf('%7.3fN',lt2);

  is1 = strncmp(nm,nmo,8);
  if is1,
    iss=iss+1;
  else
    nmo=nm;
    iss=1;
  end
  
  
  
  if 2*round(iss/2)~=iss 
    text(i1-70,j1-20,SPP1,'FontSize',11);
    text(i2-70,j2-20,SPP2,'FontSize',11);
  end
 
% Draw norms:
  nx=SEGM(ir).Norm(1);
  ny=SEGM(ir).Norm(2);
  mx=0.5*(i1+i2);
  my=0.5*(j1+j2);
  plot([mx mx+20*nx],[my my+20*ny],'g-','linewidth',2);
  
  stl=sprintf('Sections for eddy fluxes, vrs=%i',vrs);
  title(stl);
end;


bottom_text(txtb);
if s_fig==1
  fgnm=sprintf('%ssections_eddy_fluxes_v%i',pthfig,vrs);
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r200',fgnm);
end






if s_mat==1
  fmat=sprintf('%sSubarctic_sections_v%i',pthmat,vrs);
  fprintf('Saving %s\n',fmat);
  save(fmat,'SEGM','IND');
end

