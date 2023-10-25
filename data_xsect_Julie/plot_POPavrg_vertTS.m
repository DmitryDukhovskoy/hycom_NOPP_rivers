% Check saved mean T/S xsections

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

% Specify temp/salin and 0.04/0.08
%pfld  = 'temp';
pfld  = 'salin';  
ires = 0.04;  % 0.04 or 0.08

fld0  = pfld;
YR1 = 2017;
YR2 = 2020;  % if YR2>YR1 - average over these years
expt = 23;
expt_nm = 'dltEddTmltEAPJRA';

%xname = 'BerNatl'; % Bering - N.Atl
%xname = 'BeaufNatl'; % Beaufort Shelf - South Norway
%xname = 'LaptSea'; % Laptev Sea
xname = 'BeaufIcel';  % Beaufort Sea to Iceland via N Pole
%xname = 'BafnFram'; % section around Greenland from Baffin to Fram str 


pthout = '/nexsan/people/ddmitry/data_POP_Julie/';
icc = 0;
for YR=YR1:YR2
  for mo=1:12
    fmat_out = sprintf('%sCrossArctic_MonthlyAvg_%4.4i-%2.2i.mat',pthout,YR,mo);
    fprintf('Loading %s\n',fmat_out);
    XSCT = load(fmat_out);

    ZZi = -XSCT.z_t*1e-2;
    if strncmp(pfld,'temp',4)
      AAi = XSCT.avg_t;
    else
      AAi = XSCT.avg_s;
    end
    
    icc = icc+1;
    if icc==1
      AAsum = AAi;
    else
      AAsum = AAsum + AAi;
    end
  end
end
if icc>1
  AAi = AAsum/icc;
end

Xl  = XSCT.lon;
Yl  = XSCT.lat;

% Derive distance along the line
% Adjust POP segment distances to match HYCOM
% 1000-km segments - see: plot_sections_map.m
%Ltot = 4431.1512; % Total length of the section
load('pop_chck_pnts.mat');  % Chck points for asjuting POP
I=find(Xl>180);
Xl(I)=Xl(I)-360;
%XL = sub_distance_section(Xl,Yl,Ltot);
XL = sub_scale_POPdist(Xl, Yl, Chck);

% There are messed up grid points around the NP in POP
% Need to cut out these segments:
Xlo = Xl;
Ylo = Yl;
XLo = XL;
  
XL  = XL(:);
ii1 = 389;
ii2 = 407;
Xl  = [Xl(1:ii1); Xl(ii2:end)];
Yl  = [Yl(1:ii1); Yl(ii2:end)];
XL  = [XL(1:ii1); XL(ii2:end)];
AAi = [AAi(:,1:ii1),AAi(:,ii2:end)];


% Find bottom:
[ll,nn]=size(AAi);

for ii=1:nn
  dmm = AAi(:,ii);
  k = min(find(isnan(dmm)));
  if k>1
    Zbt(ii,1) = 0.5*(ZZi(k)+ZZi(k-1));
  else
    Zbt(ii,1) = 0.;
  end

end


stl=sprintf('POP %s %i-%i',pfld,YR1, YR2);
fprintf('Plotting %s\n',stl);
nf=1;
xdr = 1;
sub_plot_xsectZ(nf,XL,ZZi,AAi,stl,fld0,xdr,Zbt);
%sub_plot_xsectZZv2(nf,XL,ZZi,AAi,stl,fld0,xdr,Hb);

btx = 'plot_POPavrg_vertTS.m';
bottom_text(btx,'pwd',1);

