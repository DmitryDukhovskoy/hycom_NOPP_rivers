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


pthout = '/nexsan/people/ddmitry/hycom/data_extracted/';
icc = 0;
for YR=YR1:YR2
  fmat_out = sprintf('%shycom%3.3i_%3.3i_xsect_%s_%s_%i-%i.mat',...
                    pthout,ires*100,expt,xname,pfld,YR,YR);

  fprintf('Loading %s\n',fmat_out);
  load(fmat_out);

  TM  = XSCT.TM;
  XL0 = XSCT.Sect_dist; % DIstance wrt to point 0 - gives non-monotonically increasing XL
  ZZi = XSCT.Depths;
  AAi = XSCT.Mean_Field;
  Hb  = XSCT.Bott_depth;
  IIs = XSCT.Sect_I;
  JJs = XSCT.Sect_J;
  Xl  = XSCT.Sect_lon;
  Yl  = XSCT.Sect_lat;
  
  icc = icc+1;
  if icc==1
    dayS = TM(1);
    AAsum = AAi;
  else
    AAsum = AAsum + AAi;
  end
end
if icc>1
  AAi = AAsum/icc;
end

% Derive distance along the line:
XLd = [];
nn=length(Xl);
x1=Xl(1);
y1=Yl(1);
for ii=1:nn
  x2=Xl(ii);
  y2=Yl(ii);
% Distance from point 0:
  dxx=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
%  XL(ii) = dxx;
  dxx=max([dxx,0.01]);
  dXX(ii)=dxx;
  x1=x2;
  y1=y2;
  XLd(ii)=sum(dXX(1:ii));
end

% Normalize and rescale to actual distance from point 0:
XL=XLd/max(XLd)*max(XL0);


dayE = TM(end);
fprintf('Start %s, End %s\n',datestr(dayS),datestr(dayE));

stl=sprintf('%4.2f, %s %s, %s, %i-%i',...
     ires, expt_nm, pfld, xname, YR1, YR2);
fprintf('Plotting %s\n',stl);
nf=1;
xdr = 1;
sub_plot_xsectZ(nf,XL,ZZi,AAi,stl,fld0,xdr,Hb);
%sub_plot_xsectZZv2(nf,XL,ZZi,AAi,stl,fld0,xdr,Hb);

btx = 'plot_avrg_vertTS.m';
bottom_text(btx,'pwd',1);

