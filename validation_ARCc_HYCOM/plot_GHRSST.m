% NASA JPL GHRSST fields
% https://podaac.jpl.nasa.gov/GHRSST?sections=about%2Bdata
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;



url0 = 'https://podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/GDS2/L4/GLOB/UKMO/OSTIA/v2/';



dnmb = datenum(2016,08,28);
dv   = datevec(dnmb);
jDay = dnmb-datenum(dv(1),1,1)+1; 

url = sprintf('%s%4.4i/%3.3i/',url0,dv(1),jDay);
fnmend = '-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc';
fnm = sprintf('%s%4.4i%2.2i%2.2i120000%s',url,dv(1:3),fnmend);

% Get fields:
% ncdisp(fnm)

ln1=-120;
ln2=60;
lt1=55;
lt2=90;
if ~exist('lont','var');
  lont = ncread(fnm,'lon');
  D=abs(lont-ln1);
  ix1=find(D==min(D),1);
  D=abs(lont-ln2);
  ix2=find(D==min(D),1);

  latt = ncread(fnm,'lat');
  D=abs(latt-lt1);
  iy1=find(D==min(D),1);
  D=abs(latt-lt2);
  iy2=find(D==min(D),1);

end

dx=ix2-ix1+1;
dy=iy2-iy1+1;
dmm = ncread(fnm,'analysed_sst');
dmm = dmm'-273.15;
stt = dmm(iy1:iy2,ix1:ix2);


  c1=-2;
  c2=12;

  nint = 200;
  CMP = colormap_sclr2(nint,c1,c2);
  cmp = CMP.colormap;
  cnt = CMP.intervals;
  nint=length(cmp);


figure(1); clf;
pcolor(stt); shading flat;
colormap(cmp);
caxis([c1 c2]);
%axis('equal');
set(gca,'xtick',[],...
        'ytick',[],...
        'xlim',[500 3000],...
        'ylim',[1 600],...
        'Color',[0 0 0]);

stl = sprintf('0.05 L4 GHRSST, daily, %i/%2.2i/%2.2i',dv(1:3));
title(stl);

  hps = [0.1 0.06 0.8 0.02];

  hb = colorbar;
  if hps(3)>hps(4),
    cloc = 'southoutside';
  else
    cloc = 'eastoutside';
  end
  set(hb,'Location',cloc,...
  'position',hps,...
  'TickLength',0.02,...
  'Fontsize',14);

btx = 'plot_GHRSST.m';
bottom_text(btx,'pwd',1);


