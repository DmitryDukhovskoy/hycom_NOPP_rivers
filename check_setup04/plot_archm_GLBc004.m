% Plot fields from archm GLBc0.04
% used for creating restart or nest fields
% for ARCc0.04
%
% GLBc fields are converted to ARCc for plotting
% using 
% 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

hg = 1e20;
pthbin = '/Net/kronos/ddmitry/hycom/GOFS3.5/';

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);

fprintf('Reading topo: %s\n',fltopo);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn]= size(HH);
[m,n]= size(HH);
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
Ib = find(HH<0);

fina = sprintf('%sarchv_GLBc2arc04.2017_001_12.a',pthbin);
finb = sprintf('%sarchv_GLBc2arc04.2017_001_12.b',pthbin);


ftot = 1; % =0 - plot baroclinic u,v spd; =1 - plot total u,v, spd
%fld = 'v';
%fld = 'u';
%fld = 'spd'; 
%fld = 'dp';
fld = 'temp';
%fld = 'salin';
rlr = 15;
fbb = [];
if ~strncmp(fld,'spd',3)
  F = read_hycom(fina,finb,fld,IDM,JDM,'r_layer',rlr);
  F = squeeze(F);
  F(F>1e20)=nan;

  if (strncmp(fld,'u-vel.',3) | strncmp(fld,'v-vel.',3)) & ftot==1
    fbb = sprintf('%s_btrop',fld);
    B = read_hycom(fina,finb,fbb,IDM,JDM);
    B = squeeze(B);
    B(B>1e20)=nan;
    F = F+B;
  end

else
  U = read_hycom(fina,finb,'u-vel.',IDM,JDM,'r_layer',rlr);
  U = squeeze(U);
  if ftot==1
    B = read_hycom(fina,finb,'u_btrop',IDM,JDM);
    B = squeeze(B);
    B(B>1e20)=nan;

    U = B+U;
  end

  V = read_hycom(fina,finb,'v-vel.',IDM,JDM,'r_layer',rlr);
  V = squeeze(V); 
  if ftot==1
    B = read_hycom(fina,finb,'v_btrop',IDM,JDM);
    B = squeeze(B);
    B(B>1e20)=nan;

    V = B+V;
  end
  F = sqrt(U.^2+V.^2);
  clear U V
end


if strncmp(fld,'u',1) | strncmp(fld,'v',1)
  cl1 = flipud(colormap_blue(200));
  cl2 = colormap_red(200);
  cl1(end,:) = [1 1 1];
  cl2(1,:) = [1 1 1];
%  for j=0:5
%    cl1(end-j,:)=[1 1 1];
%  end
  cmp = [cl1;cl2];
else
  cmp = colormap('parula');
%  if strncmp(fld,'spd',3)
%    cmp(1,:)=[1 1 1];
%  end

end;


figure(1); clf;
pcolor(F); shading flat;
colormap(cmp);
hold on;
contour(HH,[0 0],'k');
axis('equal');
set(gca,'xlim',[1 nn],...
        'ylim',[1 mm]);

Is=max(strfind(fina,'/'));
if strncmp(fld,'u',1) | strncmp(fld,'v',1) | strncmp(fld,'spd',3)
  if ftot==0
    fldnm = sprintf('%s-bcln',fld);
  else
    fldnm = sprintf('%s(bcl+btrop)',fld);
  end
else
  fldnm=fld;
end

title(sprintf('%s, %s Layer %i',fina(Is+1:end), fldnm, rlr),'Interpreter','none');
colorbar

if strncmp(fld,'u',1) | strncmp(fld,'v',1)
  caxis([-0.5 0.5]);
elseif strncmp(fld,'spd',3);
  caxis([0 0.5]);
end


btx = 'plot_archm_GLBc004.m';
bottom_text(btx,'pwd',1);





