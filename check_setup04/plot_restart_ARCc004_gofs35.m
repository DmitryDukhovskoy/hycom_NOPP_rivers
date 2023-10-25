% Plot hycom restart fields
% for ARCc0.04 GOFS3.5
%
% One set of initial fields is prepared from my old AO HYCOM-CICE0.04 experiment 01.2
% another is from the Global 0.04 GLBc HYCOM-CICE5, GOFS3.5 
%
% Read all fields and check min/max with *b
% Check Land mask for sal temp
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

hg = 1e20;

%frst = '0.04 AO HYCOM-CICE4-012';
frst = 'GLBc0.04-GOFS3.5';

%pthbin  = '/Net/mars/ddmitry/hycom/ARCc0.04/restart/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
%fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);
if strncmp(frst,'GLBc0.04',8)
  pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.04/restart_022/';  % expt 02.2 - HYCOM-CICE5, GOFS3.5
else
  pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.04/restart_012/';  % expt 01.2 - HYCOM-CICE4, GOFS3.1
end
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

if strncmp(frst,'GLBc0.04',8)
  fina = sprintf('%sGLBc2ARCc004_restart_117a.a',pthbin);
  finb = sprintf('%sGLBc2ARCc004_restart_117a.b',pthbin);
else
  fina = sprintf('%srestart_from012_008a_116a.a',pthbin);
  finb = sprintf('%srestart_from012_008a_116a.b',pthbin);
end

fprintf('Plotting restart %s %s\n',frst,fina);

ftot = 1; % =0 - plot baroclinic u,v spd; =1 - plot total u,v, spd
%fld = 'v';
%fld = 'u';
%fld = 'spd'; 
%fld = 'dp';
fld = 'temp';
%fld = 'saln';
rlr = 15;
fbb = []; 
if ~strncmp(fld,'spd',3)
  F = read_hycom_restart(fina,finb,fld,IDM,JDM,'r_layer',rlr);

  if (strncmp(fld,'u',1) | strncmp(fld,'v',1)) & ftot==1
    fbb = sprintf('%sbavg',fld);
    B = read_hycom_restart(fina,finb,fbb,IDM,JDM);
    F = F+B;
  end 

else
  U = read_hycom_restart(fina,finb,'u',IDM,JDM,'r_layer',rlr);
  if ftot==1
    B = read_hycom_restart(fina,finb,'ubavg',IDM,JDM);
    U = B+U;
    fbb = 'ubavg';
  end

  V = read_hycom_restart(fina,finb,'v',IDM,JDM,'r_layer',rlr);
  if ftot==1
    B = read_hycom_restart(fina,finb,'vbavg',IDM,JDM);
    V = B+V;
  end
  F = sqrt(U.^2+V.^2);
  clear U V
end

		% Check land point:
IL=find(HH<0);
JL=find(F(IL)>1e20);
if ~isempty(JL), 
  fprintf('Field huge values are in the ocean points !!!! \n');
  [Jn,In]=ind2sub(size(HH),JL);
end;

I=find(F>1e20);
F(I)=nan;




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


btx = 'plot_restart_ARCc004_gofs35.m';
bottom_text(btx,'pwd',1);





