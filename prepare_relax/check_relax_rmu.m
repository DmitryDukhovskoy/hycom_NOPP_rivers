addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

%R = 'ARCc0.04';
R = 'ARCc0.08';
E = '011';
ntopo1=11;
%ntopo2=17;
TV = sprintf('%2.2i',ntopo1);

ptharc  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthout  = sprintf('/Net/mars/ddmitry/hycom/%s/relax/%s/41layers_T%s/',R,E,TV);
pthout = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/110/41layers_T11/';
pthtopo = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);
pthin1  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/110/41layers_T11/';


% Get topo and grid:
fltopo=sprintf('%sdepth_ARCc0.08_%s.nc',pthtopo,TV);
fprintf('Reading topo: %s\n',fltopo);
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
[m,n]= size(HH);


fina = sprintf('%srelax_rmu_T%s.a',pthout,TV);
finb = sprintf('%srelax_rmu_T%s.b',pthout,TV);

fprintf('Reading %s\n',fina);
fprintf('Reading %s\n',finb);

fida = fopen(fina,'r','ieee-be');
fidb = fopen(finb,'r');

aa='a';
pat='[ij]dm';
while isstr(aa)
  aa = fgetl(fidb);
  if ~isstr(aa); break; end;
  ii=regexp(aa,pat);
  if ~isempty(ii),
    ind=findstr(aa,'=');
    ch=aa(ind+1:end);
    dmm=sscanf(ch,'%i');
    IDM=dmm(1);
    JDM=dmm(2);
    break
  end
end
fclose(fidb);

IJDM = IDM*JDM;
npad = 4096-mod(IJDM,4096);
toto = ones(npad,1);
knew = 41;

fprintf('%s domain, ID=%i JD=%i\n',R,IDM,JDM);

% Read relaxaion file:
fida = fopen(fina,'r');
A = fread(fida,IJDM,'float32','ieee-be');
A = (reshape(A,IDM,JDM))';
fclose(fida);

A(A==0) = nan;
RL = 1./A*(1/3600); % sec^(-1) -> hours

figure(1);
clf;
contour(HH,[0 0],'Color',[0.8 0.8 0.8]);
hold on;
pcolor(RL); shading flat;
%  caxis([0  c2]);
axis('equal');

caxis([0 2500]);

colorbar;
%Is=max(strfind(ftrca,'/'));
spp=sprintf('%s, Min Relax Time %6.2f hrs',R,min(min(RL)));
title(spp,'interpreter','none','Fontsize',12);

txb = 'check_relax_rmu.m';
bottom_text(txb,'pwd',1);
drawnow
