% Plot area-mean vorticity 
% calculated in calc_meanVort_arc04
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


YR1=2017;
YR2=2020;
dday=5;

RR = 200; % ring size (diamter), km 
Ilr = 5;  % Z layer where vort is calcualted

%regn = 'natl'; % Natl region
%regn = 'arctA'; % 
%regn = 'arctB'; % 
regn = 'ArctOc'; % 


% This code is for archv output fields
% CHeck if archm - do not add U/V barotropic !!! - use extr_*mean_ 
s_mat=1;  % =2 - load saved and start from the last record

ixx    = 9; % experiment name and dir - check with EXPT - expt 023
%ixx    = 6;  % expt 022 original 
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;


%Zmn = 100; % average over the top Zmn m
rg  = 9806;
hgg = 1e20; % 


pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/strait_fluxes/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

% Start from:
ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

[mm,nn]=size(LON);


% Grid of the mean vorticity field - check with calc_meanVort_arc04.m
% Subset grid points:
switch(regn)
 case('ArctOc')
  IJr = [ 720        3734
         760        3167
         882        2933
        1060        2750
        1161        2541
        1273        2425
        1415        2333
        1491        2274
        1649        2124
        1806        2041
        1892        1957
        2167        1949
        2314        2233
        2405        2575
        2446        3067
        2217        3751
        1390        3868];

end

xl1 = min(IJr(:,1));
xl2 = max(IJr(:,1));
yl1 = min(IJr(:,2));
yl2 = max(IJr(:,2));

dii = 6;
HHp=HH(yl1:dii:yl2,xl1:dii:xl2);
[IIp,JJp] = meshgrid([xl1:dii:xl2],[yl1:dii:yl2]);

icc = 0;
for iyr = YR1:YR2
  fmat = sprintf('%s%3.3i_meanVort_%s_%4.4i.mat',...
           pthmat,expt,regn,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);

  TM = MVRT.TM;
  DV = datevec(TM);
  IW = find(DV(:,2)>=10 | DV(:,2)<4);
  IS = find(DV(:,2)>=4 & DV(:,2)<10);

  Ioc = MVRT.Ioc;
  if icc==0
    AW = Ioc*0;
    AS = Ioc*0;
  end

  icc = icc+1;
  dmm = MVRT.VRT(IW,:);
  AW = AW+mean(dmm,1)';  
  dmm = MVRT.VRT(IS,:);
  AS = AS+mean(dmm,1)';

end
AW = AW/icc;
AS = AS/icc;


VRTW = IIp*0;
VRTS = IIp*0;
[ms,ns]=size(IIp);
AA = HH*0;
AA(HH>=0)=nan;
AA(Ioc) = AW;

BB = HH*0;
BB(HH>=0)=nan;
BB(Ioc) = AS;
for isb=1:ns
  for jsb=1:ms
    isb0=IIp(jsb,isb);
    jsb0=JJp(jsb,isb);
    VRTW(jsb,isb)=AA(jsb0,isb0);
    VRTS(jsb,isb)=BB(jsb0,isb0);
  end
end


pgrd=45;
Vrt0=VRTW;
Hmsk=Vrt0*0;
Hmsk(abs(Vrt0)>0)=1;
Vrt = sub_fltr(Vrt0,pgrd,Hmsk);
VRTW = Vrt;

Vrt0=VRTS;
Vrt = sub_fltr(Vrt0,pgrd,Hmsk);
VRTS = Vrt;


clr1 = colormap_red(200);
clr2 = flipud(colormap_blue(200));
for ik=1:10
  clr1(ik,:)=[1 1 1];
  clr2(end-ik+1,:)=[1,1,1];
end
cmp = [clr2;clr1];
cmp = smooth_colormap(cmp,9);
cmp = smooth_colormap(cmp,9);

hmsk=HHp;
hmsk(HHp<0)=nan;

c1 = -5.e-7;
c2 =  5.e-7;

for ik = 1:2
  if ik==1
    AA = VRTW;
    ssn = 'Oct-Mar';
  else
    AA = VRTS;
    ssn = 'Apr-Spt';
  end

	figure(ik); clf;
	set(gcf,'position',[1000         483         752         833]);

	ax1=axes('Position',[0.1 0.2 0.8 0.75]);
  pcolor(hmsk); shading flat;
  colormap([0 0 0]);
  freezeColors;

  hold on

	pcolor(AA); shading flat;

  caxis([c1 c2]);
  colormap(cmp);
  axis('equal');
%	set(gca,'xlim',[xlim1 xlim2],'ylim',[ylim1 ylim2]);
	set(gca,'xtick',[],'ytick',[]);
	clr=[0.9 0.9 0.9];
%	plot_gridlines(45,10,1,clr,LON,LAT);
  stl = sprintf('Area-mean vorticity, 023 0.04 HYCOM-CICE5 %i-%i %s',YR1,YR2,ssn);
	title(stl,'Fontsize',12,'Interpreter','none');
  hb = colorbar('horizontal');
  set(hb,'Position',[0.2 0.1 0.6 0.015],'Fontsize',12,...
          'Ticks',[-5e-7:1e-7:5e-7]);
  set(gca,'Color',[0.8 0.8 0.8]);
  
%  hp = get(ax1,'Position');
%  axes('Position',hp);
%  set(gca,'Visible','off');


end





