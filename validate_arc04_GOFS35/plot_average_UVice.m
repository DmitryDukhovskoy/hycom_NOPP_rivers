% Arctic Ocean 0.04 HYCOM-CICEv5 GOFS3.5
% Greenland runoff, no passive tracers
% analyze monthly mean sea ice fields
% extract from daily instantenous output
%
% Plot time-averaged velocity vectors
% of sea ice U 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 23;
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat = pthout;

f_conc = 0;
f_uv   = 0;

mmean = 1; % = 0 - instanteneous or monthly mean fields
% Averaging time window:
dnmb1 = datenum(2017,06,01);
dnmb2 = datenum(2017,06,30);
dv1   = datevec(dnmb1);
dv2   = datevec(dnmb2);
yr    = dv1(1);
mo    = dv1(2);
mday  = dv1(3);
YR    = yr;


%pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%i_cice/',expt,yr);
pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_023/data/2017_cice_dltEddTmltEAPJRA/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

btx = 'plot_average_UVice.m';


ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;


fprintf('Average over: %s-%s\n',datestr(dnmb1),datestr(dnmb2));

nrec=0;
for dnmb=dnmb1:dnmb2
  nrec=nrec+1;
  iday=dnmb-datenum(YR,1,1)+1;

		dv0   = datevec(dnmb);
		yr    = dv0(1);
		mo    = dv0(2);
		mday  = dv0(3);
		YR    = yr;

		if mmean 
				flnm = sprintf('%s%3.3i_cice.%i-%2.2i-%2.2i.nc',...
                       pthbin,expt,YR,mo,mday);
		else
				flnm = sprintf('%s%3.3i_cice_inst.%i-%2.2i-%2.2i-00000.nc',...
                       pthbin,expt,YR,mo,mday);
		end


		fprintf('Reading %s \n',flnm);
		ui = squeeze(nc_varget(flnm,'uvel'));
		vi = squeeze(nc_varget(flnm,'vvel'));
  ci = squeeze(nc_varget(flnm,'aice'));

  if ~exist('TLON','var')
				TLON = nc_varget(flnm,'TLON');
				TLAT = nc_varget(flnm,'TLAT');
  end

  if nrec==1
    Ui=ui;
    Vi=vi;
    Ci=ci;
  else
    Ui=Ui+ui;
    Vi=Vi+vi;
    Ci=Ci+ci;
  end
end
Ui=Ui/nrec;
Vi=Vi/nrec;

Si = sqrt(Ui.^2+Vi.^2);


%
% U ice
%
cs1=0;
cs2=0.15;
CMP = create_colormap_WBYR(200,cs1,cs2);
cmu = CMP.colormap;


xl1 = 350;
xl2 = 3160;
yl1 = 1200;
yl2 = 3800;


fn=2;
fprintf('Plotting Ice U,V\n');
sub_plot_uice(Si,Lmsk,xl1,xl2,yl1,yl2,cs1,cs2,fn,cmu);
%contour(Ci,[0.15 0.15],'r','Color',[0.8 0.3 0],'Linewidth',1.5);
% Plot ice u,v vectors
II=find(Ci<1e-3);
Ui(II)=nan;
Vi(II)=nan;
[mm,nn]=size(Ui);


scl=400;
cf=0.3;
beta=20;
lwd=1.0;
v_col=[0 0 0];
dii=40;
for ii=dii:dii:nn
		for jj=dii:dii:mm
				clear u v
				u = Ui(jj,ii);
				v = Vi(jj,ii);
				if isnan(u), continue; end;

    sp = sqrt(u*u+v*v);
    

				x0=ii;
				y0=jj;
     
				x1=x0+u*scl;
				y1=y0+v*scl;
				draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);

		end
end
%stl = sprintf('CPL-TEST HYCOM-CICEv5, Uice, %2.2i/%2.2i/%4.4i',mday,mo,YR);
stl = sprintf('HYCOM-CICEv5, avrg Uice, %2.2i/%2.2i/%4.4i-%2.2i/%2.2i/%4.4i',...
              dv1(3),dv1(2),dv1(1),dv2(3),dv2(2),dv2(1));
title(stl);

bottom_text(btx,'pwd',1,'fontsize',8);










