% Arctic Ocean 0.04 HYCOM-CICEv5 GOFS3.5
% Greenland runoff, no passive tracers
%
% Plot time-averaged atmospheric winds velocity vectors
% from the atm forcing fields interpolated onto HYCOM domain
% wndewd, wndnwd_YYY[a-l].[ab] 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 22;

f_extr = 0;
% Averaging time window:
hYR = 117;
AB  = 'g';
dltT = 6; % hr, skip step for time averaging


% 1hr data output
% direct access data
dS  = 1; % day to start
hrS = 0; % hr to start
dt  = 1; % 1hr data
dE  = 31;
hrE = 23;

recS = (dS-1)*24/dt+hrS+1; % record to start, hr, day=dS
recE = (dE-1)*24/dt+hrE+1; % End record to read


btx = 'plot_average_UVwind.m';


pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/atm_force04/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthout   = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';

fmatUV = sprintf('%savrg_UVwind_%3.3i%s.mat',pthout,hYR,AB);


ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

ID = nn;
JD = mm;
IJDM=ID*JD;
lrec = ID*JD*8;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

if f_extr==1
		finUa = sprintf('%swndewd_%3.3i%s.a',pthbin,hYR,AB);
		finUb = sprintf('%swndewd_%3.3i%s.b',pthbin,hYR,AB);
		fidU1 = fopen(finUa,'r');  %

		finVa = sprintf('%swndnwd_%3.3i%s.a',pthbin,hYR,AB);
		finVb = sprintf('%swndnwd_%3.3i%s.b',pthbin,hYR,AB);
		fidV1 = fopen(finVa,'r');  %


		usum = zeros(JD,ID);
		vsum = zeros(JD,ID);
		nrec = 0;
		dltR = dltT/dt;   % recrods to skip

		for iRec = recS:dltR:recE
				tic;
		% Skip to start record:
				statu=fseek(fidU1,(iRec-1)*(IJDM+npad)*4,-1);
				statv=fseek(fidV1,(iRec-1)*(IJDM+npad)*4,-1);

				fprintf('Reading record %i, last rec %i\n',iRec,recE);
				dmm=fread(fidU1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
				dm1=fread(fidU1,npad,'float32','ieee-be');  % read npad 
		%  fprintf('ii=%i, dmm=%ix%i\n',ii,size(dmm));
		%keyboard
				dmm=reshape(dmm,ID,JD);
				A=dmm';
				usum = usum+A;

				dmm=fread(fidV1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
				dm1=fread(fidV1,npad,'float32','ieee-be');  % read npad 
		%  fprintf('ii=%i, dmm=%ix%i\n',ii,size(dmm));
		%keyboard
				dmm=reshape(dmm,ID,JD);
				A=dmm';
				vsum = vsum+A;

				nrec = nrec+1;

				fprintf('Nrec=%i, Read 1 rec, %8.6f sec\n',nrec,toc);
		end

		Uw = usum/nrec;
		Vw = vsum/nrec;

		fprintf('Saving mean atm UV %s\n',fmatUV);
		save(fmatUV,'Uw','Vw');

else
  fprintf('Loading %s\n',fmatUV);
  load(fmatUV);
end

Uw(HH>=0)=nan;
Vw(HH>=0)=nan;
Sw = sqrt(Uw.^2+Vw.^2);



cs1=0;
cs2=6;
CMP = create_colormap_WBYR(200,cs1,cs2);
cmu = CMP.colormap;


xl1 = 350;
xl2 = 3160;
yl1 = 1200;
yl2 = 3800;


fn=1;
fprintf('Plotting Atm U,V\n');
sub_plot_uice(Sw,Lmsk,xl1,xl2,yl1,yl2,cs1,cs2,fn,cmu);
%contour(Ci,[0.15 0.15],'r','Color',[0.8 0.3 0],'Linewidth',1.5);
% Plot ice u,v vectors
[mm,nn]=size(Ua);


scl=12;
cf=0.3;
beta=20;
lwd=1.0;
v_col=[0 0 0];
dii=40;
for ii=dii:dii:nn
  for jj=dii:dii:mm
    clear u v
    u = Uw(jj,ii);
    v = Vw(jj,ii);
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
  finUa = sprintf('%swndewd_%3.3i%s.a',pthbin,hYR,AB);

stl = sprintf('avrg Uatm(m/s) from wndewd/wndnwd_%3.3i%s.[ab], days: %i-%i, dt=%ihrs',...
              hYR,AB,dS,dE,dt);

title(stl,'interpreter','none');

bottom_text(btx,'pwd',1,'fontsize',8);






