% Plot sea ice volume daily 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 22;
% Test simulations with corrected CICE-HYCOM coupling
% Sept 2016, note different lenght of the runs
%expt = 0; % restart 0.04 GLBc GOFS3.5 created for 0.04 ARCc
%expt = 221; % fixed CICE-HYCOM coupled, restart from old run 0.04ARCc-CICEv5 Sept 1, 2016
%expt = 222; % fixed CICE-HYCOM coupled, CICE restart 0.04GLBc 01/09/2017
%expt = 223; % fixed CICE-HYCOM coupled, CICE&HYCOM restarts 0.04GLBc 01/09/2017 - CICE error T cond 9 days
%expt = 224; % dixed CICE-HYCOM coupled, CICE&HYCOM restarts 0.04GLBc 01/09/2017 - reduced dt=60 sec


pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat = pthout;
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

btx = 'seaice_vol_daily.m';


ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

c1=0;
c2=4;
pos=[0.12 0.1 0.82 0.82];
xl1= 50;
xl2= nn;
yl1= 500;
yl2= 4000;
%CMP = colormap_PBYR(200,c1,c2);
CMP = create_colormapBGY(200,c1,c2);
cmp = CMP.colormap;

% HYCOM-CICE 0.04 old run with corrected coupling (ARCc0.04 GOFS3.5) using
% CICE restart created from GLBb GOFS2.5 reanalysis Sept 1 2017
%pth1 = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_rest_ciceGOFS3.5_0901/';
%pth2 = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_rest_allGOFS3.5_0901/';


YR1=2016;
YR2=2017;
yr=YR1;
cc=0;
TM=[datenum(YR1,9,1):datenum(YR2,8,31)];
ntm=length(TM);
for itm=1:ntm
  dnmb=TM(itm);
  dv=datevec(dnmb);
  yr=dv(1);
  YR=yr;
  imo=dv(2);
  mday=dv(3);

	switch(expt)
		case(22);
			pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i_cice/',yr);
			fin = sprintf('%s022_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',pthbin,yr,imo,mday);

   if yr>2016  % mean fields
     fin = sprintf('%s022_cice.%4.4i-%2.2i-%2.2i.nc',pthbin,yr,imo,mday);
   end

		case(0); % restart
			pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/restart_022/';
			fina = sprintf('%srestart_116i.a',pthbin);
			finb = sprintf('%srestart_116i.b',pthbin);
		case(221)
% Only last day of the run is saved
			pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data_coupledCICE_restartOLDrun/';
			if mday==1
					pthbin='/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2016_cice/';
			end
			fin = sprintf('%s022_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',pthbin,yr,imo,mday);
		case(222)
			pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_rest_ciceGOFS3.5_0901/';
			fin = sprintf('%s022_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',pthbin,yr,imo,mday);
		case(223)
			pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_rest_allGOFS3.5_0901/';
			fin = sprintf('%s022_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',pthbin,yr,imo,mday);
		case(224)
			pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_dtCICE_rest_allGOFS_0901/'; 
			fin = sprintf('%s022_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',pthbin,yr,imo,mday);
	
		end

		cc=cc+1;

% Missing Dec. 2016 - get from surface fields
  if yr==2016 & imo==12
    iday=dnmb-datenum(yr,1,1)+1;
    pthsrf=sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%4.4i_surf/',yr);
    fina=sprintf('%s022_arche.%4.4i_%3.3i_12.a',pthsrf,yr,iday);
    finb=sprintf('%s022_arche.%4.4i_%3.3i_12.b',pthsrf,yr,iday);
     
    fld='sih';
    [F,n,m,l] = read_hycom_arche(fina,finb,fld);
    F(F>1e20)=nan;
    F(F<1e-6)=0;
    Hi=squeeze(F);

    fld='sic';
    [F,n,m,l] = read_hycom_arche(fina,finb,fld);
    F(F>1e20)=nan;
    F(F<1e-6)=0;
    Ci=squeeze(F);

    dmm = nansum(nansum(Hi.*Ci.*Acell*1e-9));  % km3
    Vol(cc,1)=dmm;
    HiL=Hi;
    yrL=yr;
    mdayL=mday;
    continue 

  end

		if ~exist(fin,'file')
				Vol(cc,1)=nan;
				fprintf('Not found %s\n',fin);
				continue
		end

		fprintf('Reading %s \n',fin);
		Hi = squeeze(nc_varget(fin,'hi'));
%		Si = squeeze(nc_varget(fin,'hs'));
		Ci = squeeze(nc_varget(fin,'aice'));
%		Ui = squeeze(nc_varget(fin,'uvel'));
%		Vi = squeeze(nc_varget(fin,'vvel'));
		dmm = nansum(nansum(Hi.*Ci.*Acell*1e-9));  % km3
		Vol(cc,1)=dmm;
		HiL=Hi;
		yrL=yr;
		mdayL=mday;

		if (cc==1), Hi1=Hi; end;
		if (cc==2), Hi2=Hi; end;
end

%P=T

xlbl=[];
mold=0;
kk=0;
DV=datevec(TM);
for ik=1:length(TM);
  dnmb=TM(ik);
  YR=DV(ik,1);
  mdays=ndays_month(dnmb);
  imo=DV(ik,2);
  mday=DV(ik,3);
%  xtime(ik)=imo+(mday-1)/mdays;

  Jd1=datenum(YR,1,1);
  Lday=datenum(YR,12,31)-Jd1;
  iday=dnmb-Jd1;
  xtime(ik)=YR+iday/Lday;

  if imo~=mold
    kk=kk+1;
    xtck(kk)=xtime(ik);
    mold=imo;
    xlbl{kk}=sprintf('%2.2i',imo);
  end
end

fprintf('Plotting results ...\n');

figure(1); clf;
axes('Position',[0.1 0.5 0.84 0.42]);
plot(xtime,Vol,'linewidth',2);
set(gca,'tickdir','out',...
        'xlim',[xtime(1)-0.001 xtime(end)+0.001],...
        'xtick',xtck,...
        'xticklabel',xlbl,...
        'ylim',[0 1.05*max(Vol)],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
xlabel('Months');

stl=sprintf('%3.3i, %i, sea ice vol, km^3',expt,yr);
title(stl);

btx='seaice_vol_daily.m';
bottom_text(btx,'pwd',1,'position',[0.05 0.36 0.4 0.05]);


% 
%
f_hice=0;
if f_hice==1
		figure(2); clf;
		pcolor(Hi1); shading flat;
		caxis([0 4]);
		colorbar
		title('GLBc CICE restart: 022_cice_inst.2016-09-01-00000.nc','Interpreter','none');
		bottom_text(btx,'pwd',1);

		figure(4); clf;
		pcolor(HiL); shading flat;
		caxis([0 4]);
		colorbar
		stl=sprintf('%3.3i 004AO HYCOM-CICE, Hice, %4.4i/%2.2i/%2.2i',expt,yr,imo,mdayL);
  title(stl,'Interpreter','none');
		bottom_text(btx,'pwd',1);

% Plot Hice difference:
  dH = HiL-Hi1;
  cmp1=colormap_red(200);
  cmp1(1,:)=[1 1 1];
  cmp1(2,:)=[1 1 1];
  cmp2=flipud(colormap_blue(200));
  cmp2(end,:)=[1 1 1];
  cmp2(end-1,:)=[1 1 1];
  cmp=[cmp2;cmp1];
  cmp=smooth_colormap(cmp,5);
  
  figure(5); clf;
  pcolor(dH); shading flat;
  colormap(cmp); 
  caxis([-0.8 0.8]);
  colorbar

  set(gca,'Color',[0.7 0.7 0.7]);

  stl=sprintf('%3.3i 004AO HYCOM-CICE, \\DeltaHice, %2.2i/%2.2i-09/01',expt,imo,mdayL);
  title(stl);
%
% Change in 1 day
  dH = Hi2-Hi1;
  figure(6); clf;
  pcolor(dH); shading flat;
  colormap(cmp);
  caxis([-0.8 0.8]);
  colorbar

  set(gca,'Color',[0.7 0.7 0.7]);

  stl=sprintf('%3.3i 004AO HYCOM-CICE, \\DeltaHice, 09/02-09/01',expt);
  title(stl);
  

end



