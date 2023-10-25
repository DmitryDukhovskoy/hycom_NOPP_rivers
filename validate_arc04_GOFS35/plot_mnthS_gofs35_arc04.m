% Monthly mean T/S
% Plot only 1 layer at a time
% if need to layer-average - modify the code
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2017;
YR2 = 2020;
pfld  = 'salin';
plr = 1;    % layer to extract
dday = 5;   % day stepping
f_extr = 0; % =1 extract and average fields, =0 - load in saved fields

% Choose experiment:
ixx    = 9; % experiment name and dir - check with EXPT - expt 023
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;



regn = 'ARCc0.04';
rg = 9806;
hgg=1e20;


s_fig = 0;

%mplt = [12,1,2]; % average over N months, if not, mplt=1month
mplt = [6,7,8]; % average over N months, if not, mplt=1month
nav = length(mplt);
im1=mplt(1);
im2=mplt(2);
dnmb1 = datenum(YR1,im1,1);
iday1 = dnmb1-datenum(YR1,1,1)+1;
mday = datenum(YR1,im2,15);
endD = ndays_month(mday);
dnmb2 = datenum(YR1,im2,endD);
iday2 = dnmb2-datenum(YR1,1,1)+1;

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/atl_water/',expt);


ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

SM={'J','F','M','A','M','J','J','A','S','O','N','D'};


fprintf('S fields, layer=%i, %i-%i days: %i-%i\n',plr,YR1,YR2,iday1,iday2);

fmat = sprintf('%shycom004_%3.3i_%s_Lr%2.2i_%2.2i-%2.2i_%i-%i.mat',...
               pthmat,expt,pfld,plr,im1,im2,YR1,YR2);
if f_extr==1
  icc = 0;  % overall counter
  for yr=YR1:YR2
    for iday=iday1:dday:iday2
      dnmb = datenum(yr,1,1)+iday-1;
			DV=datevec(dnmb);

			pthbin=sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%4.4i_%s/',...
											expt,yr,texpt);

			fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
			finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);

			if yr==2017 & DV(2)<4
				pthbin='/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_BL99Tfrz/';
				fina = sprintf('%s022_archv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
				finb = sprintf('%s022_archv.%4.4i_%3.3i_00.b',pthbin,yr,iday);
			end
					
			if ~exist(fina,'file') | ~exist(finb,'file')
				fprintf('Not found %s or %s, skipping ...\n',fina,finb);
				continue;
			end

			tic;
			[F,n,m,nlr] = read_hycom(fina,finb,pfld,'r_layer',plr);
			F(F>hgg)=nan;
			SS=squeeze(F);
		 
			icc=icc+1;
			if icc==1
				AA = SS;
			else
				AA = AA+SS;
			end 

    end  % month
  end  % year

  SS = AA/icc;
  fprintf('Total records = %i\n',icc);
  fprintf('Min/max Value = %6.2f / %6.2f\n',nanmin(nanmin(SS)),nanmax(nanmax(SS)));
  fprintf('Saving averaged fields to %s\n',fmat);
  save(fmat,'SS')
else
  fprintf('Loading %s\n',fmat);
  load(fmat);
end

S = SS;

% Subp. N. Atl
xlim1 = 740;
xlim2 = 2500;
ylim1 = 300;
ylim2 = 2200;

% Greenland
%xlim1 = 900;
%xlim2 = 2200;
%ylim1 = 600;
%ylim2 = 2200;
%xlim1 = 
% SE Greenland:
%xlim1 = 1200;
%xlim2 = 1720;
%ylim1 = 700;
%ylim2 = 1280;

% Zoom in S Greenland
fzoom = 0;
if fzoom==1
	xlim1 = 1045;
	xlim2 = 1580;
	ylim1 = 805;
	ylim2 = 1460;
end

fprintf('Plotting ...\n');
nf = 1;
sb=33;
stl = sprintf('ARCc0.04-%3.3i, Mean S (cntr0=%3.1f) %i-%i, mo:%i-%i',expt,sb,YR1,YR2,im1,im2);
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
pfld='salin';
hps = [0.93 0.1 0.02 0.8];
Fpos = [1350  255  1034  1085]; % Figure gcf position
sub_plot_scalar_v0(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
		'c1',30,'c2',35,'cmp',2,'clbpos',hps,'figpos',Fpos);
contour(S,[30:0.5:36],'k');
contour(S,[sb sb],'k','Linewidth',1.8);

txtb = 'plot_mnthS_gofs35_arc04.m';
bottom_text(txtb,'pwd',1);



