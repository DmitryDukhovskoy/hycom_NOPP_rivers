% Test runs with 0.08 HYCOM-CICEv5
% Plot cice.*.r fields
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

regn = 'ARCc0.08';
expt = 122;
s_fig = 0;

%fld = 'glbrad';
%fld = 'lwdflx';
fld = 'airtmp';

f_plt = 0; % plot instant 1-hr fields
f_avrg = 1; % plot average fields

% Record to start:
iyr = 2017;
imo = 6;
mday = 7;
hr = 12; 
nrec = 1;  % # of recrods to read/average
% 
drc1 = datenum(iyr,imo,mday,hr,0,0);  % start record to read
drc2 = drc1+nrec-1;  % last record
% 1-hr fields
dt = 1; 
Ntot = 24/dt*31; % all records for 31 days

DN = [datenum(iyr,imo,1):dt/24:datenum(iyr,imo,1)+31];
dnn = abs(DN-drc1);
irc1 = find(dnn == min(dnn));
if min(dnn)>dt/24, error('Could not find record'); end;
dnn = abs(DN-drc2);
irc2 = find(dnn == min(dnn));

btx = 'plot_cice_forcing008.m';

%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

ID = nn;
JD = mm;
IJDM=ID*JD;
lrec = ID*JD*8;


%pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/incoming/';
pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/atm_force008/';

finr = sprintf('%scice.%s.r',pthbin,fld);
fid1 = fopen(finr,'r','ieee-be');  %


%k0=nrec;
%stat=fseek(fid1,k0*lrec,-1);
% For CICE forcing fields are not fixed-length records
% as for HYCOM

switch(fld)
 case('lwdflx');
  c1=220;
  c2=360;
 case('glbrad');
  c1=100;
  c2=400;
 case('airtmp');
  c1=-10;
  c2=20;
end

CMP = create_colormap6(400,c1,c2);
cmp = CMP.colormap;

% Daily average:
Asum = HH*0;
cntr=0;

frewind(fid1);
for ii=1:irc2
  dmm=fread(fid1,[ID, JD],'float64');  % 

		if (isempty(dmm));
						fprintf('E-o-F, month = %i\n',ii);
		end
  AA=dmm';
  maxA = max(max(AA(500:2000,50:1575)));
  minA = min(min(AA(500:2000,50:1575)));
  fprintf('Reading  %i (%i), min/max = %8.5g %8.5g\n',ii,nrec,minA,maxA);

  if ii<irc1, continue; end;
  if ii>irc2, break; end;

		Asum = Asum+AA;
		cntr = cntr+1;

% Plot instanteneous hourly field: 
		if f_plt==1; 
				fprintf('Plotting ...\n');
%				AA(HH>0)=nan;

				nDay = floor(ii/24)+1;
				nHr  = ii-(nDay-1)*24-1;
				stl = sprintf('%i/%2.2i/%2.2i, hr=%i, max=%5.4g W/m2',iyr,imo,nDay,nHr,maxA);

				figure(1); clf;
				set(gcf,'Position',[1750 519 777 810]);
				pcolor(AA); shading flat;
    
				colormap(cmp);
				caxis([c1 c2]);

				hold on;
				contour(HH,[0 0],'k');
				axis('equal');
				title(stl);

				set(gca,'xlim',[50 1575],...
												'ylim',[250 2000],...
												'fontsize',12);

				clb = colorbar;
				set(clb,'Position',[0.89 0.11 0.025 0.81],...
												'Fontsize',12);

				bottom_text(btx,'pwd',1);

				keyboard
		
		end

end
%  fprintf('ii=%i, dmm=%ix%i\n',ii,size(dmm));
%keyboard

Asum=Asum/cntr;
%Asum(HH>0)=nan;

if strncmp(fld,'airtmp',6)
  Asum=Asum-273.15;
end

fprintf('Plotting ...\n');
%if f_avrg==1
		figure(10); clf;
		set(gcf,'Position',[1750 519 777 810]);
		pcolor(Asum); shading flat;
		colormap(cmp);
		caxis([c1 c2]);

  nhrs = (cntr-1)*dt;  
  maxAav = max(max(Asum));

  dv1=datevec(drc1);
  dv2=datevec(drc2);
  stl = sprintf('cice.%s.r, Avrg: %i/%2.2i/%2.2i:%2.2i-%2.2i/%2.2i:%2.2i, nrec=%i, max=%5.4g W/m2',...
                fld,dv1(1:4),dv2(2:4),cntr,maxAav);



		hold on;
		contour(HH,[0 0],'k');
		axis('equal');
		title(stl);

		set(gca,'xlim',[50 1575],...
										'ylim',[250 2000],...
										'fontsize',12);

		clb = colorbar;
		set(clb,'Position',[0.89 0.11 0.025 0.81],...
										'Fontsize',12);

		bottom_text(btx,'pwd',1);


%end



