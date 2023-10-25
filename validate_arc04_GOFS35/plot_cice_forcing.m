% Plot cice.*.r fields
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

regn = 'ARCc0.04';
expt = 22;
s_fig = 0;

f_plt = 0; % plot instant 1-hr fields
f_avrg = 1; % plot average fields
%plt_fld = 'meltb'; 
%plt_fld = 'frzmlt'; 
%plt_fld = 'fswdn';

btx = 'plot_cice_forcing.m';

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;

ID = nn;
JD = mm;
IJDM=ID*JD;
lrec = ID*JD*8;


pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/incoming/';

fld = 'glbrad';
finr = sprintf('%scice.%s.r',pthbin,fld);
fid1 = fopen(finr,'r','ieee-be');  %

% 1-hr fields
nrec = 24*31; % all records for 31 days
t1 = 1;
t2 = nrec; 
dt = 1; 

iyr = 2017;
imo = 7;


frewind(fid1);
%k0=nrec;
%stat=fseek(fid1,k0*lrec,-1);
% For CICE forcing fields are not fixed-length records
% as for HYCOM

c1=0;
c2=500;
CMP = create_colormap6(400,c1,c2);
cmp = CMP.colormap;

% Daily average:
Asum = HH*0;
cntr=0;

for ii=1:t2
  dmm=fread(fid1,[ID, JD],'float64');  % 
		if (isempty(dmm));
						fprintf('E-o-F, month = %i\n',ii);
		end
  AA=dmm';
  maxA = max(max(AA(1000:4000,100:3150)));
  minA = min(min(AA(1000:4000,100:3150)));
  fprintf('Reading  %i (%i), min/max = %8.5g %8.5g\n',ii,nrec,minA,maxA);

  if ii>=t1 & ii<=t2 & mod(ii-1,dt)==0
    Asum = Asum+AA;
    cntr = cntr+1;
 
% Plot instanteneous hourly field: 
    if f_plt==1; 
						fprintf('Plotting ...\n');
						AA(HH>0)=nan;

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

						set(gca,'xlim',[100 3150],...
														'ylim',[500 4000],...
														'fontsize',12);

						clb = colorbar;
						set(clb,'Position',[0.89 0.11 0.025 0.81],...
														'Fontsize',12);

						bottom_text(btx,'pwd',1);

						keyboard
				
    end

  end

end
%  fprintf('ii=%i, dmm=%ix%i\n',ii,size(dmm));
%keyboard

Asum=Asum/cntr;
Asum(HH>0)=nan;

if f_avrg==1
		figure(10); clf;
		set(gcf,'Position',[1750 519 777 810]);
		pcolor(Asum); shading flat;
		colormap(cmp);
		caxis([c1 c2]);

  nDay=floor(t1/24)+1;
  nhrs = (cntr-1)*dt;  
  maxAav = max(max(Asum));

  stl = sprintf('Average: %i/%2.2i/%2.2i, nhrs=%i, max=%5.4g W/m2',iyr,imo,nDay,nhrs,maxAav);



		hold on;
		contour(HH,[0 0],'k');
		axis('equal');
		title(stl);

		set(gca,'xlim',[100 3150],...
										'ylim',[500 4000],...
										'fontsize',12);

		clb = colorbar;
		set(clb,'Position',[0.89 0.11 0.025 0.81],...
										'Fontsize',12);

		bottom_text(btx,'pwd',1);


end



