% Average over specified time 
% Plot monthly-mean, depth-integrated 
% tracer concentrations 
% The Data are extracted in anls_BG/extr_trcr_mnth.m
% Option: plot anomalies relative to reference Year/month
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_fig = 0;
cc=0; % last saved frame - if need to restart from frame # 
        % = 0 - start from beginning

% Years and months to average
YR1=2016;
YR2=2016;
%IMO=[5,6,7,8,9,10];  % list of months for averaging
IMO=[1,2,3,4,5,6,7,8,9,10,11,12];  % list of months for averaging
im1=IMO(1);
im2=IMO(end);


fprintf('Mean Tr Conc m3 in 1 m3 of sea water, averaged for years: %i-%i, mo: %i-%i\n',...
  YR1,YR2,im1,im2);
%dref = datenum(2006,1,15); % reference year/month, = 0 - no reference 
ilev = 1; %=1 - mixed layer, 0-50; =2 - 50-150m

regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/frames_trcrGreen/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthriv = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_mat/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


hmsk=HH;
hmsk(HH<0)=nan;

nint = 320;
c1 = -4;
c2 = 4;
%CMP = create_colormap2_3(nint,c1,c2);
CMP = colormap_sclr2(nint,c1,c2);
cmp = CMP.colormap;
cmp(1,:) = [1 1 1];
cmp(2,:) = [1 1 1];
cmp(3,:) = [1 1 1];
cmp(4,:) = [1 1 1];
cmp(5,:) = [1 1 1];
cmp(6,:) = [1 1 1];
cmp(7,:) = [1 1 1];
cmp(8,:) = [1 1 1];
cmp(9,:) = [1 1 1];
cmp(10,:) = [1 1 1];
cmp = smooth_colormap(cmp,18);
cmp = smooth_colormap(cmp,18);
cmp(1,:) = [1 1 1];
cnt = CMP.intervals;

ck = 0;
for iyr = YR1:YR2
  for nmm=1:length(IMO)
    imo=IMO(nmm);
    ck=ck+1;
    YRPLT(ck,1)=iyr;
    YRPLT(ck,2)=imo;
  end
end
nrc=ck;

frv = sprintf('%sGreenland_cumFWFlux.mat',pthriv);

f_griv = 0; %=1 - rederive cumulative Greenland FWFlux from Bamber
if f_griv==1
		sub_get_GrRunoff(friv)
end
fprintf('f_griv %i, Loading %s\n',f_griv,frv);
load(frv);

sumGr=HH*0;
%cc=0;
for ic=1:nrc
  iyr = YRPLT(ic,1);
  imo = YRPLT(ic,2);
  
%  pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_mnth%i/',...
%		   regn,expt,iyr);
%  if ~exist(pthfig,'dir')
%    scm = sprintf('mkdir -pv %s',pthfig);
%    system(scm);
%  end
		cc=cc+1;
		fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
		if ~exist(fmat,'file')
				fprintf('File is missing %s\n',fmat);
				continue
		end
		
		fprintf('cc=%i, Loading %s\n',cc,fmat);
		load(fmat);

	
		nTr = 1; % Greenland
		ilv=1;
%    for ilv=1:1
		fprintf('Reading: %i/%2.2i, Tracer %i, Lev %i\n',...
			iyr,imo,nTr,ilv);
									Tr = squeeze(TRCR(ilv).TR(nTr,:,:));
		Tr(Tr<=0)=nan;

%   % Convert kg/m3 of tracer into m3/ m3 of GFWA:
%   1e-3 m3 GFWA in m3 of sea water = 1 L of freshwater
%
		pthmat2=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
		fmat2 = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat2,nTr,iyr,imo);
		fprintf('Loading %s\n',fmat2);
		load(fmat2);
% find whole-depth layer
		nlr = length(TRCR);
		ibtm=5; % whole-depth tracer mass
		Tr_dom=squeeze(TRCR(ibtm).MassTr_kg); % tracer integrated ove whole water depth
		amm=Tr_dom./(abs(HH).*DX.*DY);  % kg/m3
		Tr_dom(amm<=0.125)=nan;
		MTr_dom = nansum(nansum(Tr_dom)); % overall mass of the tracer in the domain

		ism = find(Ygr==iyr);
		if isempty(ism), ism=length(cFWF); end; % after 2016 - same 
		fwf0 = cFWF(ism); % km3

		MTr_dom(MTr_dom==0)=nan;
		rr=Tr./MTr_dom;
		vGr=rr*fwf0*1e9;      % m3 per 1 m3 of sea water
%keyboard
		sumGr=sumGr+vGr;
end
vGr=sumGr/cc;

lTr=log(vGr);
lTr(isnan(lTr))=-1000;
lTr(HH>=0)=nan;

close all
%    ff=figure('visible','off');
%    figure(1); clf;
POS(1,:) = [0.08 0.08 0.85 0.85];
POS(2,:) = [0.06 0.52 0.45 0.45];
POS(3,:) = [0.48 0.52 0.45 0.45];
POS(4,:) = [0.06 0.03 0.45 0.45];
POS(5,:) = [0.48 0.03 0.45 0.45];

if ilv==1
%	  pst = [0.03 0.08 0.4 0.85];
		lvl='0-50m';
else
%	  pst = [0.47 0.08 0.4 0.85];
		lvl='50-150m';
end
%keyboard
c1=-10;
c2=-4;
nint=300;
CMP=create_colormap8(nint,c1,c2);
cmp=CMP.colormap;
cnt=CMP.intervals;


% Plot Mean Tr concentration inside Gr Shelf
GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section
IG = GC.cntr_Iindx;
JG = GC.cntr_Jindx;
[IDX,JDX]=meshgrid([1:nn],[1:mm]);

IN=inpolygon(IDX,JDX,IG,JG);
lTr(IN==0)=nan;


figure(1); clf;
set(gcf,'Position',[1358 451 832 891]);
pst = POS(nTr,:);
axes('Position',pst);
hold on;
% Plot Greenl
Lmsk=HH*0;
Lmsk(HH<0)=1;
pcolor(Lmsk); shading flat;
lcmp=[0.8 0.8 0.8; 1 1 1];
colormap(lcmp);
freezeColors;

pcolor(lTr); shading flat;
%						contour(HH,[0 0],'k','Linewidth',1);
contour(HH,[-5000:1000:-10],'Color',[0.7 0.7 0.7],'Linewidth',0.8);
caxis([c1 c2]);
colormap(cmp);
axis('equal');
% For Greenland:
xlim1 = 450;
xlim2 = 1020;
ylim1 = 350;
ylim2 = 1100;
set(gca,'xlim',[xlim1 xlim2],...
	'ylim',[ylim1 ylim2]);
set(gca,'xtick',[],'ytick',[]);
stl = sprintf('ARCc0.08-%3.3i, Mean YR: %i - %i, Mo: %i-%i GFWA - m3/m3 sea water, %i, %s',...
							expt,YR1,YR2,im1,im2,nTr,lvl);
title(stl,'Fontsize',8);

clb=colorbar;
set(clb,'TickLength',0.02,...
						'Position',[0.88, 0.1, 0.02, 0.8],...
						'Fontsize',11);



%keyboard	  
%				end
btx = 	'avrg_trcrGreen.m';
bottom_text(btx,'pwd',1);


%end   % time loop




