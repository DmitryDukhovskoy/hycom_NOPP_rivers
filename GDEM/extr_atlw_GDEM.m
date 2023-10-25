% Extract T and layer thickness of Atl. Water 
% from GDEM climatology
% GDEM is interpolated onto HYCOM 0.08 grid in interp_GDEM2HYCOM08.m

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 1; % 

pthmat = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/023/atl_water/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthin = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/data_mat/';  % GDEM fields


% Extract data for:
phi0 = 50;  % southmost latitude
zmin = -200; % depth isobath - within which perform the analysis

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo);
HH  = nc_varget(ftopo,'Bathymetry'); 
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH); 
[X,Y]=meshgrid([1:nn],[1:mm]);


% Consider deep Arctic ocean:
Iocn = find(LAT>65 & HH<-500.);
nI = length(Iocn);

f_bth=0;
if f_bth==1
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:1000:-100],'b');
  plot(X(IN),Y(IN),'r.');
  axis('equal');
end

cc = 0;

ATL = struct;
ATL.Iocn = Iocn;
ATL.Zmax  = [];
ATL.LrThck= [];
ATL.Tmax  = [];
ATL.Smax  = [];
ATL.HtCnt = [];

fmat=sprintf('%sGDEM_AtlLr.mat',pthmat);
im1 = 1;
if s_mat==2
	fprintf('Loading %s\n',fmat);
	if ~exist(fmat,'file')
		fprintf('File does not exist\n');
	else
		load(fmat);
		im1=length(ATL)+1;
	end
	fprintf('Continue from month = %i\n',im1);
end

for imo=im1:12
  ftin = sprintf('%sTave_GDEM_%2.2i.mat',pthin,imo);
  fsin = sprintf('%sSave_GDEM_%2.2i.mat',pthin,imo);

  if ~exist(ftin,'file')
    fprintf('Input fields for %i missing, skipping ...\n',imo);
    continue;
  end

  fprintf('Processing month = %i\n',imo);

  tic;
  fprintf('Loading %s\n',ftin);
  TAV=load(ftin);
  fprintf('Loading %s\n',fsin);
  SAV=load(fsin);
  
  T = TAV.Tav;
  S = SAV.Sav;
%
% Bug in interpolation code - smaller T,S domain
  [nlr,mt,nt] = size(T);
  if mt~=mm | nt~=nn
    T(:,mt+1:mm,:) = nan;
    S(:,mt+1:mm,:) = nan;
    T(:,:,nt+1:nn) = nan;
    S(:,:,nt+1:nn) = nan;
  end 
      
  ZZ = -1*SAV.zz;
  dz = abs(diff(ZZ));
  nlr = length(ZZ);

	Tmax = Iocn*nan; % Tmax
	Smax = Iocn*nan; % S at Tmax
	Zmax = Iocn*nan; % depth of T max
	dHL  = Iocn*nan; % thickness Atl. layer
	HAtl = Iocn*nan;  % heat content in Atl L.
	for ii=1:nI % points inside domain
		i0 = Iocn(ii);

    [jx,ix] = ind2sub([mm,nn],i0);

		tt  = squeeze(T(:,i0));
		tt0 = squeeze(T(:,i0));
		z0  = ZZ;
    zz  = ZZ;
%
% Delete below bottom T:
    hb = HH(i0);
    izb = max(find(zz>=hb));
    tt(izb+1:end)=nan;
    tt0(izb+1:end)=nan;

		iz = min(find(tt>0));
		if isempty(iz), continue; end;

% Assuming AtlLayer in the Arctic Ocean is between ~200 and 700 m
		iz0 = max(find(z0>-100));
%    iz0 = min(find(tt<=0));
		iz2 = max(find(z0>-800));
    if isempty(iz0); iz0=1; end; 

		tt(1:iz0) = -1.e-6;
		if ~isempty(iz2)
			tt(iz2:end) = -2.222;
		end

		tm = nanmax(tt);  % max T within the layer:
		Tmax(ii) = tm;

		iz = find(tt==max(tt) & tt>0.,1);
    if isempty(iz); continue; end;
		Zmax(ii) = ZZ(iz);

%		zm = ZM(:,i0);
%keyboard
% Thickness of warm Atl. Layer
% Interpolate to get exact thickness
    ifrz = min(find(tt<=0));
    if isempty(ifrz) | ifrz>iz  % no negative t above max T
      zup = 0;  % should not happen because surf T inforced to be <0
    else
			Ipz = find(tt>0);
			iz1=Ipz(1)-1;
			iz2=Ipz(1);
			z1=ZZ(iz1);
			z2=ZZ(iz2);
			t2=tt(iz2);
			t1=tt(iz1);
			t0=0;
			dtdz = (t2-t1)/(z2-z1);
			dz0 = (t0-t1)/(dtdz);
			zt0 = z1+dz0;
			zup = zt0;
			dzup = abs(dz0);
    end

		iz1=Ipz(end);
		iz2=iz1+1;
		if iz2<=nlr & iz2<izb
			z1=ZZ(iz1);
			z2=ZZ(iz2);
			t2=tt(iz2);
			t1=tt(iz1);
			dtdz = (t2-t1)/(z2-z1);
			dz0 = (t0-t1)/(dtdz);
			zbt = z1+dz0;
			dzbt= abs(dz0);
		else
			zbt = HH(i0);
		end

		dh = abs(zbt-zup); % Atl.L. thickness

		dHL(ii)=dh;

% Heat content within the Atl. Layer:
%		tt  = squeeze(T(:,i0));
		iz1 = Ipz(1);
		iz2 = Ipz(end);
		dzz = abs(diff(zz));

		zIndx=[];
		dZ = [];
		ii1 = max(find(zz>=zup));
		dZ(1) = abs(zup-zz(ii1+1));
		ii2 = max(find(zz>=round(zbt)));
		if ii2-ii1<2, ii2=ii1+2; end;
		if ii2>nlr, continue; end;
		dZbt = abs(zz(ii2)-zbt);
		zIndx=(ii1+1:ii2-1)';
		dZ=[dZ;dzz(zIndx);dZbt];

		ss = S(:,i0);
		Cp = 4200; % J/kg K
		Tref1 = -1.8; % Ref T to calc. H flux
		rhow = sw_dens0(ss,tt);
		hc = Cp*rhow.*(tt-Tref1);
		itop = zIndx(1)-1;
		ibot = zIndx(end)+1;
		ibot = min([ibot,nlr]);
		dmm = [hc(itop);hc(zIndx);hc(ibot)];
		HAtl(ii) = dmm'*dZ;  % J/m2 - heat content integrated in Atl Layer
		Smax(ii) = ss(iz);

% Plot T;
		chckT=0;
		if chckT==1
			figure(10); clf;
			plot(tt0,ZZ,'.-');
			hold on;
			plot([0 0],[min(ZZ) 0],'r--');
      plot([-2 1],[hb hb],'k-');
      set(gca,'ylim',[hb 0]);
			plot(tt0(Ipz),ZZ(Ipz),'go');
			plot(t0,zt0,'m*'); % located 0
			plot(t0,zbt,'m*'); % located 0
		end

		if mod(ii,20000)==0
			fprintf(' Atl.Layer, %4.1f%% done ...\n',ii/length(Iocn)*100);
			fprintf('     Tmax=%4.1f Smax=%4.1f Zmax=%6.1f dHatl=%6.1f HtCont=%6.4d GJ/m2\n\n',...
								Tmax(ii), Smax(ii),Zmax(ii),dHL(ii),HAtl(ii)*1e-9);
		end
	end  % points

	ATL(imo).Zmax  = Zmax;
	ATL(imo).LrThck= dHL;
	ATL(imo).Tmax  = Tmax;
	ATL(imo).Smax  = Smax;
	ATL(imo).HtCnt = HAtl;

	if s_mat>0
		fprintf('saving %s\n',fmat);
		save(fmat,'ATL');
	end
	
	fprintf('Processed 1 record %6.4f min\n\n',toc/60);
	fprintf('-----------------------------------\n')
		
%keyboard
end    % months




