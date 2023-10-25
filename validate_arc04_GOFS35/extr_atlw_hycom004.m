% Extract T and Layer thicknes of Atl. Water
% in specified regions (Canada Basin, BG, ...)
% for 0.04 and 0.08 HYCOM-CICE5 and save fields for python
% 
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps/
startup

format long g
clear all
close


YR1 = 2020;
YR2 = 2020;
dday = 30;   % 1 a month

s_mat=1;

ixx    = 9; % experiment name and dir - check with EXPT - expt 023
%ixx    = 6;  % expt 022 original 
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;


rg = 9806;
hgg=1e20;

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/atl_water/',expt);

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
btx = 'extr_atlw_hycom004.m';

fprintf('arc04-%3.3i Atlantic Water %i-%i save mat=%i', expt,YR1,YR2,s_mat);


% Get topo:
ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY;
[X,Y]=meshgrid([1:nn],[1:mm]);

Adpth=zeros(mm,nn)*nan;
nxx=0;
%mm=2520;
%nn=1600;
%ll=33;

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

icc=0;
for yr=YR1:YR2
  cc = 0;

  ATL = struct;
  ATL.Iocn = Iocn;
  ATL.Zmax  = [];
  ATL.LrThck= [];
  ATL.Tmax  = [];
  ATL.Smax  = [];
  ATL.HtCnt = [];

  fmat=sprintf('%shycom004_%3.3i_%s_AtlLr_%4.4i.mat',...
                    pthmat,expt,texpt,yr);
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
    mo=imo;
    md=15;
    dnmb=datenum(yr,mo,md);
    DV=datevec(dnmb);
    iday=dnmb-datenum(yr,1,1)+1;

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
      fprintf('Not found %s or %s\n\n',fina,finb);
      if dday==1; continue; end;
%
% Search for close days:
      for ik=1:dday
        iday = iday-1;
        dnmb = datenum(dJ1)+iday-1;
        DV   = datevec(dnmb);
        fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
        finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
        if exist(fina,'file');
          fprintf(' Found closest file: %s\n',fina);
          break;
        end;
      end
% If went back to iday0 - skip
      if ik==dday & ~exist(fina,'file');
        continue;
      end
    end

    cc=cc+1;

    tic;
    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
    F(F>hgg)=nan;
    T=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'salin');
    F(F>hgg)=nan;
    S=F;

    [ZM,ZZ] = sub_zz_zm(fina,finb,HH);

    Tmax = Iocn*nan; % Tmax
    Smax = Iocn*nan; % S at Tmax
    Zmax = Iocn*nan; % depth of T max
    dHL  = Iocn*nan; % thickness Atl. layer
    HAtl = Iocn*nan;  % heat content in Atl L.
		for ii=1:nI % points inside domain
			i0 = Iocn(ii);

			tt  = squeeze(T(:,i0));
			tt0 = squeeze(T(:,i0));
			z0  = ZZ(:,i0);
			iz0 = max(find(z0>-100));
			iz2 = max(find(z0>-800));
			if ~isempty(iz0)
				tt(1:iz0) = -999;
				tt0(1:iz0) = -999;
			end
      if ~isempty(iz2)
				tt(iz2:end) = -999;
      end
	
			tm = max(tt);  % max T within the layer:
			Tmax(ii) = tm;

			iz = find(tt==max(tt),1);
			Zmax(ii) = ZM(iz,i0);

			iz = min(find(tt>=0));
      if isempty(iz), continue; end;
      zm = ZM(:,i0);
      zz = ZZ(:,i0);
	%keyboard
	% Depth max T
	% Interpolate to find exact z
	% upper 0 interface
			iz1=iz-1;
			iz2=iz;
			t1 = tt(iz1);
			t2 = tt(iz2);
			z2 = ZZ(iz2,i0);
			z1 = ZZ(iz1,i0);
			t0 = 0;
			dtdz = (t2-t1)/(z2-z1);
			dz0 = (t0-t1)/(dtdz);
			zt0 = z1+dz0;
%			Zt0(ii) = zt0;
			
	% Thickness of warm Atl. Layer
	% Interpolate to get exact thickness
			Ipz = find(tt0>0);
			iz1=Ipz(1)-1;
			iz2=Ipz(1);
			z1=ZM(iz1,i0);
			z2=ZM(iz2,i0);
			t2=tt0(iz2);
			t1=tt0(iz1);
			t0=0;
			dtdz = (t2-t1)/(z2-z1);
			dz0 = (t0-t1)/(dtdz);
			zt0 = z1+dz0;
			zup = zt0;
			dzup = abs(dz0);
			
			iz1=Ipz(end);
			iz2=iz1+1;
      if iz2<=nlr
				z1=ZM(iz1,i0);
				z2=ZM(iz2,i0);
				t2=tt0(iz2);
				t1=tt0(iz1);
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
%			tt0 = squeeze(T(:,i0));
      tt  = squeeze(T(:,i0));
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
				plot(T(:,i0),ZM(:,i0),'.-');
				hold on;
				plot([0 0],[min(ZM(:,i0)) 0],'r--');
				plot(tt0(Ipz),ZM(Ipz,i0),'go');
				plot(t0,zt0,'m*'); % located 0
				plot(t0,zbt,'m*'); % located 0
			end
	
      if mod(ii,100000)==0
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
    ATL(imo).TM    = dnmb;
	
		
		if s_mat>0
			fprintf('saving %s\n',fmat);
			save(fmat,'ATL');
		end
		
		fprintf('Processed 1 record %6.4f min\n\n',toc/60);
		fprintf('-----------------------------------\n')
			
	%keyboard
  end    % months
end




    








