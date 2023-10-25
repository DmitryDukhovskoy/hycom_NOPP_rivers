% Analyze GFWA fluxes in the SPNA
% majors straits and SPNA contour
% Analyze flow structure (2D) and 1D 
%
% data extracted in extr_Trcr_daily_straits08.m
% Tracer is converted into GFWA
% using fract of Tr in a cell given overall Tracer mass in the domain

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear

get_data=2;  % =1 read data and save; =2 load saved

expt=112;
TV=11;
% Years to average:
YR1=2016;
YR2=2016;

%  Extract data 
% Average years should be within the time intervals 
% extracted 
YRe1=2006;
YRe2=2019;

if YR1<YRe1 | YR2>YRe2
  error('Averged time is outside saved mat years %i %i',YRe1,YRe2);
end


dday=7;
nTr=1;  % extract 1 tracer at a time, 1 - Greenland


hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation
ptharch= '/nexsan/people/ddmitry/';
pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthriv = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_mat/';
pthtopo= '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='anls_TrFlux2D_SPNA.m';



ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


% Fram Section is close to Moorings ~79N
SCT = sub_define_SPNA_sections(HH,LON,LAT);
nsct = length(SCT);


CLR=[0 0.4 0.8; ...
     0.8 0.4 0; ...
     1 0.2 0; ...
     0 1 0; ...
     0.8 0 0.6; ...
     0.8 1 0; ...
     0.3 0.7 1; ...
     1 0.7 0.5; ...
     0.9 0.3 0.6; ...
     0.2 0.7 0.4; ...
     0.1 0.4 0.3; ...
     0.8 0.3 0.9; ...
     1 0.4 0.6];

f_map=0;
if f_map==1
  fprintf('Drawing map with segments\n');
  fn=1;
%  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:500:-100],'Color',[0.9 0.9 0.9]);

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%        'Linewidth',2.5,'Color',[1. 0.6 0]);
   clr=CLR(ip,:);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
         'Linewidth',2.5,'Color',clr);
    Ip=SCT(ip).IJPR(1);
    Jp=SCT(ip).IJPR(2);

    plot(Ip,Jp,'.','Markersize',14,'Color',clr);
  end

  axis('equal');
  set(gca,'xlim',[300 1300],...
          'ylim',[100 1300]);

  bottom_text(btx,'pwd',1);

%keyboard
end


ftrspna=sprintf('%sanls_trc_str%i-%i.mat',pthmat,YRe1,YRe2);

if get_data == 1
		frv = sprintf('%sGreenland_cumFWFlux.mat',pthriv);

		f_griv = 0; %=1 - rederive cumulative Greenland FWFlux from Bamber
		if f_griv==1
				sub_get_GrRunoff(friv)
		end
		fprintf('f_griv %i, Loading %s\n',f_griv,frv);
		load(frv);


		% Get tracer fraction in grid cells
		IntSrf = 1; % = 1 - integrates over all layers from 0 - ilv, for ilv=1 is the same
		%ilv = 1; % 0-50m
		%ilv = 2; % 50-150m
		%ilv = 3; % 150-300m
		%ilv = 4; % 300-500 m
		ilv = 5; % whole depth <- do not use this

		% Vertical layers
		LRS = load('LRS.dat');
		nlrs= length(LRS);


		zz1 = LRS(ilv,1);
		zz2 = LRS(ilv,2);

		if IntSrf==1; zz1=0; end;

		dz=abs(zz2-zz1);
		hZ = abs(LRS(ilv,2)-LRS(ilv,1));

		fprintf('FW Volume Subpolar Gyre, ilv=%i, %i - %i, nTr=%i, 1993-2016\n',ilv,zz1,zz2,nTr);


		cc=0;
		for YR=YRe1:YRe2
				fmatout=sprintf('%shycom008_%3.3i_Trcr%2.2i_StraitFluxesDay_%4.4i.mat',...
																						pthmat,expt,nTr,YR);
				fprintf('Loading %s \n',fmatout);
				load(fmatout);

				cc=cc+1;
				TM=SCT(1).Time;
				DV=datevec(TM);
				Is=find(DV(:,2)>4 & DV(:,2)<10);
				Iw=find(DV(:,2)<=4 | DV(:,2)>=10);

		%

				imo0=-1;
				MTr_dom=[];
				for it=1:length(TM)
						DV = datevec(TM(it));
						iyr = DV(1);
						imo = DV(2);
						if imo~=imo0
								expt0=110;
        if YR<1999
								  pthmat2=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt0);
        else
								pthmat2=sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/%3.3i/data_matTr/',expt);
        end
								fmat2 = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat2,nTr,iyr,imo);
								fprintf('Loading %s\n',fmat2);
								load(fmat2);
								imo0=imo;
						end
		% find whole-depth layer
						nlr = length(TRCR);
						ibtm=5; % whole-depth tracer mass
						Tr_dom=squeeze(TRCR(ibtm).MassTr_kg); % tracer integrated ove whole water depth
						amm=Tr_dom./(abs(HH).*DX.*DY);  % kg/m3
						Tr_dom(amm<=0.125)=nan;
						MTr_dom(it,1) = nansum(nansum(Tr_dom)); % overall mass of the tracer in the domain
				end

		%
		% Get tracer fluxes across straits
				nsct = length(SCT);
				for isc=1:nsct
						nm=SCT(isc).Name;
						Nrm=SCT(isc).Norm;
						Xn=Nrm(:,1);
						Yn=Nrm(:,2);

		% Obtain time series of overall Mass tracer for all time in this year:
						TM = SCT(isc).Time;
						DVm = datevec(TM);
						if DVm(1,1)~=YR
							fprintf('Section %s date error %i %i !!!!!!!! \n\n',isc,YR,DVm(1,1));
						end
%
% GFWA fluxes across straits:
						TFlx=SCT(isc).TrFlx;  % kg/s of tracer
		%
		% Convert kg/s of tracer into m3/s of GFWA:
		% Greenland runoff anomaly for given date
						ism = find(Ygr==DV(1,1));
						if isempty(ism), ism=length(cFWF); end; % after 2016 - same 
						fwf0 = cFWF(ism); % km3

						MTr_dom(MTr_dom==0)=nan;
						rr=TFlx./MTr_dom;
						vflx=rr*fwf0;      % km3/sec
						nrc=length(vflx);
						FLX(isc).GrFlx_mean(cc,1)=sum(vflx)/nrc*3600*24*365; % km3/yr

		%
		% Compute mean vol flux, Sv
						ZZ=SCT(isc).ZZintrp;
						nlr=length(ZZ);
						dZ=abs(diff(ZZ));
						dZ(nlr)=dZ(nlr-1);
						dl=SCT(isc).segm_dL;
						[DL,DZ]=meshgrid(dl,dZ);

						vflx = SCT(isc).VolFlx_m3s;
						VFlx = mean(vflx,1);
						FLX(isc).VolFlux_Sv(cc,1) = VFlx*1e-6; % Sv 

% Mean 2D maps of GFWA concentration
		%     
						dmm=SCT(isc).Unrm;
						Un=squeeze(nanmean(dmm,1)); % normal flow
						[a1,a2]=size(Un);
		% Project on norm: one of the components is 0
		%    inrm=sign(Xn+Yn);
		%    inrm=inrm(:);
		%    for ikk=1:a1
		%      Un(ikk,:)=Un(ikk,:).*inrm';
		%    end
						
						dmm=SCT(isc).Trcr;
						Tr=squeeze(nanmean(dmm,1)); % tracer kg/m3 - tracer concentration in grid cell

%						nlr=a1;
%						ZZ=SCT(isc).ZZintrp;
%						dZ=abs(diff(ZZ));
%						dZ(nlr)=dZ(nlr-1);
%						dl=SCT(isc).segm_dL;
%						[DL,DZ]=meshgrid(dl,dZ);
						
		%
		% Convert kg of tracer into m3 of GFWA in 1 m3 of ocean:
		% Greenland runoff anomaly for given date
						ism = find(Ygr==DV(1,1));
						if isempty(ism), ism=length(cFWF); end; % after 2016 - same 
						fwf0 = cFWF(ism); % km3

						MTr_dom_mn=mean(MTr_dom); 
						rr=Tr./MTr_dom_mn;
						GrVol=rr*fwf0*1e9;      % m3 of GFWA in 1 m3 of ocean, 1e-3 m3 = 1 L
						
						FLX(isc).Name = nm;
						FLX(isc).YR(cc,1)=YR;
						FLX(isc).GrVol_m3m3(cc,:,:) = GrVol;
						FLX(isc).Umn_ms(cc,:,:) = Un;

% Compute Tr Flux depth-intgr:
% annual mean
						TrFlx=nansum(Tr.*Un.*DL.*DZ);
      GrFlx=TrFlx/MTr_dom_mn*fwf0*3600*24*365; % km3/yr per 1 grid pnt
      FLX(isc).GrFlx_km3yr(cc,:)=GrFlx;
%
				end

		end

  fprintf('Saving %s\n',ftrspna);
  save(ftrspna,'FLX');
else
  fmatout=sprintf('%shycom008_%3.3i_Trcr%2.2i_StraitFluxesDay_%4.4i.mat',...
                      pthmat,expt,nTr,YR2);
  fprintf('Loading %s \n',fmatout);
  load(fmatout);

  fprintf('Loading %s\n',ftrspna);
  load(ftrspna);
end

%c1=0;
%c2=0.15;
%nint=300;
%CMP=create_colormap8(nint,c1,c2);
%cmp=CMP.colormap;
%cnt=CMP.intervals;

%
% Select years to average
YRS = FLX(1).YR;
iy1  = find(YRS==YR1);
iy2  = find(YRS==YR2);


% Plotting
fprintf('Plotting GFWA tracer-based sections\n');
iSCT=[1,2,3];
nplt=length(iSCT);

FLX(1).xlm=[];
FLX(1).ylm=[];
FLX(2).xlm=[];
FLX(2).ylm=[];
FLX(3).xlm=[0 350];
FLX(3).ylm=[-400 0];

for ii=1:nplt
  isc=iSCT(ii);
  VTr=FLX(isc).GrVol_m3m3(iy1:iy2,:,:);
  mVTr=squeeze(mean(VTr,1));
  [a1,a2]=size(mVTr);

  xx=SCT(isc).long(1:a2);
  yy=SCT(isc).latd(1:a2);

  DX=distance_spheric_coord(yy,xx,yy(1),xx(1))*1e-3;
  ZZ=SCT(isc).ZZintrp;
  Js=SCT(isc).J(1:a2);
  Is=SCT(isc).I(1:a2);
  Hb=[];
  for ipp=1:length(Is);
    Hb(ipp,1)=HH(Js(ipp),Is(ipp));
  end;
  DX=DX(:);
  xBtm=[DX(1); DX; DX(end)];
  yBtm=[-6000; Hb; -6000];

% fill nans for plotting:
  for ipp=1:length(Is)
    dmm=mVTr(:,ipp);
    i1=max(find(~isnan(dmm)));
    if ~isempty(i1)
      dmm(i1+1:end)=dmm(i1);
    end
    mVTr(:,ipp)=dmm;
  end

% Normal V - project on section plane 
% to get rid of zig-zags
  UU=FLX(isc).Umn_ms;
  mUU=squeeze(mean(UU,1));
  xn=SCT(isc).Norm(:,1);
  yn=SCT(isc).Norm(:,2);
  Ix0=find(xn==0);
  Ix1=find(abs(xn)==1);
  if length(Ix1) > length(Ix0)  % projecting on section along Y-axis with Xnorm=1
    Jzz=Ix0;
  else  % projecting on X-axis section
    Jzz=Ix1;
  end

  mUU(:,Jzz)=-999;
  II=find(mUU(1,:)~=-999);
  mUU=mUU(:,II);
  xU=DX(II);
  
% Plot log
  lmVTr = log(mVTr);

  zmin=min(Hb);
  yl1=1.02*zmin;
  xl1=0;
  xl2=max(DX);

  if ~isempty(FLX(isc).xlm), xl1=FLX(isc).xlm(1); xl2=FLX(isc).xlm(2); end;
  if ~isempty(FLX(isc).ylm), yl1=FLX(isc).ylm(1); end;

  c1=-7;
  c2=-4.5;
  CTR=[]; 
  sub_plot_TrSct(isc,DX,ZZ,lmVTr,xBtm,yBtm,xU,mUU,c1,c2,yl1,xl1,xl2,CTR);

  nm=SCT(isc).Name;
  sttl=sprintf('%s, ln(GFWA) m3 in 1 m3, log, %i-%i',nm,YR1,YR2);
  title(sttl);

  bottom_text(btx,'pwd',1);
%keyboard

end  


% Report
fprintf('GFWA fluxes through SPNA straits\n');
fprintf('  Years  %i - %i\n',YR1,YR2)
fprintf('Positive flux is north and east\n');
% Check overall transp in SPNA
ICH = [-1; -2; -6; -7; -8; -9; 10; 11; 12];
Vspna_dg=0;
for isc=1:nsct
  nm=SCT(isc).Name;
  flx=FLX(isc).GrFlx_mean(iy1:iy2);
  mflx=mean(flx);
  nn=length(flx);
  serr=std(flx)/sqrt(nn);

  vflx=FLX(isc).VolFlux_Sv(iy1:iy2);
  mVflx=mean(vflx);
  sVerr=std(vflx)/sqrt(nn);

		jpp = find(abs(ICH)==isc);
		if ~isempty(jpp)
				fsgn=sign(ICH(jpp));
				Vspna_dg=Vspna_dg+fsgn*mVflx;  % Sv
		end

% Compute only inflows at Davis and Denamrk only near coast:
% to avoid returning tracers with other currents
  DL=SCT(isc).segm_dL;
		xdst=cumsum(DL);
  flx=0;
  TRin=0;
		if isc==1
				jj1=1;
%				jj2=max(find(xdst<=250000));
    jj2=length(xdst);
				GF = FLX(isc).GrFlx_km3yr(iy1:iy2,jj1:jj2);
%    GF(GF>0)=nan;
    TRin=nansum(GF,2);
  end
  if isc==2
    jj1=1;
    jj2=max(find(xdst<=210000));
    GF = FLX(isc).GrFlx_km3yr(iy1:iy2,jj1:jj2);
    GF(GF>0)=nan;
    TRin=nansum(GF,2);
  end
  flx=TRin;
  mflx=mean(flx);
  nn=length(flx);
  serr=std(flx)/sqrt(nn);




  fprintf('%2i %s: VolFlx  %6.4g +/- %6.4g Sv\n',isc,nm,mVflx,sVerr);
  if isc==1 | isc==2
    fprintf('GFWA    %6.4g +/- %6.4g km3/yr\n', mflx,serr);
  end

end

fprintf(' SPNA overall transport = %6.4g Sv\n',Vspna_dg);

keyboard

% Plot 1D diagrams of depth-integrated GFWA flux km3/yr
%for isc=1:nsct
IPLT=[1; 2];
for ill=1:length(IPLT);
  isc=IPLT(ill);
  TrFlx = nanmean(FLX(ill).GrFlx_km3yr);
  dmm=TrFlx;
		nav=5;
		for jj=1:2
				for ik=1:length(TrFlx);
						i1=max([ik-nav,1]);
						i2=min([ik+nav,length(TrFlx)]);
						dmm(ik)=nanmean(TrFlx(i1:i2));
				end;
				dmm(isnan(TrFlx))=nan;
				TrFlx=dmm;
		end

  [a1,a2]=size(TrFlx);
  xx=SCT(isc).long(1:a2);
  yy=SCT(isc).latd(1:a2);

  DX=distance_spheric_coord(yy,xx,yy(1),xx(1))*1e-3;

  dd=mean(diff(DX)); 
  nm=SCT(isc).Name;

  figure(12+ill); 
  axes('Position',[0.09 0.6 0.85 0.32]);
  hold on;
  plot(DX,TrFlx/dd,'Linewidth',2.5);  % km3 yr per 1 km of section
  set(gca,'tickdir','out',...
          'xlim',[0 max(DX)],...
          'xgrid','on',...
          'ygrid','on',...
          'Fontsize',14);

  stl = sprintf('%s Depth-integrated GFWA flux, km3/yr per 1 km of section');
  title(stl);
end    
  









