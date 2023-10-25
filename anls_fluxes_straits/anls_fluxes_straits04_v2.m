% Strait fluxes extracted in 
% extr_TSVdaily_straits04.m
% Combine several years of fluxes
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear

%
% Plot several fluxes
iEX = 9;
EXPT   = sub_cice_experiments;

YR1=2017;
YR2=2020;

%pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.04/data_straits/';
%pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.04/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

btx='anls_fluxes_straits04_v2.m';

Tref1 = -1.8; % Ref T to calc. H flux
Tref2 = 0;    % Ref T to calc. H flux
Sref1 = 34.8;
Sref2 = 34.9;

ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
%fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

% Fram Section is close to Moorings ~79N
%SCT = sub_define_sections04(HH,LON,LAT);
SCT = sub_define_AO_NA_sections04(HH,LON,LAT);
nsct = length(SCT);


CLR=[0 0.4 0.8; ...
     0.8 0.4 0; ...
     1 0.2 0; ...
     0 1 0; ...
     0.8 0 0.6; ...
     1 0.3 0; ...
     0.3 0.7 1; ...
     1 0.7 0.5; ...
     0 0.5 0.9; ...
     0.2 0.7 0.4; ...
     0.1 0.4 0.3; ...
     0.8 0.3 0.9; ...
     1 0.4 0.6; ...
     0.45 0.2 0.85; ...
     0.8 0.7 0.4];


f_map=0;
if f_map==1
  fprintf('Drawing map with segments\n');
  fn=1;
%  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:500:-100],'Color',[0.9 0.9 0.9]);

  for ip=1:14
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%        'Linewidth',2.5,'Color',[1. 0.6 0]);
    clr=CLR(ip,:);
    clr = [0.9,0.2,0];

    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
         'Linewidth',3.5,'Color',clr);
    Ip=SCT(ip).IJPR(1);
    Jp=SCT(ip).IJPR(2);

%    plot(Ip,Jp,'.','Markersize',14,'Color',clr);
  end

  axis('equal');
  set(gca,'xlim',[600 2800],...
          'ylim',[200 4800]);

  bottom_text(btx,'pwd',1);

keyboard
end

YR=YR1;
sv2km = 1e6*3600*24*365*1e-9; % Sv -> km3/yr
msv2km = 1e3*3600*24*365*1e-9; % mSv -> km3/yr


%for ii = 1:length(iEX)
%  ix1 = iEX(ii);
ix1 = iEX;
nmex = EXPT(ix1).cice_opt;
expt = EXPT(ix1).Nmb;
res  = EXPT(ix1).res;
pthmat=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/strait_fluxes/',...
							 expt);
%fmatout=sprintf('%shycom%3.3i_%3.3i_%s_StraitFluxes_%4.4i.mat',...
%									pthmat,res*100,expt,nmex,YR);

FLX = combine_Uflx(pthmat,YR1,YR2,expt,nmex,res);
%FLX = sub_readUflx(fmatout); 


nsct = length(FLX);
ii=1;
%for ik = 1:nsct
for ik = 1:3
	nm = FLX(ik).Name;
	VFlx = FLX(ik).Vol*1e-6;
	N = length(VFlx);
	tm1 = FLX(ik).TM;
% Discard Jan 1 - some suspecious values, restart? 
%    if tm1(1)==datenum(YR,1,1);
%      VFlx(1) = nan;
%    end

	mvol = nanmean(VFlx);
	verr=nanstd(VFlx)/sqrt(N); % mean error

	sinfo = sprintf('%3.2f%s: VFlx=%5.2f+/-%5.2f Sv, %6.1f+/-%6.1f km3/y3',...
	res,nmex, mvol, verr, mvol*sv2km, verr*sv2km);

	SINFO{ii,ik} = sinfo;

	figure(ik); 
	if ii==1
		clf;
%		set(gcf,'Position',[1538         894        1001         424]);
		set(gcf,'Position',[298         779        1057         560]);
		axes('Position',[0.09 0.4 0.84 0.4]);
		hold on;
	end
	d1 = datenum(YR,1,1);
	tmo1 = days2months(tm1);
%    vol1 = FLX(ik).Vol*1.e-6;
	vol1 = VFlx;

%    clr = CLR(ix1,:);
	clr = [0,0.4,0.8];
	plot(tmo1,vol1,'-','Color',clr,'Linewidth',2.5);

  xl2=ceil(max(tmo1));
  dxl=1;
  if xl2>12, dxl=3; end;
	if ii==1
		set(gca,'Tickdir','out',...
					'xlim',[1 xl2],...
					'xtick',[0:dxl:xl2],...
					'Fontsize',14,...
					'xgrid','on',...
					'ygrid','on');
		stl = sprintf('0.04HYCOM-CICE, VolFlx (archv), %i, %s',YR,nm);
		title(stl,'Interpreter','none');
	end

	lgd_txt{ii}=sprintf('%3.2f-%s',res,nmex);

	if ii==length(iEX)
		axes('Position',[0.05,0.1,0.4,0.1]);
		for ipp=1:length(iEX)
			info_txt{ipp} = SINFO{ipp,ik};
		end
		text(0.02,0.5,info_txt,'Fontsize',10);
		set(gca,'ylim',[0.4 0.6],...
						'xlim',[0 0.5],...
						'visible','off');
%      lgd = legend(lgd_txt);
%      set(lgd,'Position',[0.7 0.18 0.2 0.11])
		axes('Position',[0.6,0.1,0.35,0.2]);
		hold on;
		for ipp=1:length(iEX)
			iclr = iEX(ipp);
			yy0=length(iEX)-ipp+1;
			xx0=0.1;
			dxx=0.1;
			clr=CLR(iclr,:);
			plot([xx0, xx0+dxx],[yy0 yy0],'-','Color',clr,'linewidth',2.5);
			text(xx0+1.2*dxx,yy0,lgd_txt{ipp},'Fontsize',10);
		end
		set(gca,'xlim',[xx0, xx0+5*dxx],...
						'ylim',[0.5 length(iEX)+0.5],...
						'visible','off');
		
		bottom_text(btx,'pwd',1);
	end


	drawnow

	fprintf('%s %s\n',nm,sinfo);
	fprintf(' --------------------------------- \n\n');
end
	











