% Strait heat fluxes extracted in 
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
iEX = 9;  % experiment
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
ncf = 'TW';
switch(ncf)
 case('TW')
  cff = 1e-12;  % TW = 1e12 W
 case('PW')
  cff = 1e-15;
end


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

FLX = combine_HTflx(pthmat,YR1,YR2,expt,nmex,res);

nsct = length(FLX);
ii=1;
%for ik = 1:nsct
for ik = 1:3
	nm = FLX(ik).Name;
	Hflx1 = FLX(ik).Hflx1*cff;
  Hflx2 = FLX(ik).Hflx2*cff;
	N = length(Hflx1);
	tm1 = FLX(ik).TM;

  mhflx1 = nanmean(Hflx1);
  verr1  = nanstd(Hflx1)/sqrt(N); % mean error
  mhflx2 = nanmean(Hflx2);
  verr2  = nanstd(Hflx2)/sqrt(N); % mean error

  sinfo = sprintf('%3.2f%s: Hflx1=%5.2f+/-%5.2f%s, Tref1=%5.2f^oC, Hflx2=%5.2f+/-%5.2f%s, Tref2=%5.2f^oC',...
	        res, nmex, mhflx1, verr1, ncf, Tref1, mhflx2, verr2, ncf, Tref2);

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
%    clr = CLR(ix1,:);
	clr1 = [0,0.4,0.8];
  clr2 = [1,0.6,0];
	plot(tmo1,Hflx1,'-','Color',clr1,'Linewidth',2.5);
  hold on
  plot(tmo1,Hflx2,'-','Color',clr2,'Linewidth',2.5);


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
		stl = sprintf('0.04HYCOM-CICE, HeatFlx (%s), %i, %s',ncf,YR,nm);
		title(stl,'Interpreter','none');
	end

  axes('Position',[0.7,0.12,0.25,0.1])
  hold on
  plot([0 0.5],[1 1],'-','Color',clr1,'Linewidth',2.5);
  text(0.7,1,sprintf('HFlux1 Tref=%3.2fC',Tref1))
  plot([0 0.5],[0.5 .5],'-','Color',clr2,'Linewidth',2.5);
  text(0.7,0.5,sprintf('HFlux2 Tref=%3.2fC',Tref2))
  set(gca,'ylim',[0. 1.2],...
            'xlim',[0 2],...
            'visible','off');
  


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
%		axes('Position',[0.6,0.1,0.35,0.2]);
%		hold on;
%		for ipp=1:length(iEX)
%			iclr = iEX(ipp);
%			yy0=length(iEX)-ipp+1;
%			xx0=0.1;
%			dxx=0.1;
%			clr=CLR(iclr,:);
%			plot([xx0, xx0+dxx],[yy0 yy0],'-','Color',clr,'linewidth',2.5);
%			text(xx0+1.2*dxx,yy0,lgd_txt{ipp},'Fontsize',10);
%		end
%		set(gca,'xlim',[xx0, xx0+5*dxx],...
%						'ylim',[0.5 length(iEX)+0.5],...
%						'visible','off');
		
		bottom_text(btx,'pwd',1);
	end


	drawnow

	fprintf('%s %s\n',nm,sinfo);
	fprintf(' --------------------------------- \n\n');
end
	











