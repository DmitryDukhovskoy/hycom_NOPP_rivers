% Experiments with CICE5
% Select experiments to plot (or all of them)
% Plot Extracted time series of sea ice area & volume,
% from test experiments
% extracted in extr_ice_area_vol.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

pthout   = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat   = pthout;

%pEXPT = [2,8,9];
pEXPT = 8;

YR1=2017;
YR2=2020;



CLR = [0.8 0.6 0; ...
       0.5 1   0.2; ...
       0   0.6 0.9; ...
       0.7 0   0.9; ...
       1  0.1  0.1; ...
       0  0.7  0.2; ...
       0.7 0.  0.2; ...
       0  0.9  0.7;
       0  0    1];
lwd = 3;
%
btx = 'plot_ice_area_vol_yrs.m';

% Plot NSIDC Ice extent
f_nsidc = 1;
if f_nsidc>0
% Monthly mean from NSIDC NOAA - extent (including holes inside the contour)
% Note that calculated CICE area is ice area (excluding holes)
% Get Satellite-based ice area for total arctic for specified time period
	d1 = datenum(2017,1,1);
	d2 = datenum(2020,12,31);

	[ArObs,TMobs] = sub_get_obs_icearea(d1,d2,'TotArc');
	%[ArObs,TMobs] = sub_get_obs_icearea(d1,d2,'ArcOc');

end

% ---------------------------------------

fmatout = sprintf('%sseaice_VolArea_tests%i.mat',pthmat,2018);
load(fmatout);
nexpt = length(ICE);
for ixx=1:nexpt
	nm = ICE(ixx).Name;
	i0=find(pEXPT==ixx);
	if isempty(i0);
		fprintf('%i %s\n',ixx,nm);
	else
		fprintf('%i %s <---- Ploting \n',ixx,nm);
	end
end


for ifg=2:3   
	figure(ifg); clf;
	set(gcf,'Position',[925         775        1109         556]);
	axes('Position',[0.09 0.5 0.85 0.4]);
	hold on;

	vmx = 0;
	vmn =1e10;
	amx = 0;
	amn = 1.e10;

	for ixx=1:nexpt
		i0=find(pEXPT==ixx);
		if isempty(i0); continue; end;
		TM=[];
		Vi=[];
		Ai=[];
    Ei=[];   % ice extent
		for YR=YR1:YR2
			fmatout = sprintf('%sseaice_VolArea_tests%i.mat',pthmat,YR);
			fprintf('Loading %s\n',fmatout);
			load(fmatout);

			TM = [TM;ICE(ixx).TM'];
			Vi = [Vi;ICE(ixx).Vol_km3'];
			Ai = [Ai;ICE(ixx).Area_km2'];
      Ei = [Ei;ICE(ixx).IExt_km2'];
		end

		DV  = datevec(TM);
		YR  = DV(1,1);
		dJ1 = datenum(YR1,1,1);
		jD  = TM-dJ1+1;

	%  clr = CLR(ixx,:);
		clr = [0., 0.3, 0.8];
    clrE= [0,0.7, 1];

		if ifg==2
			plot(jD,Vi,'-','Color',clr,'Linewidth',lwd);
		else
			plot(jD,Ai,'-','Color',clr,'Linewidth',lwd);
%      plot(jD,Ei,'-','Color',clrE,'Linewidth',lwd); does not work
		end

		vmx = max([vmx,max(Vi)]);
		vmn = min([vmn,min(Vi)]);
		amx = max([amx,max(Ai)]);
		amn = min([amn,min(Ai)]);

	end


  if ifg==3 & f_nsidc==1
		for ii=1:length(TMobs)
				d1 = TMobs(ii,1)-dJ1+1;
				d2 = TMobs(ii,2)-dJ1+1;
				ai = ArObs(ii);

				plot([d1 d2],[ai ai],'k-','Linewidth',2);
		end
  end


	if ifg==2
		yy2=vmx;
    dyy = 5000;
		title('Ice Volume, km^3');
	else
		yy2=amx;
    dyy = 1e6;
		title('Ice Area, km^2');
	end

	icc=0;
	imo1=1;
	imo2=12;
	for YR=YR1:YR2
		for imm=imo1:imo2
			dnmb = datenum(YR,imm,1);
			jd = dnmb-dJ1+1;
			if imm==1
				plot([jd jd],[0 1.2*yy2],'--','Color',[0.8 0.8 0.8]);
			else
				plot([jd jd],[0 1.2*yy2],':','Color',[0.9 0.9 0.9]);
			end
			icc=icc+1;
			xtck(icc)=jd;
		end
	end



	set(gca,'tickdir','out',...
					'xtick',xtck,...
					'xlim',[xtck(1) xtck(end)],...
					'ylim',[0 1.1*yy2],...
					'ytick',[0:dyy:1.2*yy2],...
					'ygrid','on',...
					'fontsize',12);


	% Legend
	axes('Position',[0.1 0.1 0.6 0.3]);
	hold on;
	yy0=length(ICE)+1;
	xx1=0.05;
	xx2=xx1+0.2;
	xx3=xx2+0.1;

	for ixx=1:length(ICE);
		i0=find(pEXPT==ixx);
		if isempty(i0); continue; end;

		nm = ICE(ixx).Name;
	%  clr = CLR(ixx,:);

		iy = yy0-ixx+1;
		plot([xx1 xx2],[iy iy],'-','Color',clr,'Linewidth',lwd);
		text(xx3,iy,nm,'Fontsize',12,'Interpreter','none');
	end

%  if ifg==3
%    iy = iy-1;
%    plot([xx1 xx2],[iy iy],'-','Color',clrE,'Linewidth',lwd);
%    text(xx3,iy,'HYCOM pseudo IceExt','Fontsize',12);
%  end
     
  if ifg==3 & f_nsidc==1
	  iy = iy-1;
	  plot([xx1 xx2],[iy iy],'-','Color',[0 0 0],'Linewidth',lwd);
	  text(xx3,iy,'NOAA IceExt','Fontsize',12);
  end

	set(gca,'xlim',[0 5],...
					'ylim',[1 yy0+1],...
					'visible','off')

	bottom_text(btx,'pwd',1);

end




