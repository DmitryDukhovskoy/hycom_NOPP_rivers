% The updated Greenland runoff data set 
% From J. Bamber 2018
% Estimate Gr FW anomaly
% for different sectors
% Interested in the GFW anomaly flux to Subpolar Gyre
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

pthG = '/nexsan/people/ddmitry/Net_ocean/arctic_AOregimes/data/GreenlandRunoffv3/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat/';
btx = 'GreenlandRunoff_anomalySPG.m';


s_mat = 0;

fnm = sprintf('%sFWF17.v3.nc',pthG);
tmm = double(nc_varget(fnm,'TIME'));
Xgr = nc_varget(fnm,'lon');
Ygr = nc_varget(fnm,'lat');
Rt  = double(nc_varget(fnm,'runoff_tundra')); % tundra runoff, km3/mo
Rg  = double(nc_varget(fnm,'runoff_ice')); % GrIS runoff - meltwater
D   = nc_varget(fnm,'solid_ice'); % solid ice
LGr = nc_varget(fnm,'LSMGr'); % Greenland mask
TM  = datenum(1958,1,1)+tmm; 
DV  = datevec(TM);
nrc = length(TM);

[mm,nn]=size(LGr);
[II,JJ]=meshgrid((1:nn),(1:mm));

% Define 4 sectors: SPG (Labr + Irm from Denmark to Davis)
% Baffin, AO, Greenland+Iceland 
IJ(1).Name='LabrIrm';
IJ(1).IJ=[251, 476; ...
	  255, 690; ...
	  476, 690; ...
	  499, 469; ...
	  433, 417; ...	  
	  329, 479];
IJ(2).Name='BaffinBay';
IJ(2).IJ=[251, 476; ...
	  214, 171; ...
	  289, 146; ...
	  325, 172; ...
	  329, 479];
IJ(3).Name='GreenlIcel';
IJ(3).IJ=[499, 469; ...
	  517, 465; ...
	  547, 391; ...
	  529, 271; ...
	  473, 149; ...
	  456, 147; ...
	  410, 215; ...
	  433, 417];
IJ(4).Name='AO';
IJ(4).IJ=[289, 146; ...
	  325, 172; ...
	  410, 215; ...
	  456, 147; ...
	  479, 117; ...
	  385, 84];

iR=length(IJ);
for ik=1:iR
  Iv=IJ(ik).IJ(:,1);
  Jv=IJ(ik).IJ(:,2);
  dmm=inpolygon(II,JJ,Iv,Jv);
  IN=find(dmm==1);
  IJ(ik).IN=IN;
end


CLR=[1 0.5 0;...
     0 0.6 1; ...
     0.5 1 0;...
     1 0.2 0.8];
cmsk=[1 1 1; 0.8 0.8 0.8];
f_map=0;
if f_map==1
  figure(10); clf;
  pcolor(LGr); shading flat;
  colormap(cmsk);
  hold on
  for ik=1:iR
    clr=CLR(ik,:);
    IN=IJ(ik).IN;
    plot(II(IN),JJ(IN),'.','Color',clr);
  end
  contour(LGr,[0.9 0.9],'k');
  contour(Xgr,[-180:10:150],'Color',[0.5 0.5 0.5]);
  contour(Xgr,[-30 -30],'Color',[1 0.5 0.5]);
  contour(Ygr,[40:10:89],'Color',[0.5 0.5 0.5]);
  contour(Ygr,[66:2:68],'Color',[1 0.5 0.5]);
  caxis([0 1]);
  axis('equal');
  set(gca,'ydir','reverse',...
	  'xlim',[180 600],...
	  'ylim',[100 700]);
   
  title('Regions Greenland FW Fluxes');
  bottom_text(btx,'pwd',1);
  
end


% Total FW flux = D+Rg+Rt (no CAA, Svalbard here);
% Reproduce Bamber's Fig.3:
clear FWF iD iRt iRg
for it=1:nrc
  dv = datevec(TM(it));
  fprintf('Greenland Runoff, %i/%2.2i/%2.2i\n',dv(1:3));
  dmm = squeeze(D(it,:,:));
  Di = dmm.*LGr; % Greenland mask applied
  dmm = squeeze(Rt(it,:,:));
  Rti = dmm.*LGr; % Gr tundra
  dmm = squeeze(Rg(it,:,:));
  Rgi = dmm.*LGr; % Greenland meltwater
  FWF(it,1) = nansum(nansum(Di))+...
      nansum(nansum(Rti))+...
      nansum(nansum(Rgi));
  iD(it,1) = nansum(nansum(Di)); % total ice discharge, km3/mo
  iRt(it,1)= nansum(nansum(Rti));
  iRg(it,1)= nansum(nansum(Rgi));
  
% FWF by regions
  for ik=1:iR
    IN=IJ(ik).IN;
    IJ(ik).Disch(it,1)=nansum(Di(IN));   % solid discharge
    IJ(ik).Rtdr(it,1)=nansum(Rti(IN));   % tundra runoff
    IJ(ik).Rmlt(it,1)=nansum(Rgi(IN));   % meltwater
    IJ(ik).FWF(it,1)=nansum(Di(IN))+...
	nansum(Rti(IN))+...
	nansum(Rgi(IN));                   % total FW flux
  end
  
end  

%keyboard

dmm=0;
for ik=1:iR
  fwf(ik,:)=IJ(ik).FWF;
  dmm=dmm+IJ(ik).FWF;
end




TY = [0:nrc-1]/12+DV(1,1) ;
iyr=36;  % reference mean 1958-1993, anomaly plot from 1993

% Do yearly fluxes:
YR = [DV(1,1):DV(end,1)];
nyr = nrc/12;
for ik=1:iR
  iD=IJ(ik).Disch;
  iRt=IJ(ik).Rtdr;
  iRg=IJ(ik).Rmlt;
  A = reshape(iD,[12 nyr]);
  sD = sum(A);
  A = reshape(iRt,[12 nyr]);
  sRt = sum(A);
  A = reshape(iRg,[12 nyr]);
  sRg = sum(A);
  FWFt = sD+sRt+sRg;

  IJ(ik).Disch_yr=sD;
  IJ(ik).Rtdr_yr=sRt;
  IJ(ik).Rmlt_yr=sRg;
  IJ(ik).FWF_yr=FWFt;
end

A=reshape(FWF,[12 nyr]);
FWFyr=sum(A);


lnm=[];
for ik=1:iR
  lnm{ik}=IJ(ik).Name;
end




figure(2); clf
axes('Position',[0.09 0.58 0.85 0.35]);
hold on;
for ik=1:iR
  sD=IJ(ik).Disch_yr;
  sRt=IJ(ik).Rtdr_yr;
  sRg=IJ(ik).Rmlt_yr;
  clr=CLR(ik,:);
  plot(YR,sD,'Color',clr);
end
set(gca,'tickdir','out',...
	'xlim',[YR(iyr) YR(end)],...
	'xtick',[1950:2016],...
	'xgrid','on',...
	'ygrid','on');
title('Greenland FWF: Solid Disch, km3/yr');
pl=legend(lnm);
set(pl,'Position',[0.75 0.3 0.2 0.16]);

bottom_text(btx,'pwd',1,'Position',[0.08 0.1 0.4 0.1]);


figure(3); clf
axes('Position',[0.09 0.58 0.85 0.35]);
hold on;
for ik=1:iR
  sRt=IJ(ik).Rtdr_yr;
  clr=CLR(ik,:);
  plot(YR,sRt,'Color',clr);
end
set(gca,'tickdir','out',...
	'xlim',[YR(iyr) YR(end)],...
	'xtick',[1950:2016],...
	'xgrid','on',...
	'ygrid','on');
title('Greenland FWF: Tundra Runoff, km3/yr');
pl=legend(lnm);
set(pl,'Position',[0.75 0.3 0.2 0.16]);

bottom_text(btx,'pwd',1,'Position',[0.08 0.1 0.4 0.1]);


figure(4); clf
axes('Position',[0.09 0.58 0.85 0.35]);
hold on;
for ik=1:iR
  sRg=IJ(ik).Rmlt_yr;
  clr=CLR(ik,:);
  plot(YR,sRg,'Color',clr);
end
set(gca,'tickdir','out',...
	'xlim',[YR(iyr) YR(end)],...
	'xtick',[1950:2016],...
	'xgrid','on',...
	'ygrid','on');
title('Greenland FWF: Meltwater, km3/yr');
pl=legend(lnm);
set(pl,'Position',[0.75 0.3 0.2 0.16]);

bottom_text(btx,'pwd',1,'Position',[0.08 0.1 0.4 0.1]);



figure(5); clf
axes('Position',[0.09 0.58 0.85 0.35]);
hold on;
for ik=1:iR
  Ft=IJ(ik).FWF_yr;
  clr=CLR(ik,:);
  plot(YR,Ft,'Color',clr);
end
plot(YR,FWFyr,'Color',[0 0 0]);
lnm{ik+1}=sprintf('%s ','GreenlTot');
set(gca,'tickdir','out',...
	'xlim',[YR(iyr) YR(end)],...
	'xtick',[1950:2016],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',12);
title('Greenland FWF total, km3/yr');
pl=legend(lnm);
set(pl,'Position',[0.75 0.3 0.2 0.16],'Fontsize',14);

bottom_text(btx,'pwd',1,'Position',[0.08 0.1 0.4 0.1]);



% Plot anomaly
fwf0=mean(FWFyr(1:36));
dFWF=FWFyr-fwf0;

figure(6); clf
axes('Position',[0.09 0.58 0.85 0.35]);
hold on;
for ik=1:iR
  Ft=IJ(ik).FWF_yr;
  F0=mean(Ft(1:36));
  dF=Ft-F0;
  clr=CLR(ik,:);
  plot(YR,dF,'Color',clr);
  plg{ik}=sprintf('%s (%3i)',IJ(ik).Name,round(F0));
end
plot(YR,dFWF,'Color',[0 0 0]);
plg{ik+1}=sprintf('%s (%3i)','GreenlTot',round(fwf0));

set(gca,'tickdir','out',...
	'xlim',[YR(iyr) YR(end)],...
	'xtick',[1950:2016],...
	'xgrid','on',...
	'ygrid','on');
title('Greenland FWF anomaly (1958-1993 ref), km3/yr');
pl=legend(plg);
set(pl,'Position',[0.6 0.3 0.25 0.18],'Fontsize',14);

%bottom_text(btx,'pwd',1,'Position',[0.08 0.1 0.4 0.1]);

% Plot anomaly as fraction of total overall FWF anomaly
%figure(7); clf
%axes('Position',[0.09 0.58 0.85 0.35]);
%hold on;

clear plg
for ik=1:iR
  Ft=IJ(ik).FWF_yr;
  F0=mean(Ft(1:36));
  dF=Ft-F0;
  DF(ik,:)=dF; % anomaly by regions
%  clr=CLR(ik,:);
%  plot(YR,dF./dFWF,'Color',clr);
  plg{ik}=sprintf('%s (%3i)',IJ(ik).Name,round(mean(dF(36:end))));
end
plg{ik+1}=sprintf('%s (%3i)','GreenlTot',round(mean(dFWF(36:end))));


axes('Position',[0.1 0.3 0.25 0.18]);
text(0.1, 0.5,plg,'Fontsize',14);
text(0.1,0.75,'Mean Anomalies km3/yr for 1993-2016','Fontsize',14);
set(gca,'xlim',[0 0.43],...
	'ylim',[0.25 0.7],...
	'visible','off');

%set(gca,'tickdir','out',...
%	'xlim',[YR(36) YR(end)],...
%	'xtick',[1950:5:2016],...
%	'xgrid','on',...
%	'ygrid','on');
%title('Greenland FWF anomaly (1958-1993 ref), km3/yr');
%pl=legend(plg);
%set(pl,'Position',[0.75 0.3 0.2 0.16]);

bottom_text(btx,'pwd',1,'Position',[0.08 0.1 0.4 0.1]);










