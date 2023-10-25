% The updated Greenland runoff data set 
% From J. Bamber
% Dear Dmitry,
%
% Here is the link to the new FWF data which is in a netcdf file:
% 
% https://drive.google.com/file/d/1uDryUZqD0g4w6_3BH99-l4ezY7hljvqj/view?usp=sharing
%
%I also attach the latest version of paper I. 
% Figure 3 in particular shows the individual flux 
% components and when you read the data, 
% you should make sure your time series has the same values as in the plot. 
% To do that, you will need to sum the monthly data from Jan-Dec. 
% Remember that the dashed lines are annual 
% and the solid lines are 5 year running means.

%Also, it is important to note that, in the paper, 
% we only plot and include tundra runoff from Greenland 
% but we include tundra runoff from all the other 
% land areas shown in Fig 4. 
% To exclude these (if you want to) you will need to use the 
% Greenland land/sea mask that is provided in the ncdf file.
%
% Hope that all makes sense. Any questions let me know.
%
%Regards, Jonanthan
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

pthG = '/Net/ocean/ddmitry/arctic_AOregimes/data/GreenlandRunoffv3/';
pthmat = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat/';

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

% Total FW flux = D+Rg+Rt (no CAA, Svalbard here);
% Reproduce Bamber's Fig.3:
clear FWF iD iR*
for it=1:nrc
  dv = datevec(TM(it));
  fprintf('Greenland Runoff, %i/%2.2i/%2.2i\n',dv(1:3));
  dmm = squeeze(D(it,:,:));
  Di = dmm.*LGr; % Greenland mask
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
end  

TY = [0:nrc-1]/12+DV(1,1) ;

% Plot fraction of solid discharge vs runoff
fr = iD./(iRt+iRg);
Is = find(DV(:,2)>=5 & DV(:,2)<9); % summer
Iw = find(DV(:,2)>=9 | DV(:,2)<5); % winter
bxs=[0:0.5:10];
hxs=hist(fr(Is),bxs);
bxw=[0:50:700];
hxw=hist(fr(Iw),bxw);

figure(1); clf;
subplot(2,2,1); 
boxplot(fr(Is));
title('Greenland SolidDisch/(Melt+Tundra), Summer');
subplot(2,2,2); 
boxplot(fr(Iw));
title('Fraction, Winter');

btx = 'read_GreenlandRunoffv3.m';
bottom_text(btx,'pwd',1,'Position',[0.02 0.5 0.9 0.1]);

% Do yearly fluxes:
YR = [DV(1,1):DV(end,1)];
nyr = nrc/12;
A = reshape(iD,[12 nyr]);
sD = sum(A);
A = reshape(iRt,[12 nyr]);
sRt = sum(A);
A = reshape(iRg,[12 nyr]);
sRg = sum(A);

FWFt = sD+sRt+sRg;


figure(2); clf
axes('Position',[0.09 0.58 0.85 0.35]);
plot(YR,sD);
hold on;
plot(YR,sRt,'g');
plot(YR,sRg,'r');
set(gca,'tickdir','out',...
	'xlim',[YR(1) YR(end)],...
	'xtick',[1950:5:2016],...
	'xgrid','on',...
	'ygrid','on');
title('Greenland FW Flux Components, km3/yr');
legend('Disch','Rtundra','Rmelt');

axes('Position',[0.09 0.1 0.85 0.35]);
plot(YR,FWFt);
set(gca,'tickdir','out',...
	'xlim',[YR(1) YR(end)],...
	'xtick',[1950:5:2016],...
	'xgrid','on',...
	'ygrid','on');
title('Greenland total FW Flux, km3/yr');

btx = 'read_GreenlandRunoffv3.m';
bottom_text(btx,'pwd',1);

% save integrated fluxes
if s_mat ==1
  fmat = sprintf('%sGreenlFWFv3.mat',pthmat);
  GrF.TM = TM;
  GrF.Years_mo = TY;
  GrF.Years = YR;
  GrF.IceDisch_km3mo = iD;
  GrF.Tundra_Runoff_km3mo = iRt;
  GrF.GrMeltWater_km3mo = iRg;
  GrF.TotalFWF_km3mo = FWF;
  GrF.IceDisch_km3yr = sD;
  GrF.Tundra_Runoff_km3yr = sRt;
  GrF.GrMeltWater_km3yr = sRg;
  GrF.FWFtotal_km3yr = FWFt;
  
  fprintf('Saving %s\n',fmat);
  save(fmat,'GrF');
end




