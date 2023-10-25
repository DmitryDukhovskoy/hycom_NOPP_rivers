% Compare river runoff in the 
% hycom river_*.[ab] files
% between experiments with runoff climatology
% expt 11.0 and with NCAR rivers
% also
% Compare with AOMIP climatology
% In forcing HYCOM files, Units:  m/s
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

PTH.data  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.data2 = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/riversNCAR/';
PTH.topo  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%flt=[pth,'regional.grid.b'];

fmat = sprintf('%sncar_rivers_Arctic_1993-2015.mat',PTH.data);
fprintf('Loading %s\n',fmat);
load(fmat);

YY=1993;
%flriva=sprintf('%srivers_09.a',PTH.data);
%flrivb=sprintf('%srivers_09.b',PTH.data);
%flriva=sprintf('%srivers_09_corrected.a',PTH.data);
%flrivb=sprintf('%srivers_09_corrected.b',PTH.data);
%flriva=sprintf('%srivers_09_Greenland_%i.a',PTH.data,YY);
%flrivb=sprintf('%srivers_09_Greenland_%i.b',PTH.data,YY);
flriva=sprintf('%srivers_11.a',PTH.data);  % expt 11.0 clim. rivers, no Gr.
flrivb=sprintf('%srivers_11.b',PTH.data);
flriva2=sprintf('%srivers_11_NCAR_Gr_%4.4i.a',PTH.data2,YY); % NCAR riv+Green.
flrivb2=sprintf('%srivers_11_NCAR_Gr_%4.4i.b',PTH.data2,YY); % NCAR riv+Green.
flriva3=sprintf('%srivers_11_NCAR_Gr_%4.4i.a',PTH.data2,2015); % NCAR riv+Green.
flrivb3=sprintf('%srivers_11_NCAR_Gr_%4.4i.b',PTH.data2,2015); % NCAR riv+Green.
fltopo=sprintf('%sdepth_ARCc0.08_09.nc',PTH.topo);

%fltxt=[pth1,fvrs,'.b'];
%flbth=[pth1,fvrs,'.a'];
%flrv=[pth,'rivers_03.a'];

% Get topo and grid:
HH   = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn]= size(HH);
[m,n]= size(HH);

[DX,DY]=sub_dx_dy(LON,LAT);
ACell=DX.*DY;

% River data with NCAR runoff and correct 
% river locations:
%pthmat = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
%fmat   = sprintf('%sncar_rivers_Arctic_1993-2015.mat',pthmat);
%load(fmat); % <-- RVR

% Arctic Domain with Greenland
IJarc=[  380         353
         505         355
	 705         376
	 841         558
	 937         575
        1184         449
        1298         686
        1483         612
        1594         664
        1594        1323
	1235        1916
         443        1916
         188        1398
          51        1005
          44         531
         147         234];
IJarc(end+1,:)=IJarc(1,:);
[II,JJ]=meshgrid((1:nn),(1:mm));
dmm = inpolygon(II,JJ,IJarc(:,1),IJarc(:,2));
IN = find(dmm==1);
OUT = find(dmm==0);

% Greenland: 
IJgr =[  602         372
         912         645
         982         960
         917        1080
         785        1076
         583        1084
         558         923
         503         488
         539         401];
IJgr(end+1,:)=IJgr(1,:);
dmm = inpolygon(II,JJ,IJgr(:,1),IJgr(:,2));
INg = find(dmm==1);
OUTg = find(dmm==0);






% Reading rivers & bathymetry:
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);

RV = sub_define_rvr_clim;
NR=length(RV);

% Rivers:
riv_fid=fopen(flriva,'r','ieee-be');
rivfid2=fopen(flriva2,'r','ieee-be');
rivfid3=fopen(flriva3,'r','ieee-be');

for k=1:12  % time - 12mo
  mo = k;
  fprintf('check_river_runoff %i:    Reading month %i\n',YY,k);

  A  = fread(riv_fid,IJDM,'float32'); % read 2D field
  dm1= fread(riv_fid,npad,'float32');  % Padding = size(toto)

  I=find(A>1e10);
  A(I)=NaN;
  A=reshape(A,IDM,JDM)';
  A(A==0)=nan;


  A2 = fread(rivfid2,IJDM,'float32'); % read 2D field
  dm1= fread(rivfid2,npad,'float32');  % Padding = size(toto)

  A3 = fread(rivfid3,IJDM,'float32'); % read 2D field
  dm1= fread(rivfid3,npad,'float32');  % Padding = size(toto)

  I=find(A2>1e10);
  A2(I)=NaN;
  A2=reshape(A2,IDM,JDM)';
  A2(A2==0)=nan;
  
  I=find(A3>1e10);
  A3(I)=NaN;
  A3=reshape(A3,IDM,JDM)';
  A3(A3==0)=nan;
  
% Convert m/s -> m3/s  
  RFlx  = A.*ACell;   % climatology river 
  RFlx2 = A2.*ACell; % NCAR rivers year YY
  RFlx3 = A3.*ACell; % NCAR rivers year 2015
  
% Check total river runoff, m3/s, m3/mo
% Arctic Ocean+ Greenland:
  dmm = A;
  dmm(OUT) = nan;
  II = find(dmm>0);
  R = sum(A(II).*ACell(II));
  dmm = A;
  dmm(OUTg) = nan;
  II = find(dmm>0);
  Rgr = sum(A(II).*ACell(II));
  clear dmm

  Rarc = R-Rgr;
  RarcM(k,1) = Rarc*3600*24*30*1e-9; % Arctic only, no Gr., km3/mo  
%  RgrM  = Rgr*3600*24*30*1e-9;  % Greenland, km3/mo it is zero in old experiment

% New experiment with NCAR rivers, year=YY
  dmm = A2;
  dmm(OUT) = nan;
  II = find(dmm>0);
  R = sum(A2(II).*ACell(II));
  dmm = A2;
  dmm(OUTg) = nan;
  II = find(dmm>0);
  Rgr = sum(A2(II).*ACell(II));
  clear dmm
  Rarc = R-Rgr;
  RarcM2(k,1) = Rarc*3600*24*30*1e-9; % Arctic only, no Gr., km3/mo
  RgrM2(k,1)  = Rgr*3600*24*30*1e-9;  % Greenland, km3/mo it is zero in old experiment
  
% New experiment with NCAR rivers, 2015
  dmm = A3;
  dmm(OUT) = nan;
  II = find(dmm>0);
  R = sum(A3(II).*ACell(II));
  dmm = A3;
  dmm(OUTg) = nan;
  II = find(dmm>0);
  Rgr = sum(A3(II).*ACell(II));
  clear dmm
  Rarc = R-Rgr;
  RarcM3(k,1) = Rarc*3600*24*30*1e-9; % Arctic only, no Gr., km3/mo
  RgrM3(k,1)  = Rgr*3600*24*30*1e-9;  % Greenland, km3/mo it is zero in old experiment


  
% Get Runoff by rivers
  for ir=1:NR
    rnm=RV(ir).Name;
    i1=RV(ir).I(1);
    i2=RV(ir).I(2);
    j1=RV(ir).J(1);
    j2=RV(ir).J(2);
    dmm=RFlx(j1:j2,i1:i2);
    Rtot=nansum(nansum(dmm));
    RV(ir).R_hycom_m3s(k)=Rtot;
    RV(ir).rname = rnm;
    
    [iC,jC,Qr] = sub_get_loc_rivername(rnm,LAT,LON,HH,YY,mo);    
    RV2(ir).rname = rnm;
    RV2(ir).R_hycom_m3s(k) = Qr;

    [iC,jC,Qr] = sub_get_loc_rivername(rnm,LAT,LON,HH,2015,mo);    
    RV3(ir).rname = rnm;
    RV3(ir).R_hycom_m3s(k) = Qr; % river runoff in 2015
  end

  
end

%keyboard

% Plot monthly runoff for the whole Arctic:
figure(1); clf;
axes('Position',[0.08 0.55 0.88 0.37]);
hold on;
plot(RarcM,'b-','Linewidth',2);
plot(RarcM2,'r-','Linewidth',2);
plot(RarcM3,'g-','Color',[0 0.8 0.4],'Linewidth',2);
set(gca,'tickdir','out',...
	'xlim',[0.9 12.1],...
	'xtick',[1:12],...
	'ylim',[0 1.12*max(RarcM3)],...
	'xgrid','on',...
	'ygrid','on');
qsm = sum(RarcM);
qsm2= sum(RarcM2);
qsm3= sum(RarcM3);
LG=legend('Clim ',sprintf('NCAR %i',YY),'NCAR 2015');
tts = sprintf('ARCc Arctic Runoff, km3/mo');
title(tts);

spp{1}=sprintf('Clim: %5.1f km3/yr',qsm);
spp{2}=sprintf('%i:  %5.1f km3/yr',YY,qsm2);
spp{3}=sprintf('2015: %5.1f km3/yr',qsm3);
text(1.5,700,spp);

axes('Position',[0.08 0.08 0.88 0.37]);
hold on;
plot(RgrM2,'r-','Linewidth',2);
plot(RgrM3,'g-','Color',[0 0.8 0.4],'Linewidth',2);
set(gca,'tickdir','out',...
	'xlim',[0.9 12.1],...
	'xtick',[1:12],...
	'ylim',[0 1.12*max(RgrM3)],...
	'xgrid','on',...
	'ygrid','on');
qgr2= sum(RgrM2);
qgr3= sum(RgrM3);
tts = sprintf('Greenland Runoff, km3/mo');
title(tts);
srr{1}=sprintf('%i:  %5.1f km3/yr',YY,qgr2);
srr{2}=sprintf('2015: %5.1f km3/yr',qgr3);
text(1.5,200,srr);
txtbtm='compare_NCAR_runoff.m';
bottom_text(txtbtm,'pwd',1);



% Plot individual rivers
mdays = [31,28,31,30,31,30,31,31,30,31,30,31];
figure(2); clf;
for ir=1:NR
  nm=RV(ir).Name;
  rcl=RV(ir).Clim_m3s;
  rhycom=RV(ir).R_hycom_m3s;
  rhycom2=RV2(ir).R_hycom_m3s;
  rhycom3=RV3(ir).R_hycom_m3s;
  
  rtot = rhycom*mdays'*3600*24*1e-9; % km3/yr
  
  subplot(3,3,ir);
%  plot(rcl,'c--','linewidth',2);
  hold on;
  plot(rhycom,'b','linewidth',2);
  plot(rhycom2,'r','linewidth',2);
  plot(rhycom3,'g','Color',[0 0.8 0.4],'linewidth',2);
  stt=sprintf('%s, m3/s, Rtot=%4.1f km3/yr',nm,rtot);
  title(stt,'Fontsize',9);
  set(gca,'tickdir','out',...
	  'xlim',[1 12],...
	  'xtick',[1:12],...
	  'xgrid','on',...
	  'ygrid','on');
end;


bottom_text(txtbtm,'pwd',1);







