% Analysis of FW/Vol fluxes 
% in Fram Strait, on Greenland shelf
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

Sref=34.9;
YR1=2002;
YR2=2016;


pthdat='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_Myers/FRAM/';

btx = 'anls_trnsp_Fram_myers.m';

fout=sprintf('%smonthly_FW_vol_fluxes.mat',pthdat);
fprintf('Loading %s\n',fout);
load(fout);

ixW=VF1D.ixW;
ixE=VF1D.ixE;

% 1D depth-integrated fluxes:
% average by years
nyr=length(VF1D);
mV=zeros(1,12);
mV1=mV;
mV2=mV;
mV3=mV;
mF=mV;
mF1=mV;
mF2=mV;
mF3=mV;
% FW fluxes:
for ik=1:nyr
  for imo=1:12
    aa=VF1D(ik).FWF_1D(imo,:)*1e-3;  % mSv
    mF(imo)=mF(imo)+nansum(aa);
    mF1(imo)=mF1(imo)+nansum(aa(1:ixW-1));  % Gr Shelf - F17 segment S1
    mF2(imo)=mF2(imo)+nansum(aa(ixW:ixE-1));
    mF3(imo)=mF3(imo)+nansum(aa(ixE:end));
  end
end
mF=mF/nyr;
mF1=mF1/nyr;
mF2=mF2/nyr;
mF3=mF3/nyr;

FF=-[abs(mF1)', abs(mF2)', abs(mF3)'];  % keep both +/- fluxes
pp=mF2./(abs(mF1)+abs(mF2)+abs(mF3)); % fraction of observed FW flux

cmp1=[0. 0. 0.4; ...
	  0. 0.4 0.8; ...
	  0  0.8 1];


for ik=1:nyr
  for imo=1:12
    aa=VF1D(ik).VolF_1D(imo,:)*1e-6;  % Sv
    mV(imo)=mV(imo)+nansum(aa);
    mV1(imo)=mV1(imo)+nansum(aa(1:ixW-1));
    mV2(imo)=mV2(imo)+nansum(aa(ixW:ixE-1));
    mV3(imo)=mV3(imo)+nansum(aa(ixE:end));
  end
end
mV=mV/nyr;
mV1=mV1/nyr;
mV2=mV2/nyr;
mV3=mV3/nyr;

FVV=-[abs(mV1)',abs(mV2)',abs(mV3)'];
ppv=abs(mV2)./(abs(mV1)+abs(mV2)+abs(mV3)); % fraction of observed volume

cmp2=[0 0.3 0; ...
      0 0.7 0.3;...
      0 1 0.6];

    
figure(10); clf;
axes('Position',[0.1 0.55 0.8 0.35]);
hold on;
HB=bar(FF,'stacked');
colormap(cmp1);
for ib=1:3
  HB(ib).EdgeColor='none';
  HB(ib).FaceColor=cmp1(ib,:);
  HB(ib).BarWidth=0.95;
end
for imo=1:12
  stx=sprintf('%3.2f',pp(imo));
  text(imo-0.3,0.7*mF(imo),stx,'Fontsize',14);
end

set(gca,'tickdir','out',...
	'Fontsize',14,...
	'xlim',[0.5 12.5],...
	'ylim',[-150 0],...
	'xtick',[1:12],...
	'ytick',[-150:25:0]);


hcb=colorbar;
set(hcb,'Position',[0.91 0.55 0.014 0.35],...
	'Ticks',[1/6 1/2 5/6],...
	'TickLabels',{'S1','S2','S3'},...
	'Fontsize',14);
stl=sprintf('NEMO FW Flux %i-%i',YR1,YR2);
title(stl);

bottom_text(btx,'pwd',1,'position',[0.08 0.4 0.4 0.05]);


% Volume flux - Note there is 
% an error in the flux calculation
% need to recompute
% transport should be ~-2 Sv
figure(11); clf;
axes('Position',[0.1 0.55 0.8 0.35]);
hold on;
HB=bar(FVV,'stacked');
colormap(cmp2);

for ib=1:3
  HB(ib).EdgeColor='none';
  HB(ib).FaceColor=cmp2(ib,:);
  HB(ib).BarWidth=0.95;
end
for imo=1:12
  stx=sprintf('%3.2f',ppv(imo));
  text(imo-0.3,0.7*mV(imo),stx,'Fontsize',14);
end

set(gca,'tickdir','out',...
	'Fontsize',14,...
	'xlim',[0.5 12.5],...
	'ylim',[-14 0],...
	'xtick',[1:12],...
	'ytick',[-20:2:0]);


hcb=colorbar;
set(hcb,'Position',[0.91 0.55 0.014 0.35],...
	'Ticks',[1/6 1/2 5/6],...
	'TickLabels',{'S1','S2','S3'},...
	'Fontsize',14);
stl=sprintf('NEMO VolTrt %i-%i',YR1,YR2);
title(stl);

bottom_text(btx,'pwd',1,'position',[0.08 0.4 0.4 0.05]);








