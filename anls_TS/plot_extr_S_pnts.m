% Plot  Extracted S at some locations
% 
% Monthly mean T/S averaged within the layers
% for 2 experiments (with & without Greenland runoff)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.08';
pfld = 'salin';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
%expt = 110; % no Greenland runoff  
%expt = 112;  % Greenland runoff

LRS = [0; -52; -152; -300; -500; -1000];
nlr = length(LRS)-1;


for ixp=1:2
  if ixp==1
    expt=110;
  else
    expt=112;
  end
 
  pthmat=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
  fmat=sprintf('%sarc08_%3.3i_%s_sections.mat',pthmat,expt,pfld);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  if ixp==1
    S1=SCT;	
  else
    S2=SCT;
  end
  
end

% Butterworth filter
% output freq. 1/12 mo^-1
% 12 mo cutoff = 12/12=1 -> in Matlab Wn=1/6
% 
Wn = 1/6; % cutoff freq 6 mo: 2/6, 1yr=1/6
[Bf,Af] = butter(9,Wn,'low');


nsct=length(SCT);

TM=S1(1).TM;
DV=datevec(TM);
yrs=[DV(1,1):1/12:DV(end,1)+0.99];
nyrs=DV(end,1)-DV(1,1)+1;
YY=[DV(1,1):DV(end,1)];

CLR=[0 0 0;...
     0 0.3 1;...
     0 0.7 0.2; ...
     0.8 0 0;...
     1 0.7 0.3];

btx = 'plot_extr_S_pnts.m';
for ik=1:nsct
  ss1=S1(ik).SS;
  ss2=S2(ik).SS;
  nm =S2(ik).Name;
  
  figure(ik); clf;
  axes('position',[0.09 0.38 0.87 0.55]);
  hold on;
  
  for ilr=1:nlr
    s1=ss1(ilr,:);
    s2=ss2(ilr,:);
    ds=s2-s1;
    I=find(isnan(ds));
    mnn=nanmean(ds);
    ds(I)=mnn;
    clr=CLR(ilr,:);
% Filter - 3 mo 
    dmm = filtfilt(Bf,Af,ds);
    plot(yrs,dmm,'Color',clr,'linewidth',1.8);
% Plot annual mean ds
%    dmm=reshape(ds,12,nyrs);
%    mnds=nanmean(dmm); % mean
%    mnds=min(dmm);             % min ds
%    plot(yrs,ds,'Color',clr,'linewidth',1.8);
%    plot(YY,mnds,'Color',clr,'linewidth',1.8);
%    plot(s1,'Color',clr,'linewidth',1.8);
  end
  
  stl = sprintf('arc0.08 dltS btw expts, %s, July',nm);
  title(stl,'Interpreter','none');
  
  set(gca,'tickdir','out',...
	  'xtick',[1990:1:2016],...
	  'xlim',[1993 2016],...
	  'xgrid','on',...
	  'ygrid','on');
	  
  axes('position',[0.1 0.1 0.2 0.15]);
  hold on
  for ilr=1:nlr
    clr=CLR(ilr,:);
    plot([0 1],[nlr-ilr nlr-ilr],'Color',clr,'linewidth',1.8);
    zz = LRS(ilr+1);
    stl = sprintf('%4i m',zz);
    text(1.1,nlr-ilr,stl,'Fontsize',12);
  end
  set(gca,'xlim',[-0.2 3.5],...
	  'ylim',[-0.2 nlr+0.3],...
	  'visible','off');
  
  bottom_text(btx,'pwd',1);
  
end



  



