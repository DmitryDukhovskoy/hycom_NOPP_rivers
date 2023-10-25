% Compute volume integrated mass of the tracer
% within the specified domain IN
% If IN is empty - the whole domain
% nTr - tracer # that is integrated
% Acell - grid cell area, m2
% plr - layer to integrate, plr = 0 - all layers
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.04';
expt = 011;  


plr   = 0;  % layer to integrate, =0 - over whole depth, use this
s_plt = 1;  % plot spatial mean Mass tracer 2D
s_mat = 1; % =1 - save mat file flag
           % =2 - load saved fields

nTr=1;   % tracer to plot
TrMn=1e-10; % Threshold value to plot
rg = 9806;
fprintf('Tracer #: %i, Threhold value: %8.5d\n',nTr,TrMn);


rg = 9806;

YRPLT=[2006,172];

% Experiments:
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);
pthmat8 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/110/data_mat/');


ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

ip=1;
yr   = YRPLT(ip,1);
iday = YRPLT(ip,2);
dnmb = datenum(yr,1,1)+iday-1;
DV   = datevec(dnmb);
pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data/%i/',yr);  
fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
    
fprintf('\n:::: Analyzing %i/%2.2i/%2.2i   ::::\n\n',DV(1:3));

if s_mat<2
  fprintf('Computing vol-integrated Tracer mass ...\n');

  if plr>0
    [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr,'r_layer',plr);
    F(F>1e6)=nan;
    Tr=F(plr,:,:);

    [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
  %  F=squeeze(F(plr,:,:));
    F=F./rg;
    F(F>1e10)=0;
    F(F<1e-2)=0;
    dH=F; 

  else
    [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
    F(F>1e6)=nan;
    Tr=F;

    [F,n,m,l] = read_hycom(fina,finb,'thknss');
    F=F./rg;
    F(F>1e10)=0;
    F(F<1e-2)=0;
    dH=F; 

    TRI = sub_intgr_tr(dH,Tr,HH,Acell);

    if s_mat==1
      fmat = sprintf('%s%s_%3.3i_VolIntgr_%i%2.2i%2.2i.mat',...
		     pthmat,regn,expt,DV(1:3));
      fprintf('Saving %s\n',fmat);
      save(fmat,'TRI');
    end

  end
  
else
  fmat = sprintf('%s%s_%3.3i_VolIntgr_%i%2.2i%2.2i.mat',...
		 pthmat,regn,expt,DV(1:3));
  fprintf('Loading %s\n',fmat);
  load(fmat);
end

if s_plt==0, fprintf('No plotting flag f_plt=%i\n',s_plt); end;
if s_plt==0; return; end;

Mtot = TRI.OverallMass_kg;
Vtot = TRI.OverallVol_m3;
Mtr  = TRI.DpthIntgrMass_kg_m2;
Ctot = Mtot./Vtot; % mean trac. concentration
[mm,nn] = size(Mtr);

% Mean vol-integrated:
fprintf('Overall tr mass = %18.15d GT\n',Mtot*1e-12); % kg->tonn->GT
fprintf('Mean overall conc. = %18.15d kg/m3\n',Ctot);

% 2D mean conc per 1 m2:
Mtr(Mtr<=0)=nan;
lTr = log(Mtr);

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;
nf = 1;
ifx=max(strfind(fina,'/'));

stl=sprintf('%s Depth-intgr. Mass kg/m2, Log2, Tr=%i',...
	    fina(ifx+1:end),nTr);
sub_plot_tracers(lTr,nf,HH,xlim1,...
		 xlim2,ylim1,ylim2,...
		 LON,LAT,stl,nTr,...
		 'c1',2,'c2',8);

stx{1} = sprintf('Overall Tr Mass, GT = %18.15d',Mtot*1e-12);
stx{2} = sprintf('Mean DIntgr  conc. = %10.8d kg/m3\n',Ctot);
axes('Position',[0.27 0.83 0.27 0.04]);
text(0.05,0.5,stx);
set(gca,'xtick',[],'ytick',[],'box','off');

txtb='vol_intgr_tracr004.m';
bottom_text(txtb,'pwd',1,'position',[0.02 0.05 0.8 0.1]);

% Mass tracers by layers:
dzm=TRI.Mean_LThkn;
zzm=-cumsum(dzm);
Mv =TRI.Vertical_Mean_Mass_kg;
Mv = Mv*1e-12; % kg->tonn (1e-3)->GT (1e-9)




% Load ARCc0.08:
% Note: need to run vol_intgr_trcr008.m for the same date
% to generate *.mat file
fmat8 = sprintf('%sARCc0.08_110_VolIntgr_%i%2.2i%2.2i.mat',...
		pthmat8,DV(1:3));
if (~exist(fmat),'file'),
  fprintf('Missing for ARCc0.08 %s\n',fmat);
  fprintf('need to run vol_intgr_trcr008.m for the same date\n');
end

A=load(fmat8);
TRI8=A.TRI;
clear A
dzm8=TRI8.Mean_LThkn;
zzm8=-cumsum(dzm);
Mv8 =TRI8.Vertical_Mean_Mass_kg;
Mv8 = Mv8*1e-12; % GT


figure(2); clf;
axes('position',[0.08 0.08 0.35 0.85]);
plot(Mv,zzm);
hold on;
plot(Mv8,zzm8,'r');
stl=sprintf('%s Layer-intgr. Mass GT, Tr=%i',...
	    fina(ifx+1:end),nTr);
title(stl,'Interpreter','none');
set(gca,'tickdir','out',...
	'xgrid','on','ygrid','on',...
	'xlim',[0 1.02*max([max(Mv8),max(Mv)])],...
	'xtick',[0:20:300],...
	'ylim',[-1120 0],...
	'ytick',[-2000:100:0]);

legend('ARCc0.04','ARCc0.08');

dM=Mv-Mv8;
mxd=max(abs(dM));
axes('position',[0.5 0.08 0.35 0.85]);
plot(dM,zzm8);
title('ARCc0.04-ARCc0.08');
set(gca,'tickdir','out',...
	'xgrid','on','ygrid','on',...
	'xlim',[-1.01*mxd 1.01*mxd],...
	'xtick',[-100:10:100],...
	'ylim',[-1120 0],...
	'ytick',[-2000:100:0]);

txtb='vol_intgr_tracr004.m';
bottom_text(txtb,'pwd',1,'position',[0.02 0.03 0.8 0.05]);

