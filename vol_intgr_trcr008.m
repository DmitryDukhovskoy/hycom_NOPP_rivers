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

regn = 'ARCc0.08';
expt = 110;  

nTr   = 1   % tracer to plot
plr   = 0;  % layer to integrate, =0 <- over whole depth, use this
s_f1  = 0;  % plot spatial mean Mass tracer 2D
s_f2  = 0;  % vertical prof., individual years
s_f3  = 1;  % vertical prof. all years
s_fig = 0; % saves all active figures

s_mat = 2; % =1 - save mat file flag
           % =2 - load saved fields

TrMn=1e-10; % Threshold value to plot
rg = 9806;
fprintf('Tracer #: %i, Threhold value: %8.5d\n',nTr,TrMn);


rg = 9806;


%YRPLT=[2006,172];
YRPLT=[1998,172;...
       2003,172;...
       2008,172;...
       2015,365];
nyy = size(YRPLT,1);

np=size(YRPLT,1);

% Experiments:
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);


%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2
MVOL = [];

for ip=1:nyy;
  yr   = YRPLT(ip,1);
  iday = YRPLT(ip,2);
  dnmb = datenum(yr,1,1)+iday-1;
  DV   = datevec(dnmb);
  
  fmat = sprintf('%s%s_%3.3i_VolIntgr_Tr%2.2i_%i%2.2i%2.2i.mat',...
		 pthmat,regn,expt,nTr,DV(1:3));
  
  
  %pthbin = sprintf('/nexsan/hycom/ARCc0.08_011/data/%i/',yr);  
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
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
%	fmat = sprintf('%s%s_%3.3i_VolIntgr_Tr%2.2i_%i%2.2i%2.2i.mat',...
%		       pthmat,regn,expt,nTr,DV(1:3));
	fprintf('Saving %s\n',fmat);
	save(fmat,'TRI');
      end

    end
  else
%    fmat = sprintf('%s%s_%3.3i_VolIntgr_%i%2.2i%2.2i.mat',...
%		   pthmat,regn,expt,DV(1:3));
    fprintf('Loading %s\n',fmat);
    load(fmat);

    
  end

  Mtot = TRI.OverallMass_kg;
  Vtot = TRI.OverallVol_m3;
  Mtr  = TRI.DpthIntgrMass_kg_m2;
  Ctot = Mtot./Vtot; % mean trac. concentration
  [mm,nn] = size(Mtr);

% Mass tracers by layers:
  dzm=TRI.Mean_LThkn;
  zzm=-cumsum(dzm);
  Mv =TRI.Vert_LrMass_kg;
  Mv = Mv*1e-12; % kg->tonn (1e-3)->GT (1e-9)
  Mv=Mv(:);
  MVOL(:,ip) = Mv; % GT (1e-9)

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
  nf = -1; % <0 - visible off
  ifx=max(strfind(fina,'/'));

  if s_f1==1
    stl=sprintf('%s Depth-intgr. Mass kg/m2, Log2, Tr=%i',...
		fina(ifx+1:end),nTr);
    sub_plot_tracers(lTr,nf,HH,xlim1,...
		     xlim2,ylim1,ylim2,...
		     LON,LAT,stl,nTr,...
		     'c1',2,'c2',7);

    stx{1} = sprintf('Overall Tr Mass, GT = %18.15d',Mtot*1e-12);
    stx{2} = sprintf('Mean DIntgr  conc. = %10.8d kg/m3\n',Ctot);
    axes('Position',[0.27 0.85 0.4 0.05]);
    text(0.05,0.5,stx);
    set(gca,'xtick',[],'ytick',[],'box','off');

    txtb='vol_intgr_tracr008.m';
    bottom_text(txtb,'pwd',1,'position',[0.02 0.05 0.8 0.1]);

    if s_fig>0
      fgnm=sprintf('%sDpthIntgrMass_tr%2.2i_%i%2.2i%2.2i',pthfig,nTr,DV(1:3));
      fprintf('Saving %s\n',fgnm);
      print('-dpng','-r250',fgnm);
    end
  end


  if s_f2==1
    figure(2); clf;
    axes('Position',[0.12 0.1 0.35 0.85]);
% Note different layers - different thickness
% thus need to normalize by dz
    Mvm=Mv./(abs(dzm)); % mass per 1 m of depth 
    plot(Mvm,zzm,'r');
    stl=sprintf('%s Area-intgr. Mass GT/1m depth in v.layers, Tr=%i',...
		fina(ifx+1:end),nTr);
    title(stl,'Interpreter','none');
    set(gca,'tickdir','out',...
	    'xgrid','on','ygrid','on',...
	    'xlim',[0 1.02*max(Mvm)],...
	    'ylim',[-1520 0],...
	    'ytick',[-2000:100:0]);

    %legend('ARCc0.04','ARCc0.08');
    txtb='vol_intgr_tracr008.m';
    %bottom_text(txtb,'pwd',1,'position',[0.02 0.05 0.8 0.05]);
    bottom_text(txtb,'pwd',1);

    if s_fig>0
      fgnm=sprintf('%svrt_LrIntgr_tr%2.2i_%i%2.2i%2.2i',pthfig,nTr,DV(1:3));
      fprintf('Saving %s\n',fgnm);
      print('-dpng','-r250',fgnm);
    end
  end
  

end; % time loop

if s_f3==1
  CLR = [0, 0.4, 0.5;...
	 0, 0.8, 1.0;...
	 0, 0.8, 0.3;...
	 0.6, 0.1, 0; ...
	 1,  0,  0; ...
	 1, 1, 0];
	 
  mxM = 0;
  figure(3); clf;
  axes('Position',[0.12 0.1 0.35 0.85]);
  hold on
  stl = sprintf('LayerIntgr Mass/1m, GT/1m, Tr=%i, ',nTr);
  for iy=1:nyy
    Mv = MVOL(:,iy);
% Note different layers - different thickness
% thus need to normalize by dz
    Mvm=Mv./(abs(dzm)); % mass per 1 m of depth 
    clr = CLR(iy,:);
    plot(Mvm,zzm,'Color',clr);
    yr = YRPLT(iy,1);
    
    stl = sprintf('%s %i',stl,yr);
    mxM = max([mxM, max(Mvm)]);
  end
  
  title(stl,'Interpreter','none');
  set(gca,'tickdir','out',...
	  'xgrid','on','ygrid','on',...
	  'xlim',[0 1.02*mxM],...
	  'ylim',[-1300 0],...
	  'ytick',[-2000:100:0]);

  %legend('ARCc0.04','ARCc0.08');
  txtb='vol_intgr_tracr008.m';
  %bottom_text(txtb,'pwd',1,'position',[0.02 0.05 0.8 0.05]);
  bottom_text(txtb,'pwd',1);

  if s_fig>0
    fgnm=sprintf('%svrt_LrIntgr_tr%2.2i_AllYrs',pthfig,nTr);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end

  
  
end



