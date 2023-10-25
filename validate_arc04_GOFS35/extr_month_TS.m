% Extract T or S fields from ARCc0.04 
% For 1 layer at a time
% if need sev layers - modify the code to take into account dH(t)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;

regn = 'ARCc0.04';
%expt = 011;  
%expt = 012;  
expt = 022;

plr   = 1;  % layer to plot
%plr   = 24;  % la
pfld  = 'temp';
%pfld  = 'salin';


rg = 9806;
% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

%pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/frames_TS/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/t_mnth/';
pthmat = pthout;

%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);


ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);



% Average fields:
cnc=0;
ip1=1;
for YR=2016:2016
 
  TMO = struct;
  TMO.Info = 'Monthly surface T 0.04 HYCOM-CICEv5 expt 022'; 
  for im=12:12
    fmatout = sprintf('%s%s_monthly_lr%3.3i_%i%2.2i.mat',pthmat,pfld,plr,YR,im);
    fprintf('Output file %s\n',fmatout);

    d1=datenum(YR,im,1);
    d2=d1+32;
    dv=datevec(d2);
    dmm=datenum(dv(1),dv(2),1);
    d2=dmm-1;
    id1=1;
    id2=d2-d1+1;

    Tsum = [];
    icc  = 0;
    for mday=id1:id2
      dnmb = datenum(YR,im,mday);
      iday = dnmb-datenum(YR,1,1)+1;
      yr = YR;

						switch(expt)
							case(12)
								pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i/',expt,yr);  
								fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
								finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
							case(22);
								pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i/',yr);
								fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
								finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
						end
  
						if ~exist(fina,'file');
								fprintf('Not found: %s\n\n',fina);
								continue;
						end
						
						DV=datevec(dnmb);

						fprintf('%s, V.Lr=%i, %4.4i_%2.2i_%2.2i: %s\n',...
							pfld,plr,DV(1:3),fina);
 
% Layer thickness:
						[F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
						F=squeeze(F);
						F=F./rg;
						F(F>1e10)=0;
						F(F<1e-2)=0;
						dH=squeeze(F); 

						tic;
						[F,n,m,l] = read_hycom(fina,finb,pfld,'r_layer',plr);

						F(F>1e6)=nan;

				%  F = squeeze(F);
				%  lTr = log(F);
						TRC = squeeze(F);      

% Plot sea ice conc:
      pthice = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%i_cice/',expt,DV(1));
      flice = sprintf('%s022_cice_inst.%i-%2.2i-%2.2i-00000.nc',pthice, DV(1:3));
      Ci = squeeze(nc_varget(flice,'aice'));

      f_chck=0;
      if f_chck==1
        figure(10); clf;
        pcolor(TRC); shading flat;
        colorbar
        hold on;
        caxis([-10 -2]);
        contour(Ci,[0.2 0.2],'r');

      end


      icc=icc+1;
      if isempty(Tsum)
        Tsum = TRC;
      else
        Tsum=Tsum+TRC;
      end

      Tmin=min(min(Tsum))/icc;
      Tmax=max(max(Tsum))/icc;  

      fprintf('%s  min/max = %8.4f %8.4f C\n',pfld,Tmin,Tmax);
      fprintf('1 day processed %6.4f min\n',toc/60);
    end

    Tmp = Tsum/icc;

    TMO.Temp = Tmp;
  
    if icc>0
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'TMO');
    end

keyboard

  end
end;  % day loop





