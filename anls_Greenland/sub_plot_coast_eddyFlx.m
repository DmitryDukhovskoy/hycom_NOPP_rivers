function sub_plot_coast_eddyFlx(fmatFlx, SEGM, s_fig, pthfig);
% Plot results from tracer_eddy_sectionsNagai.m
% Plot eddy fluxes across coastal sections
% around Greenland

f_plt=2;  % =1 - plot fluxes by sections
          % =2 - plot net cumulative flux

%pthfig  = '/Net/tholia/ddmitry/hycom/ARCc0.08/fig_Greenland_river/';

fprintf('Loading %s\n',fmatFlx);
load(fmatFlx);

TM  = TRFLX(1).TM;  % 
CYC = TRFLX(1).CYC;
CYC(TM==0)=2;
TM(TM==0)=732844;
Ndays = TRFLX(1).Nav_averaged_days; % # of days, time averaging

% Fix a bug in a time array
%dtm=diff(TM);
%I=find(dtm<=0);


TD = sub_TM(TM,CYC);
TY = TD./365.25;
xl1=TY(1);
xl2=TY(end);

%TY=TY(2:end);
%TD=TD(2:end);

% Check noraml sign convention:
% in SEGM:
% vrs=3
%fmat = sprintf('%sSubarctic_sections_v%i',pthmat,vrs); % closed Lab and Baffin

xp1=0.08;
yp1=0.7;
xax=0.4;
yax=0.23;
dyax=0.09;

POS=[xp1,      yp1, xax, yax; ...
     0.97-xax, yp1, xax, yax;...
     xp1,      yp1-(dyax+yax), xax, yax;...
     0.97-xax, yp1-(dyax+yax), xax, yax;...
     xp1,      0.05, xax, yax;...
     0.97-xax, 0.05, xax, yax];

CLR=[1,0,0;0,0,0;1,1,0;0,0,1;0,1,0;0,1,1];

cc=0;
Nr=2; % 2 contours around Greenland
nsgm=length(SEGM)/Nr; % # segments
ifg=0;
for ir=1:Nr  % Contours
  ifg=ifg+1;
  in=0;
  figure(ifg); clf;
  for isg=1:nsgm   % segments
    nisg=(ir-1)*nsgm+isg;
    rgnm=sprintf('Contour %i, segm %i',ir,isg);
%    subplot(4,3,isg);  

% Get fluxes:
% surface layer
% Depth intgr
%    UC   = TRFLX(isg).TotalFLux_UC_z2; % <U*C>
    UmCm = TRFLX(nisg).RegFlux_UmCm; % <U>*<C>
    UpCp = TRFLX(nisg).EddyFlux_UpCp; % <U'*C'>

    
    mn1=min(UmCm);
    mn2=min(UpCp);
    mx1=max(UmCm);
    mx2=max(UpCp);
    yl1=min([mn1,mn2,0]);
    yl2=max([mx1,mx2,0]);
    
%    rt=UpCp./abs(UmCm+UpCp)*100; % Ratio eddy to total flux
%    yl1=min(rt);
%    yl2=max(rt);

    if f_plt==1
      in=in+1;
      if (in>6);
	in=1;
	ifg=ifg+1;
	figure(ifg); clf;
      end
      
      pos=POS(in,:);
      axes('Position',pos);
  %    plot(rt);
      plot(TY,UmCm,'b','LineWidth',2);
      hold on;
      plot(TY,UpCp,'r','LineWidth',2);
      set(gca,'xlim',[xl1 xl2],...
	    'xtick',[1:1/2:TY(end)],...
	    'xminortick','on',...
	    'ylim',[1.05*yl1 1.05*yl2],...
	    'tickdir','out');
      stt=sprintf('<U><C> & <U"*C"., Contour %i, Segm %i',ir,isg);
      title(stt);

      if s_fig==1
	ffg=sprintf('%seddy_flx_sections_HYCOM_rg%2.2i',pthfig,ir);
	fprintf('Saving %s\n',ffg);
	print('-depsc2',ffg);
      end
    end

% Time-integrated tracer flux: kg 
    dT=3600*24*Ndays;
    cumRflx=cumsum(UmCm)*dT;  % kg
    cumEflx=cumsum(UpCp)*dT;

    CUMFLX(ir).Rflx(isg)=cumRflx(end);
    CUMFLX(ir).Eflx(isg)=cumEflx(end);
    
    if f_plt==2  % plot cumulative fluxes
      in=in+1;
      if (in>6);
	in=1;
	ifg=ifg+1;
	figure(ifg); clf;
      end
      
      pos=POS(in,:);
      axes('Position',pos);

%      axes('Position',[0.08 0.45 0.85 0.44]);
      plot(TY,cumRflx,'linewidth',1.6);
      hold on
      plot(TY,cumEflx,'r-','linewidth',1.6);
      set(gca,'xlim',[round(xl1) xl2],...
	      'xtick',[1:1/2:TY(end)],...
	    'xminortick','on',...
	    'tickdir','out',...
	    'yminortick','on',...
	     'xgrid','on','ygrid','on');

      sttl=sprintf('Tracer cumflx (kg), Cntr %i, sgm %i',ir,isg);

      title(sttl);

      if in==6 | isg==nsgm
	hg=legend('RegFlx','EddyFlx');
	set(hg,'Position',[0.02 0.08 0.18 0.1]);

	btmtxt='hycom_NOPP_rivers/anls_Greenland/sub_plot_sectFlx.m';
	bottom_text(btmtxt);
	if s_fig==1
	  ffg=sprintf('%seddyNagai_cumFlx_cntr%2.2i-sgm%2.2i',pthfig,ir,isg);
	  fprintf('Saving %s\n',ffg);
	  print('-depsc2',ffg);
	end
      end	
      
      
    end  % f_plt==2
  end  % segment  
end    % contour 1 or 2 -
    
%  keyboard
% Plot bar diagram for all fluxes:
for ir=1:Nr
  for isg=1:nsgm
    BB(ir,isg)=CUMFLX(ir).Eflx(isg);
    RR(ir,isg)=CUMFLX(ir).Rflx(isg);
  end
end

figure(11); clf;
cmp=[1,0.3,0; 0.6,0,0];
ylm1=min(min(BB));
ylm2=max(max(BB));

axes('Position',[0.08 0.45 0.65 0.45]);
hb=bar(BB',1,'grouped');
colormap(cmp);
%set(hb,'edgecolor','none');
set(gca,'tickdir','out',...
	'fontsize',16,...
	'ylim',[1.02*ylm1 1.1e13],...
	'xlim',[0.5 nsgm+0.5],...
	'yminortick','on');
hl=legend('Cntr1','Cntr2');
set(hl,'Position',[0.523 0.51 0.136 0.088]);
title('Time-integrated Eddy Tracer Flux, kg');
xlabel('Segments');
axes('Position',[0.08 0.1 0.1 0.1]); % fix for glitch with legend, disappears with btm_text
set(gca,'visible','off');
%freezeColors;
btmtxt='ddmitry@mars: hycom_NOPP_rivers/anls_Greenland/sub_plot_sectFlx.m';
bottom_text(btmtxt);
if s_fig==1
  ffg=sprintf('%scumEddyFlx_bardgrm_coastGreen',pthfig);
  fprintf('Saving %s\n',ffg);
%  print('-depsc2',ffg);
  print('-dpng','-r250',ffg);
end

figure(12); clf;
cmp=[0,0.4,1; 0,0,0.6];
ylm1=min(min(RR));
ylm2=max(max(RR));

%axes('Position',[0.08 0.5 0.6 0.35]);
axes('Position',[0.08 0.45 0.65 0.45]);
hb=bar(RR',1,'grouped');
colormap(cmp);
%set(hb,'edgecolor','none');
set(gca,'tickdir','out',...
	'fontsize',16,...
	'ylim',[1.02*ylm1 1.5e14],...
	'xlim',[0.5 nsgm+0.5],...
	'yminortick','on');
title('Time-integrated Mean Tracer Flux, kg');
xlabel('Segments');
hl=legend('Cntr1','Cntr2');
%set(hl,'Position',[0.523 0.11 0.136 0.088]);
set(hl,'Position',[0.523 0.51 0.136 0.088]);

axes('Position',[0.08 0.1 0.1 0.1]); % fix for glitch with legend, disappears with btm_text
set(gca,'visible','off');

btmtxt='ddmitry@mars: hycom_NOPP_rivers/anls_Greenland/sub_plot_sectFlx.m';
bottom_text(btmtxt);
if s_fig==1
%  ffg=sprintf('%seddyNagai_Bars_cumFlx_coastGreen',pthfig,ir,isg);
  ffg=sprintf('%scumMeanFlx_bardgrm_coastGreen',pthfig);
  fprintf('Saving %s\n',ffg);
%  print('-depsc2',ffg);
  print('-dpng','-r250',ffg);
end



return

