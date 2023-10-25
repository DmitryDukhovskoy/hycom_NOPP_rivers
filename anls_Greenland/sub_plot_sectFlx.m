function sub_plot_sectFlx(fmatFlx, SEGM, s_fig);
% Plot results from tracer_eddy_sectionsNagai.m

f_plt=1;  % =1 - plot fluxes by sections
          % =2 - plot net cumulative flux

pthfig  = '/Net/tholia/ddmitry/hycom/ARCc0.08/fig_Greenland_river/';

fprintf('Loading %s\n',fmatFlx);
load(fmatFlx);

TM  = TRFLX(1).TM;  % 
CYC = TRFLX(1).CYC;
CYC(TM==0)=2;
TM(TM==0)=732844;
Ndays = TRFLX(1).Nav_averaged_days; % # of days, time averaging

% Fix a bug in a time array
dtm=diff(TM);
I=find(dtm<=0);


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
Nr=6; % 
for ir=4:6
  switch(ir),
   case(1);  % Nordic Seas:
    rgnm='NordicSeas';
    sgm=[1,2,3,4,5];
    fpos=[1,1,1,-1,-1]; % sign conversion, fluxes - given + norm
   case(2);
    rgnm='LabradorSeas';
    sgm=[7,8,21];
    fpos=[-1,1,1]; % sign conversion, fluxes - given + norm
   case(3),
    rgnm='BaffinBay';
    sgm=[6,7,22];
    fpos=[-1,1,1]; % sign conversion, fluxes 
   case(4),   
    rgnm='GreenlandBox';
    sgm=[9,10,11,12];
    fpos=[1,1,1,1]; % sign conversion, fluxes - given + norm
   case(5),
    rgnm='IcelandBox';,
    sgm=[13,14,15,16];
    fpos=[1,1,1,1]; % sign conversion, fluxes - given + norm
   case(6),
    rgnm='LabradorBox';
    sgm=[17,18,19,20];
    fpos=[1,1,1,1]; % sign conversion, fluxes - make sure the fluxes 
  end
  
  cc=cc+1;
  
  figure(ir); clf;
  ns=length(sgm);
%keyboard
  for in=1:ns
    isg=sgm(in); % section #
% Get fluxes:
% surface layer
% Depth intgr
%    UC   = TRFLX(isg).TotalFLux_UC_z2; % <U*C>
    UmCm = TRFLX(isg).RegFlux_UmCm; % <U>*<C>
    UpCp = TRFLX(isg).EddyFlux_UpCp; % <U'*C'>
% Possible bug - need to check sign for iceland
    if isg==13 | isg == 14 | isg == 15 | isg == 16
        UmCm = -UmCm;% <U>*<C>
        UpCp = -UpCp; % <U'*C'>
    end

% Integrated, net  
    if in==1
      UmCmN=UmCm;
      UpCpN=UpCp;
    else
      UmCmN=UmCmN+UmCm;
      UpCpN=UpCpN+UpCp;
    end
    
    
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
      pos=POS(in,:);
      axes('Position',pos);
  %    plot(rt);
      plot(TY,UmCm,'b','LineWidth',2);
      hold on;
      plot(TY,UpCp,'r','LineWidth',2);
      set(gca,'xlim',[xl1 xl2],...
	    'xtick',[1:1/12:TY(end)],...
	    'ylim',[1.05*yl1 1.05*yl2],...
	    'tickdir','out');
      stt=sprintf('EddyFl/TotFlux%%, Section %i',isg);
      title(stt);

      if s_fig==1
	ffg=sprintf('%seddy_flx_sections_HYCOM_rg%2.2i',pthfig,ir);
	fprintf('Saving %s\n',ffg);
	print('-depsc2',ffg);
      end
    end
    
  end;

  Rflx=UmCmN;
  Eflx=UpCpN;
  rt=Eflx./abs(Rflx+Eflx)*100; % Ratio eddy to total flux

  if f_plt==2  % plot cumulative fluxes
    dT=3600*24*Ndays;
    cumRflx=cumsum(Rflx)*dT;  % kg
    cumEflx=cumsum(Eflx)*dT;
    
    axes('Position',[0.08 0.45 0.85 0.44]);
    plot(TY,cumRflx,'linewidth',1.6);
    hold on
    plot(TY,cumEflx,'r-','linewidth',1.6);
    set(gca,'xlim',[7 xl2],...
	    'xtick',[1:2/12:TY(end)],...
	  'xminortick','on',...
	  'tickdir','out',...
    	   'xgrid','on','ygrid','on');

    hg=legend('RegFlx','EddyFlx');
    set(hg,'Position',[0.08 0.18 0.19 0.12]);
    sttl=sprintf('HYCOM: Tracer gain (kg), region %s',rgnm);
    
    title(sttl);
    
    btmtxt='sub_plot_sectFlx.m';
    bottom_text(btmtxt);
      if s_fig==1
	ffg=sprintf('%seddyNagai_cumulative_flx_HYCOM_rg%2.2i',pthfig,ir);
	fprintf('Saving %s\n',ffg);
	print('-depsc2',ffg);
      end
  end
    
    
  
  Rflx=UmCmN;
  Eflx=UpCpN;
  rt=Eflx./abs(Rflx+Eflx)*100; % Ratio eddy to total flux
% Discard values when the total flux is ~0
  rt(abs(Rflx)<1e-1)=nan;
  figure(ir+10); clf;
  axes('Position',[0.08 0.45 0.85 0.44]);
  plot(TY,Rflx,'b','linewidth',1.6); % <U>*<C>
  hold on;
  plot(TY,Eflx,'r','linewidth',1.6);  % <U'>*<C'>
  
  set(gca,'xlim',[7 xl2],...
	  'xminortick','on',...
	  'tickdir','out');
  hg=legend('RegFlx','EddyFlx');
  set(hg,'Position',[0.08 0.18 0.19 0.12]);
  sttl=sprintf('HYCOM: Net flux, region %s',rgnm);
  title(sttl);
  
%  B=[-100:5:100];
%  [N,X]=hist(rt,B);
%  axes('Position',[0.5 0.05 0.35 0.38]);
%  hb=bar(X,N);
%  set(hb,'FaceColor',[0.8 0.8 0.8]);
%  set(gca,'xlim',[X(1)-5 X(end)+5],...%
%	  'xminortick','on',...
%	  'tickdir','out');
  
%  title('Ratio EddyFlx/TotalFlux %%');
  if s_fig>0
    ffg=sprintf('%seddyNagai_netflx_HYCOM_rg%2.2i',pthfig,ir);
    fprintf('Saving %s\n',ffg);
    print('-depsc2',ffg);
    print('-dpng','-r150',ffg);
  end

  
  btmtxt='sub_plot_sectFlx.m';
  bottom_text(btmtxt);
  
% Plot edy cum. flux:
% smooth: 
  dmm=cumEflx;
  dd=15;
  for ik=1:length(dmm);
    i1=ik-dd;
    i2=ik+dd;
    i1=max([i1,1]);
    i2=min([i2,length(dmm)]);
    a=mean(dmm(i1:i2));
    cumEflxF(ik,1)=a;
  end
  
    
%  figure(22);
%  if cc==1
%    axes('Position',[0.08 0.45 0.85 0.44]);
%    hold on;
%  end
%  
%  clr=CLR(ir,:);
%  plot(TY,cumEflxF,'k-','Linewidth',2,'Color',clr);
  
  
% keyboard
 
end    
  
%figure(22);
%title('Cumulative Eddy Flux');
%  if s_fig>0
%    ffg=sprintf('%seddyNagai_netEddyFlx_HYCOM_rg%2.2i',pthfig,ir);
%    fprintf('Saving %s\n',ffg);
%    print('-depsc2',ffg);
%    print('-dpng','-r150',ffg);
%  end





return

