% Plot fluxes in POP boxes around Greenland
% calculated in greenl_fluxes_POPboxes.m
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=110;
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';

YR=2005;

fmatout=sprintf('%shycom008_%3.3i_Greenl_flx_POPboxes_%4.4i.mat',...
                 pthmat,expt,YR);
btx='plot_greenl_fluxes_POPboxes.m';

fprintf('Loading %s\n',fmatout);
load(fmatout);

% Check Gr flux:
nbx=6;

for ibx=1:nbx
  VflxG(ibx,:)=VHFLX(ibx).S(1).Vflx_Gr;
end;

absVtot = 300; % mean overall abs flux through Gr contour, Sv = sum (abs(U)*dx*dz)
VGtot=GVHFLX.Vflx*1e-6;  % Sv
VGsum=sum(VflxG)*1e-6;
mnV=mean(VGtot);
errV=abs(mnV/absVtot); % absolute relative error

figure(1); clf;
axes('Position',[0.08 0.42 0.88 0.4]);
hold on;
plot(VGtot,'-','Color',[0 0 0],'Linewidth',4);
plot(VGsum,'--','Color',[0 0.8 1]);
title(sprintf('Greenland contour vol flux: total and sum(segments), Sv, %i',YR));
lg=legend('Tot contour','Sum segments');

t1=sprintf('Mean GrCntr flux = %6.2f+/-%6.2f Sv',mnV,std(VGtot));
t2=sprintf('Absolute rel. error=%8.6f',errV);
text(10,1.3,t1,'Fontsize',12);
text(10,1.,t2,'Fontsize',12);

nn=length(VGsum);
set(gca,'tickdir','out',...
        'xlim',[1 nn],...
        'ylim',[-2.2 2.2],...
        'ytick',[-10:0.5:10],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);
xlabel('Days');
ylabel('Sv');

bottom_text(btx,'pwd',1,'Position',[0.08 0.3 0.4 0.05]);


% Plot by boxes:
for ibx=1:nbx
  ibx0=ibx+1;
  if ibx==nbx; ibx0=1; end;
  vs1(ibx,:)=VHFLX(ibx).S(1).VolFlxGrSh_m3s*1e-6;
  vs2(ibx,:)=-VHFLX(ibx0).S(1).VolFlxGrSh_m3s*1e-6;
  vg(ibx,:)=VHFLX(ibx).S(1).Vflx_Gr*1e-6;

  h1s1(ibx,:)=VHFLX(ibx).S(1).HFlxGrSh_T1_W; 
  h1s2(ibx,:)=VHFLX(ibx0).S(1).HFlxGrSh_T1_W; 
  h1g(ibx,:)=VHFLX(ibx).S(1).Hflx1_Gr;

  h2s1(ibx,:)=VHFLX(ibx).S(1).HFlxGrSh_T2_W; 
  h2s2(ibx,:)=VHFLX(ibx0).S(1).HFlxGrSh_T2_W; 
  h2g(ibx,:)=VHFLX(ibx).S(1).Hflx2_Gr;
end
s3nm='GrContour';
s4nm='-(v1+v2)';

CLR=[0.8 0.2 0; ...  % S1
    0.3 0.9 0; ...  %S2
    0 0.5 0.8;...     %Gr
    0.7 0.7 0.7];

for ibx=1:nbx
  ibx0=ibx+1;
  if ibx==nbx; ibx0=1; end;
  s1nm=BOX(ibx).S(1).Name;
  s2nm=BOX(ibx0).S(1).Name;

  figure(ibx+1); clf;
  axes('Position',[0.09 0.68 0.85 0.25]);
  hold on;

  nm=BOX(ibx).Name; 
  v1=vs1(ibx,:);
  v2=vs2(ibx,:);
  v3=vg(ibx,:);
  v0=-(v1+v2);  % what Gr flux should be
  plot(v1,'Linewidth',2,'Color',CLR(1,:));
  plot(v2,'Linewidth',2,'Color',CLR(2,:));
  plot(v3,'Linewidth',2,'Color',CLR(3,:));
  plot(v0,'Linewidth',2,'Color',CLR(4,:));

  ts=sprintf('0.08 HYCOM, expt-110, Box %s %i Vol Flux, Sv',nm,YR);
  title(ts,'Interpreter','none');

  v1m=mean(v1);
  v2m=mean(v2);
  v3m=mean(v3);

  dF=v1m+v2m+v3m;
  absErr=abs(dF)/(abs(v1m)+abs(v2m)+abs(v3m));

  yl1=-1.05*max(abs([v1,v2,v3,v0]));
  yl2=-yl1;

  txt{1}=sprintf('%s %6.2f Sv',s1nm,v1m);
  txt{2}=sprintf('%s %6.2f Sv',s2nm,v2m);
  txt{3}=sprintf('%s %6.2f Sv',s3nm,v3m);
  txt{4}=sprintf('Tot flux: %6.2f Sv',dF);

  text(5,0.9*yl2,txt);
  if (yl2-yl1) <4, 
    dyy=0.5; 
  else
    dyy=1;
  end

  nn=length(v1);
		set(gca,'tickdir','out',...
										'xlim',[1 nn],...
										'ylim',[yl1 yl2],...
										'ytick',[-10:dyy:10],...
										'xgrid','on',...
										'ygrid','on',...
										'Fontsize',12);
%		xlabel('Days');
		ylabel('Sv');

%  lgd=legend(s1nm,s2nm,s3nm,s4nm);
%
%  set(lgd,'Position',[0.86 0.57 0.13 0.12],...
%          'Box','off',...
%          'Fontsize',10);

%
% Heat fluxes, T=-1.8
  axes('Position',[0.09 0.38 0.85 0.25]);
  hold on
  v1=h1s1(ibx,:)*1e-12;
  v2=h1s2(ibx,:)*1e-12;
  v3=h1g(ibx,:)*1e-12;
%  v0=-(v1+v2);  % 
  plot(v1,'Linewidth',2,'Color',CLR(1,:));
  plot(v2,'Linewidth',2,'Color',CLR(2,:));
  plot(v3,'Linewidth',2,'Color',CLR(3,:));
%  plot(v0,'Linewidth',2,'Color',[0.7 0.7 0.7]);

  ts=sprintf('HeatFlx, TW, Tref=-1.8C');
  title(ts,'Interpreter','none');

  yl1=-1.05*max(abs([v1,v2,v3]));
  yl2=-yl1;

  set(gca,'tickdir','out',...
          'xlim',[1 nn],...
          'ylim',[yl1 yl2],...
          'xgrid','on',...
          'ygrid','on',...
          'Fontsize',12);
%  xlabel('Days');
  ylabel('TW');

% Heat fluxes
  axes('Position',[0.09 0.08 0.85 0.25]);
  hold on
  v1=h2s1(ibx,:)*1e-12;
  v2=h2s2(ibx,:)*1e-12;
  v3=h2g(ibx,:)*1e-12;
%  v0=-(v1+v2);  % 
  plot(v1,'Linewidth',2,'Color',CLR(1,:));
  plot(v2,'Linewidth',2,'Color',CLR(2,:));
  plot(v3,'Linewidth',2,'Color',CLR(3,:));
%  plot(v0,'Linewidth',2,'Color',[0.7 0.7 0.7]);

  ts=sprintf('HeatFlx, TW, Tref=0C');
  title(ts,'Interpreter','none');

  yl1=-1.05*max(abs([v1,v2,v3]));
  yl2=-yl1;

  set(gca,'tickdir','out',...
          'xlim',[1 nn],...
          'ylim',[yl1 yl2],...
          'xgrid','on',...
          'ygrid','on',...
          'Fontsize',12);
%  xlabel('Days');
  ylabel('TW');


  axes('Position',[0.85 0.85 0.12 0.12]);
  hold on;
  x1=0.1; 
  y1=0.1;
  dy=0.1;
  dx=0.1;
  plot([x1 x1+dx],[y1 y1],'-','Linewidth',2,'Color',CLR(4,:));
  text(x1+1.2*dx,y1,s4nm,'Fontsize',11);
  y1=y1+dy;
  plot([x1 x1+dx],[y1 y1],'-','Linewidth',2,'Color',CLR(3,:));
  text(x1+1.2*dx,y1,s3nm,'Fontsize',11);
  y1=y1+dy;
  plot([x1 x1+dx],[y1 y1],'-','Linewidth',2,'Color',CLR(2,:));
  text(x1+1.2*dx,y1,s2nm,'Fontsize',11);
  y1=y1+dy;
  plot([x1 x1+dx],[y1 y1],'-','Linewidth',2,'Color',CLR(1,:));
  text(x1+1.2*dx,y1,s1nm,'Fontsize',11);
 
  set(gca,'xlim',[0 0.4],...
          'ylim',[0 0.6],...
          'Fontsize',11,...
          'visible','off');




		bottom_text(btx,'pwd',1);


end


