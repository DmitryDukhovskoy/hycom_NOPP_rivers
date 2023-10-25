function sub_plot_altmEKE(lTr,X,Y,nf,HH,xlim1,xlim2,...
			 ylim1,ylim2,LON,LAT,...
			 stl,varargin);
if (nf<0),
  nf=abs(nf);
  figure('Visible','off'); clf;
  fprintf('Figure window is turned off\n');
else  
  figure(nf);
  clf
end

% Plot altimetry-derived EKE
%
c1=[];
c2=[];
f_cmp = [];
hps = [];
clbtck=[];
nV = length(varargin);
%keyboard
for k=1:2:nV
  aa=varargin{k};
  aa=lower(aa);
  if strncmp(aa,'c1',2);
    c1=varargin{k+1};
  elseif strncmp(aa,'c2',2);
    c2=varargin{k+1};
  elseif strncmp(aa,'cmp',3);
    f_cmp=varargin{k+1};
  elseif strncmp(aa,'clbpos',6);
    hps=varargin{k+1};
  elseif strncmp(aa,'clbtick',7);
    clbtck=varargin{k+1};
  elseif strncmp(aa,'clblbl',6);
    clblbl=varargin{k+1};
  end;  
  
end

if isempty(hps)
  hps=[0.81 0.12 0.035 0.8];
end


hmsk=HH;
hmsk(HH<0)=nan;

if isempty(hps)
  posc=[0.81 0.12 0.8 0.06];
else
  posc=[hps(1) hps(2) hps(4) hps(3)];
end


switch (f_cmp)
  case(1)
  nint = 200;
  %CMP = colormap_sclr1(nint,c1,c2);
  CMP = colormap_sclr3(nint,c1,c2);
  cmp = CMP.colormap;
  cnt = CMP.intervals;
  nint=length(cmp);
 case(2)
  nint=200;
  cmp = colormap(parula(nint));
  cnt = (c1:(c2-c1)/nint:c2);
 case(3)
  nint = 200;
  CMP = colormap_sclr2(nint,c1,c2);
  cmp = CMP.colormap;
  cnt = CMP.intervals;
  nint=length(cmp);
 case(4)
  nint=200;
  cmp = colormap_red(nint);
  for ik=1:10
    cmp(ik,:)=[1 1 1];
  end
  cmp=smooth_colormap(cmp,10);
  cnt = (c1:(c2-c1)/nint:c2);
 case(5)
  nint=200;
  cmp = colormap(parula(nint));
  for ik=1:20
    cmp(ik,:)=[1 1 1];
  end
  cmp=smooth_colormap(cmp,15);
  cmp=smooth_colormap(cmp,15);
  cnt = (c1:(c2-c1)/nint:c2);
 case(6)
  nint=200;
  cmp = colormap(parula(nint));
  for ik=1:round(nint/10)
    cmp(ik,:)=[0.25 0 0.5];
  end
  for ik=nint:-1:nint-round(nint/10)
    cmp(ik,:)=[1 0.5 0];
  end
  
  cmp=smooth_colormap(cmp,15);
  cmp=smooth_colormap(cmp,15);
  cnt = (c1:(c2-c1)/nint:c2);
end

pcolor(X,Y,lTr); shading flat;
hold on;
%contour(HH,[0 0],'k','Linewidth',1);
%contour(HH,[-6000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',0.5);
contour(HH,[-500 -500],'Color',[0.7 0.7 0.7],'Linewidth',0.5);
%caxis([0 2]);
caxis([c1 c2]);
if f_cmp>0; 
  colormap(cmp); 
else
  cmp = colormap(parula(100));
  cnt = (c1:(c2-c1)/100:c2);
end;
axis('equal');
set(gca,'xlim',...
	[xlim1 xlim2],...
	'ylim',[ylim1 ylim2],...
	'xtick',[],'ytick',[],...
	'Color',[0.5 0.5 0.5]);
clr=[0.9 0.9 0.9];
plot_gridlines(45,10,1,clr,LON,LAT);
title(stl,'Fontsize',12,'Interpreter','none');

hght=[];
lngth=[];
mint=20;
mbx=mint;
fsz=14;
bxc='k';
aend=0;
[az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

%  hb = colorbar;
%set(hb,'position',hps,...
%	 'TickLength',0.04,...
%	 'Fontsize',14);

if ~isempty(clbtck)
  set(hb, 'Ticks',clbtck,...
	 'TickLabels',clblbl);
end
%keyboard

freezeColors;
pcolor(hmsk); shading flat;
colormap([0 0 0]);



%keyboard



return