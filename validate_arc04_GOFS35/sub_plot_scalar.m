function sub_plot_scalar(lTr,nf,HH,xlim1,xlim2,...
			 ylim1,ylim2,LON,LAT,...
			 stl,pfld,varargin);

c1=[];
c2=[];
f_cmp = [];
hps = [];
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
  end;  
end

if isempty(hps)
  hps=[0.85 0.12 0.03 0.8];
end


hmsk=HH;
hmsk(HH<0)=nan;

if isempty(hps)
  posc=[0.81 0.12 0.8 0.06];
else
  posc=[hps(1) hps(2) hps(4) hps(3)];
end


% Colormap:
switch(pfld),
 case('salin');
  if isempty(f_cmp), f_cmp=1; end;
%  c1=10;
%  c2=30;
%  c1 = 28.5;
%  c2 = 34.5;
  if isempty(c1)
    c1 = 16.5;
    c2 = 36.5;
  end
%  c1=1.55; % for log scale: S = ~4.7
%  c2=3.55; % log scale S~= 35 
 case('temp')
  if isempty(f_cmp), f_cmp=2; end;

  if isempty(c1)
    c1=-2;
    c2=5.6;
  end
  
 case('AtlT')
%  f_cmp=1;
  if isempty(f_cmp), f_cmp=1; end;
  if isempty(c1)
    c1=-1;
    c2=4;
  end
  posc=[0.88 0.12 0.8 0.06];
 
 case('AtlZ')
  if isempty(f_cmp), f_cmp=1; end;
%  f_cmp = 2;
  if isempty(c1)
    c1=-450;
    c2=-50;
  end
  posc=[0.88 0.12 0.8 0.06];
  
 case('speed')
  if isempty(f_cmp), f_cmp=1; end;
%  f_cmp = 2;
  if isempty(c1)
    c1=0;
    c2=0.35;
  end

  otherwise
  if isempty(f_cmp), f_cmp=1; end;
  if isempty(c1)
    c1=0;
    c2=0.35;
  end
  
end

switch (f_cmp)
  case(1)
  nint = 200;
  %CMP = colormap_sclr1(nint,c1,c2);
  CMP = colormap_sclr3(nint,c1,c2);
  cmp = CMP.colormap;
  cnt = CMP.intervals;
  nint=length(cmp);

%  for ik=1:round(nint*0.05);
%   cmp(ik,:)=[0.95 0.95 1];
%  end
%  cmp = smooth_colormap(cmp,9);
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
  nint = 200;
  CMP = colormap_sclr4(nint,c1,c2);
  cmp = CMP.colormap;
  cnt = CMP.intervals;
  nint=length(cmp);
end


if (nf<0),
  nf=abs(nf);
  figure('Visible','off'); clf;
  fprintf('Figure window is turned off\n');
else  
  figure(nf);
  clf
end
axes('Position',[0.05 0.2 0.75 0.7]);
hold on;
pcolor(hmsk); shading flat;
colormap([0.6 0.6 0.6; 1 1 1]);
freezeColors;

%keyboard


pcolor(lTr); shading flat;
%contour(HH,[0 0],'k','Linewidth',1);
%contour(HH,[-6000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',0.5);
%contour(HH,[-500 -500],'Color',[0.7 0.7 0.7],'Linewidth',0.5);
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
	'Color',[0.8 0.8 0.8]);
clr=[0.9 0.9 0.9];
%plot_gridlines(45,10,1,clr,LON,LAT);
title(stl,'Fontsize',12,'Interpreter','none');

%if f_cmp>0
%  hght=[];
%  lngth=[];
%  mint=20;
%  mbx=mint;
%  fsz=14;
%  bxc='k';
%  aend=0;
%  [az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);
%else
hb = colorbar;
set(hb,'position',hps,...
     	 'TickLength',0.02,...
	      'Fontsize',14);
%end




return
