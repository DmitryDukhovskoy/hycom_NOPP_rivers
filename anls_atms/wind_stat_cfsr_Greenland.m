% Calculate monthly climatology
% cfsr and cfsv2 winds for Greenland coast
% Note CFSv2 - monthly means, need to get hourly
% speed - mean scalar,
% direction - mean vector, i.e. predominant direction
% for different regions
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08/;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/;
startup

clear
close all


ys=1998;  % 1st year
ye=2010;  % last year (year of Oct-Dec segment)
%% To Plot wind roses for select years and regions:
% make f_roses2=1 and select years when call subfunction

% -------------
% Flags
% -------------
f_getd  = 0;  % = 1 - get data; = 0 - load extracted data
f_roses = 1;  % =1 - plot PDF, wind roses
f_roses2= 0;  % wind roses for selected years
f_histV = 1;  % plot hist of V component
s_fig   = 0;  % figure print

ndate0=datenum(ys,1,1,0,0,0);
ndate1=datenum(ye,12,31,0,0,0);     % last date

dv1=datevec(ndate0);
dv2=datevec(ndate1);

fprintf('Time range: %i/%i/%i - %i/%i/%i\n',dv1(1:3),dv2(1:3));

%pthdat = '/Net/data/ccmp/v02.0/';
pthdat1 = '/Net/kronos/ddmitry/ncep_cfsr/';
pthdat2 = '/Net/kronos/ddmitry/ncep_cfsv2/';
pthdt  = '/Net/tholia/ddmitry/ovwst/data_mat/';
pth8   = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig = '/Net/tholia/ddmitry/ovwst/fig_cycl/';
pthtopo= '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

%fmat = sprintf('%swrose_Greenland_ccmp020rss_%i_%i.mat',pthdt,ys,ye);
%fmat = sprintf('%swrose_Greenland_ccmp020rss_%i_%iv1.mat',pthdt,ys,ye);
fmat = sprintf('%swrose_Greenland_cfsr_cfsv2_%i_%i.mat',pthdt,ys,ye);


btx = 'wind_stat_cfsr_Greenland.m';

dnmb=datenum(2005,1,1);
d0vct=datevec(dnmb);
year=d0vct(1);
mnth=d0vct(2);
mday=d0vct(3);
fp=sprintf('%scfsr-sea_%4.4i_01hr_uv-10m.nc4',pthdat1,year);

%x1=210;
x1=270;  % Westernmost coord
x2=360;  % Easternmost coord.
y1=50;
y2=78.375;


% CFSR Ocean Wind vector data:
alat=nc_varget(fp,'Latitude');
elon=nc_varget(fp,'Longitude');
n=length(elon);
m=length(alat);


dx=abs(elon(2)-elon(1));
if x1<0
  x1a=x1+360;
else
  x1a=x1;
end

dst=sqrt((elon-x1a).^2);
i1=min(find(dst==min(dst)));
dst=sqrt((alat-y1).^2);
j1=min(find(dst==min(dst)));
dst=sqrt((elon-x2).^2);
i2=min(find(dst==min(dst)));
dst=sqrt((alat-y2).^2);
j2=min(find(dst==min(dst)));

[X,Y] = meshgrid(elon(i1:i2),alat(j1:j2));
subalat=alat(j1:j2);
if i2<i1
  e1=elon(i1:end);
  e2=elon(1:i2);
  subelon=[e1;e2];
else
 subelon=elon(i1:i2);
end


WRS = head_roses_Greenland(subalat,subelon);


% HYCOM grid/topo for plotting
ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);




itot=0;
crec=0;
nday=0;
if f_getd>0
  
  WCLIM=struct;
  
  for year=ys:ye
    for iday=1:366
      dnum=datenum(year,1,1)+iday-1;
      DV=datevec(dnum);
      YR=DV(1);
      m1=DV(2);
      dd=DV(3);
      mo=m1;
      if DV(1)>year, % 365-day year, jumped to next year
	continue;
      end
      

  %  year=ys;
      ndate0=datenum(year,m1,dd,0,0,0);
      dnmb=ndate0;
      ndate1=ndate0+1;     % next day
      tstep=0.25;     % step 6 hrs, 1-hr fields, 24 fields per day
      nrc1d=24;  % step 6 hrs, 1-hr fields, 24 fields per day
      ndate_crnt=ndate0-tstep;
      dj1=datenum(year,1,1,0,0,0);      % Time in CCMP is hours since 1987-1-1
      it=0;
      jcntr=0;
% 1-hr data:
      while ndate_crnt<ndate1-tstep
	jcntr=jcntr+1;
	itot=itot+1;
	ndate_crnt=ndate_crnt+tstep;
	date_s=datestr(ndate_crnt,0);
	dvv=datevec(ndate_crnt);
	hr=dvv(4);
	d0vct=dvv;

        year=d0vct(1);
        mnth=d0vct(2);
        mday=d0vct(3);
	if year<2011
          fp = sprintf('%scfsr-sea_%4.4i_01hr_uv-10m.nc4',pthdat1,year);
	else
          fp = sprintf('%scfsv2-sec2_%i_mon_uv-10m.nc',pthdat2,year);
        end	  

	if mod(jcntr-1,4)==0
	  fprintf('Day: %s, %s\n',date_s,fp);
	  disp(['Reading day: ',date_s]);
	  disp(fp)
	end

	it=nrc1d*(ndate_crnt - floor(dj1))+1;

      % ********************************
      % Check dates:
      % ********************************
	tm=nc_varget(fp,'Date',[it-1],[1]);
        yr0=floor(tm/10000);
	mn0=floor((tm-yr0*10000)/100);
	dd0=floor((tm-yr0*10000-mn0*100));
	hr0=(tm-floor(tm))*24;
	dchck=datenum(yr0,mn0,dd0,hr0,0,0);
	if dchck~=ndate_crnt,
	  datestr(dchck,0)
	  datestr(ndate_crnt,0)
	  error('*** ERR:  Check dates in CFSR/CFSv2 data and ndate_crnt  ***');
	end;
      % ********************************

	djj=j2-j1+1;
	if x1a>x2
	  dii=n-i1+1;
	  if dii<=0, error('Check indices along long. i1'); end;
	  A1=nc_varget(fp,'uwnd',[it-1 j1-1 i1-1],[1 djj dii]);
	  A2=nc_varget(fp,'uwnd',[it-1 j1-1 0],[1 djj i2]);
	  U=[A1,A2];

	  A1=nc_varget(fp,'vwnd',[it-1 j1-1 i1-1],[1 djj dii]);
	  A2=nc_varget(fp,'vwnd',[it-1 j1-1 0],[1 djj i2]);
	  V=[A1,A2];

	else

	  dii=i2-i1+1;
	  if dii<=0, error('Check indices along long. i1'); end;
	  U=nc_varget(fp,'wndewd',[it-1 j1-1 i1-1],[1 djj dii]);
	  V=nc_varget(fp,'wndnwd',[it-1 j1-1 i1-1],[1 djj dii]);

	end

	S=sqrt(U.^2+V.^2);
	I=find(~isnan(S));
	if isempty(I); 
	  fprintf('\n\n !!! Record %s is corrupted - no data, skipping ...\n\n',...
		  datestr(dchck));
	  continue; 
	end; 
	
        if jcntr==1
	  S6=S*0;
	  U6=U*0;
	  V6=V*0;
	end;
	
        S6=S6+S;
        U6=U6+U;
	V6=V6+V;
	
	if itot==1
	  for ik=1:12
	    WCLIM(ik).Title='Monthly mean CFSR/CVSF2';
	    WCLIM(ik).code='wind_stat_cfsr_Greenland.m';
	    WCLIM(ik).CC=0;
	    WCLIM(ik).Wspeed=zeros(size(S6));
	    WCLIM(ik).U=zeros(size(S6));
	    WCLIM(ik).V=zeros(size(S6));
	    WCLIM(ik).Day1=[];
          end;
	end

% Wind direction and max wind stat for boxes:
        crec=crec+1;
        WRS(1).TM(crec)=ndate_crnt;
        WRS = sub_winddir_stat(WRS,U,V,S,crec);
      
      end;  % while - 6hr 
% Daily mean
      S6=S6./jcntr;
      U6=U6./jcntr;
      V6=V6./jcntr;

      cc=WCLIM(mo).CC;
      Sm=WCLIM(mo).Wspeed;
      Um=WCLIM(mo).U;
      Vm=WCLIM(mo).V;
      
      cc=cc+1;
      if isempty(WCLIM(mo).Day1),
	WCLIM(mo).Day1=dnum;
      end
      
      
% Update monthly means:
      a1=(cc-1)/cc;
      a2=1/cc;
      Sm=a1.*Sm+a2.*S6;
      Um=a1.*Um+a2.*U6;
      Vm=a1.*Vm+a2.*V6;

      WCLIM(mo).CC=cc;
      WCLIM(mo).Wspeed=Sm;
      WCLIM(mo).U=Um;
      WCLIM(mo).V=Vm;
      WCLIM(mo).Day2=dnum;
 
    end;  % for iday 
%      keyboard
    

    fprintf('----> Saving %s\n',fmat);
    save(fmat,'WCLIM','WRS');

  end; % for year

  fprintf('Last Step:   Saving %s\n',fmat);
  save(fmat,'WCLIM','WRS');

else
  fprintf('Loading %s\n',fmat);
  load(fmat);
end


% ====================================
% Plot Wind PDF, and roses in Winter/Summer:
%IRG=1; 

% Plot regions:
figure(1); clf;
contour(HH,[0 0],'Color',[0.6 0.6 0.6]);
hold on;
contour(HH,[-500 -500],'Color',[0.8 0.8 0.8]);

nrg=length(WRS);
for irg=1:nrg
  xv=WRS(irg).Box_VX;
  yv=WRS(irg).Box_VY;
  [IC,JC] = sub_find_indxHYCOM(LON,LAT,xv,yv);
  plot(IC,JC,'r.-');
  text(IC(1),JC(1),sprintf('%i',irg),'Fontsize',12);
end
axis('equal');
bottom_text(btx,'pwd',1);


IRG=[1:nrg];
fgn=10;
ttl='ccmp020rss';
if f_roses>0
  sub_plot_wind_stat(WRS,IRG,fgn,pthfig,s_fig,ttl,btx);
%  bottom_text(btx,'pwd',1);
end

% Plot wind roses for select years and regions:
IRG=[9];
YR = 2014;
fgn=30;
ttl='Select YRS, ccmp020rss';
if f_roses2>0
  sub_wroses_YRRG(WRS,IRG,YR,fgn,pthfig,ttl,btx);
%  bottom_text(btx,'pwd',1);
end
if s_fig==1
  figure(16);
  fgnm1=sprintf('%swroses_cfsr_Greenland_slct%i_rg6',pthfig,YR);
  print('-dpng','-r300',fgnm1);

  figure(17);
  fgnm2=sprintf('%swroses_cfsr_Greenland_slct%i_rg7',pthfig,YR);
  print('-dpng','-r300',fgnm2);
end  

% Plot hist of V component 
% to check wind roses on western Greenland shelf
IRG = [8:11];
fnm=40;
if f_histV==1
  TM = WRS(1).TM;
  DV = datevec(TM);
  for ik=1:length(IRG)
    ir= IRG(ik);
    V = WRS(ik).V;
    I = find(DV(:,2)>=10 | DV(:,2)<4);
    dmm = V(I);
    
    figure(ik+fnm); clf;
    axes('position',[0.1 0.55 0.83 0.35]);
    hist(dmm,50,'Facecolor',[0.5 0.5 0.5],50);
    ttl=sprintf('Winter ONDJFM, Vcomp m/s, pnt=%i',ir);
    title(ttl);

    I = find(DV(:,2)>3 & DV(:,2)<10);
    dmm = V(I);
    axes('position',[0.1 0.08 0.83 0.35]);
    hist(dmm,50,'Facecolor',[0.5 0.5 0.5]);
    ttl=sprintf('Summer AMJJAS, Vcomp m/s, pnt=%i',ir);
    title(ttl);

    bottom_text(btx,'pwd',1);
    
  end
end
