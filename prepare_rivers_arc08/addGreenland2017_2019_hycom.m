% Preprare monthly river runoff data
% for 2017-2019
%    No need to create identical river monthly fields for 2017-2019
% simply copie existing 2016 to the following years and change year in the file name



% Using UCAR/NCAR CGD's  data set - extended through 2016 (data are till 2014)
% Dai and Trenberth data set
% http://www.cgd.ucar.edu/cas/catalog/surface/dai-runoff/index.html
% NCAR/UCAR data have been processed (gaps filled, 
% rivers with missing too many data points deleted, etc)
% in fix_NCARrivers.m
% Added N. Pacific rivers: Yucon, Kuskokiwm, Copper, etc. 
%
% Use Greenland runoff based on updated Bamber data (paper 2018)
% using last year 2016 for 2017-2019
%

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup

close all
clear

regn='ARCc0.08';
ntopo=11;      % HYCOM topo version
%YearGr=2009;  % Gr. river runoff river

%fprintf('River Runoff Year: %i\n',YearGr);

hg=2^100; 

PTH.data    = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.dataout = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/riversNCAR/';
PTH.topo    = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%PTH.river   = '/Net/ocean/ddmitry/arctic_AOregimes/data/Greenland_rivers/';
PTH.river   = '/Net/ocean/ddmitry/arctic_AOregimes/data/GreenlandRunoffv3/';

% Read Bamber's data: fluxes = km3/mo
% Total FW flux = D+Rg+Rt (no CAA, Svlabard here);
% Shown in Bamber's Fig.3, 2018
fnm = sprintf('%sFWF17.v3.nc',PTH.river);
tmm = double(nc_varget(fnm,'TIME'));
Xgr = nc_varget(fnm,'lon');
Ygr = nc_varget(fnm,'lat');
Rt  = double(nc_varget(fnm,'runoff_tundra')); % tundra runoff, km3/mo
Rg  = double(nc_varget(fnm,'runoff_ice')); % GrIS runoff - meltwater
Ds  = nc_varget(fnm,'solid_ice'); % solid ice
LGr = nc_varget(fnm,'LSMGr'); % Greenland mask
TM  = datenum(1958,1,1)+tmm; 
DV  = datevec(TM);

for iyr=2017:2019
  YearGr=2016;
  year = iyr;
  fprintf('River Runoff Year: %i, Topography: %i\n',YearGr,ntopo);

  mday=[31;28;31;30;31;30;31;31;30;31;30;31];
  if mod(iyr,4)==0,
    mday(2)=29;
  end

% HYCOM existing with climatology and no Greenland:
  flriva=sprintf('%srivers_%2.2i.a',PTH.data,ntopo);
  flrivb=sprintf('%srivers_%2.2i.b',PTH.data,ntopo);
  fltopo=sprintf('%sdepth_ARCc0.08_%2.2i.nc',PTH.topo,ntopo);
%HYCOM with monthly UCAR runoff and Greenland  
  flrivGa=sprintf('%srivers_%2.2i_NCAR_Gr_%4.4i.a',PTH.dataout,ntopo,iyr); 
  flrivGb=sprintf('%srivers_%2.2i_NCAR_Gr_%4.4i.b',PTH.dataout,ntopo,iyr);

  faold = fopen(flriva,'r','ieee-be');
  fbold = fopen(flrivb,'r');
  fanew = fopen(flrivGa,'w');
  fbnew = fopen(flrivGb,'wt');

  % Get HYCOM topo and grid:
  if ~exist('LAT','var')
    HH  = nc_varget(fltopo,'Bathymetry');
    LAT = nc_varget(fltopo,'Latitude');
    LON = nc_varget(fltopo,'Longitude');
    [m,n]=size(HH);
    [DX,DY] = sub_dx_dy(LON,LAT);
    ACell = DX.*DY;
    IDM=n;
    JDM=m;
    IJDM=IDM*JDM;
    npad=4096-mod(IJDM,4096);
    toto=ones(npad,1);

    icst=find(HH<-1 & HH>=-150);
    Xcst=LON(icst);
    Ycst=LAT(icst);
  end
  

  % Write heading:
  fprintf('Heading %s \n',flrivb);
  for nl=1:5
    aa=fgetl(fbold);
    fprintf('--> %s\n',aa);

    if nl==1
      aa=sprintf('Monthly river as precip., m/s, Arctic N.Atl, N.Pacific 81 rivers ');
    end
    if nl==2
 aa=sprintf('Yenisey,Lena,Ob,Mackenzie,Pechora,Khatanga,Kolyma,S.Dvina,Taz,Indigirka;');
    end
    if nl==3
 aa=sprintf('Mackenzie,Olenek,Yana,Mezen,Yana,Onega,Anabar,Pyasina,Anderson,Yucon;');
    end
    if nl==4
      aa=sprintf('Dai NCAR data & Greenland Runoff, Bamber, Year %4.4i',YearGr);
    end
    fprintf('<-- %s\n',aa);

    fprintf(fbnew,[aa,'\n']);
  end

% Both NCAR and Green. rivers are monthly
  for k=1:12
    mo = k; 
  % Read *.b:  
    aa=fgetl(fbold);
    disp(aa);

  % Read HYCOM runoff:
    A=fread(faold,IJDM,'float32'); % read 2D field
    dm1=fread(faold,npad,'float32');  % Padding = size(toto)
    if size(dm1) ~= size(toto)
      error('Padding in HYCOM file ???');
    end
    toto=dm1;
    clear dm1

  %  I=find(A>1e10);  % there are no land marks in rivers
  %  A(I)=NaN;
    A=reshape(A,IDM,JDM)';

% Get NCAR rivers
    Rncar = sub_ncar_rivers2arc(HH,LAT,LON,ACell,YearGr,mo);
    maxHR=max(max(Rncar));  % max Riv, m/s, should be Yenisey ~250-600 km3/yr
    A = Rncar; 
    
  % Get Greenland:
    dnmb = datenum(YearGr,k,1);
    it   = find(TM==dnmb);
    dmm  = squeeze(Ds(it,:,:));
    In = find(dmm<0);
    if ~isempty(In);
      fprintf('  !!!!    Di<0: %8.3f\n',dmm(In)); % can be bug in Bamaber's data
    end
    Di   = abs(dmm).*LGr; % Solid disch, km3/mo
    dmm  = squeeze(Rt(it,:,:));
    Rti  = abs(dmm).*LGr; % Gr tundra
    dmm  = squeeze(Rg(it,:,:));
    Rgi  = abs(dmm).*LGr; % Greenland meltwater
    
    FWF  = Di+Rti+Rgi; % disch+tundra+GrMelt
%    Rgr=GR(k).runoff;
    Rgr  = FWF*1e9/(3600*24*mday(k));  % km3/mo -> m3/s  
    Rgr(Rgr==0) = nan;

  % Collocate Greenland rivers with HYCOM coastline:
  % points that close to each other - combine in 1 source
    Rmin=50000;  % integrate runoff inside 
    II=find(~isnan(Rgr));
    xmm=Xgr(II);
    ymm=Ygr(II);
    rmm=Rgr(II);
    clear RgrInt
    nr=0;
    for ip=1:length(II);
      x0=xmm(ip);
      y0=ymm(ip);
      if ~isnan(x0),
	xmm(ip)=nan;
	ymm(ip)=nan;
	D=distance_spheric_coord(y0,x0,ymm,xmm);
	iR=find(D<=Rmin);
	V=0;
	if ~isempty(iR),
	  V=sum(rmm(iR));
	end
	V=V+rmm(ip);
	nr=nr+1;
	RgrInt.Vm3sec(nr)=V;
	RgrInt.x0(nr)=x0;
	RgrInt.y0(nr)=y0;
	xmm(iR)=nan;
	ymm(iR)=nan;
      end    
    end
  % 
  % Check total runoff: m3/s
    vtot=nansum(nansum(Rgr));
    vtot2=sum(RgrInt.Vm3sec);
    if (abs(1-vtot/vtot2)>0.001)
      error('Total river runoff not conserved');
    end

  % find HYCOM location
  % convert to m/s
  % redistribute over 9 points
    clear HRcell
    for ip=1:nr
      if mod(ip,10)==0
	fprintf('Placing Gr.runoff on HYCOM grid, ip=%i\n',ip);
      end

      x0=RgrInt.x0(ip);
      y0=RgrInt.y0(ip);
      v0m3=RgrInt.Vm3sec(ip);
  % Find closest coast point:    
      D=distance_spheric_coord(y0,x0,Ycst,Xcst);
      imc=find(D==min(D));
      xcs=Xcst(imc);
      ycs=Ycst(imc);
      D=distance_spheric_coord(ycs,xcs,LAT,LON);
      [jm,im]=find(D==min(min(D)));
  %   Area of the grid cell:
%      dx=distance_spheric_coord(LAT(jm,im),LON(jm,im),...
%			       LAT(jm,im+1),LON(jm,im+1));
%      dy=distance_spheric_coord(LAT(jm,im),LON(jm,im),...
%			       LAT(jm+1,im),LON(jm+1,im));
%      Acell=dx*dy;
      Acell = ACell(jm,im);
      v0=v0m3/Acell;  %m3/s-> m/s - river flux, as precip in HYCOM
  % Find 9 grid points around:
      ncl=1;
      Ntot=9;
      HRcell(ip).I(1)=im;
      HRcell(ip).J(1)=jm;
      HRcell(ip).X(1)=LON(jm,im);
      HRcell(ip).Y(1)=LAT(jm,im);
      HRcell(ip).v0(1)=v0/Ntot;
      HRcell(ip).Acell(1)=Acell;
      di=0;
      clear i
%keyboard      
      while ncl<Ntot
	di=di+1;  % #gr pnts from i,j - search radius
	LL=2*di;
	ic0=im+di;
	jc0=jm+di;
	ivct=-1+0*i;  % dir vector
	for kk=1:4*LL % grid pnts around i,j
  %	kk
	  if HH(jc0,ic0)<0
	    ncl=ncl+1;
	    HRcell(ip).I(ncl)=ic0;
	    HRcell(ip).J(ncl)=jc0;
	    HRcell(ip).X(ncl)=LON(jc0,ic0);
	    HRcell(ip).Y(ncl)=LAT(jc0,ic0);
	    HRcell(ip).v0(ncl)=v0/Ntot;
	    HRcell(ip).Acell(ncl)=Acell;
	  end
	  if ncl==Ntot, break; end;
	  if mod(kk-1,LL)==0 & kk>1,
	    ivct=ivct*i;
	  end
	  ic0=ic0+real(ivct);
	  jc0=jc0+imag(ivct);
	end
	if di>Ntot+1, 
	  error('Could not find %i ocean pnts Runoff',Ntot);
	end
      end
    end;  % ip
  % Check the total runoff   
    vtot3=0;
    vmax=0;
    for ip=1:nr
      Acell=HRcell(ip).Acell(ncl);
      v0m3=sum(HRcell(ip).v0(1:Ntot))*Acell;
      vtot3=vtot3+v0m3;
      v0mx=max(HRcell(ip).v0(1:Ntot))*Acell;
      vmax=max([v0mx,vmax]);
      if vmax==v0mx,
	iGrmax=ip;
      end
    end
    if ~exist('Acell_mn','var')
      Acell_mn=mean(HRcell(ip).Acell(ncl));
    end

    fprintf('Green. Runoff m3/s: orig %8.6d, in HYCOM %8.6d\n',vtot,vtot3);
    if abs(1-vtot3/vtot)>0.001,
      error('River runoff is missing in vtot3');
    end

    if vmax>2*maxHR*Acell_mn,
      fprintf('ERR: Check Greenland River #%i\n',iGrmax);
      fprintf('ERR: Gr. runoff = %7.2f m3/s\n',vmax);
      fprintf('     Arctic max ruonff = %7.2f m3/s\n',maxHR*Acell_mn);
      error('  Greenland river is too high ...');
    end



  % Add to HYCOM rivers:
    for ip=1:nr
      for jpp=1:Ntot
	ic0=HRcell(ip).I(jpp);
	jc0=HRcell(ip).J(jpp);
	v0 =HRcell(ip).v0(jpp);
	A(jc0,ic0)=A(jc0,ic0)+v0;
      end
    end

    chck_plt=logical(0);
    if chck_plt
      figure(1); clf;
      hmm=HH;
      hmm(LON<-80)=nan;
      hmm(LON>40)=nan;
      hmm(LAT<50)=nan;

      II=find(~isnan(Rgr));
      xmm=Xgr(II);
      ymm=Ygr(II);
      plot(xmm,ymm,'b.');  % Orig. Gr. runoff loc.
      hold on;
      plot(RgrInt.x0,RgrInt.y0,'ro'); % integerated sources 
      plot(RgrInt.x0(iGrmax),RgrInt.y0(iGrmax),'gd',...
	   'MarkerSize',10,'linewidth',2); % maxrunoff
      
% Plot locations of FW sources in HYCOM      
      nsr = length(HRcell); % # of FW sources in HYCOM
      for js=1:nsr
        xx = HRcell(js).X;
        yy = HRcell(js).Y;
	plot(xx,yy,'c.','MarkerSize',11); %
      end
      
      contour(LON,LAT,hmm,[0 0],'k');
%      axis('equal');
      
      stt=sprintf('%i/%2.2i, Total runoff, %6.2f Sv',year,k,vtot*1e-6);
      title(stt,'fontsize',12);
      set(gca,'xlim',[-80 0],'ylim',[57 85]);
      set(gcf,'Color',[1 1 1]);
      btx = 'addGreenland2011_2016_hycom.m';
      bottom_text(btx,'pwd',1);

  % Runoff stat:    
      figure(2); clf;
      Ina=find(A>0);
      hist(A(Ina),100);
      xlabel('m/s');
      ntr=nansum(nansum(A))*Acell_mn*1e-6;
      stt2=sprintf('All runoffs+Green., Total Arctic runoff, %6.2f Sv',ntr);
      title(stt2,'fontsize',12);
      set(gcf,'Color',[1 1 1]);

  % Runoff map:
      figure(3); clf
      lA=log10(A);
      lA(A==0)=nan;
      contour(hmm,[0 0],'k');
      hold on;
      pcolor(lA); shading flat;
      caxis([-10 -2]);
      title('Log10(runoff), m/s');
      colorbar
keyboard
    end

  % Write in HYCOM *a and *b files:
  % Note: need to change min/max in *b consistent
  % with *a
    fprintf('Writing HYCOM files %s\n',flrivGa);  
    fprintf('Writing HYCOM files %s\n\n',flrivGb);  

    dmm=A';
    dmm=reshape(dmm,IJDM,1);
    minh=min(dmm);
    maxh=max(dmm);
    is=strfind(aa,'=');
    ch=aa(1:is);
    anew=sprintf('%s%3i    %10.7E   %10.7E',ch,k,minh,maxh);
    fprintf(fbnew,[anew,'\n']);
    fwrite(fanew,dmm,'float32','ieee-be');
    fwrite(fanew,toto,'float32','ieee-be');

  end;  % for k - months

  fclose(fanew);
  fclose(fbnew);
  fclose(faold);
  fclose(fbold);

end;














