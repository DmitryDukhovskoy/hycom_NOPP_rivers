% Rivers in HYCOM are approximated as bogus precipitation
% units are m/s
% see forfun.f:
%c --- initialize input of river (precip bogas) forcing field
%c --- units of rivers are m/s    (positive into ocean)
%c --- rivers is always on the p grid.
%c --- I/O and array I/O unit 918 is reserved for the entire run.
%c --- all input fields must be defined at all grid points    
%
% Use Greenland runoff data <-- NOTE: Bamber's data are TOTAL FW
% flux, it includes ice discharge
% 
% on a 5km grid - integrate 
% runoff over some distance along the coast
%
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup

format long g
clear all
close

regn = 'ARCc0.04';
ntopo= 17;      % HYCOM topo version
TV   = sprintf('%2.2iDD',ntopo);
%YearGr=2009;  % Gr. river runoff river

%fprintf('River Runoff Year: %i\n',YearGr);

hg=2^100; 

%PTH.river = '/Net/mars/ddmitry/hycom/ARCc0.04/force/rivers/';
PTH.river='/Net/ocean/ddmitry/arctic_AOregimes/data/Greenland_rivers/';
PTH.data  = '/Net/mars/ddmitry/hycom/ARCc0.04/force/rivers/';
PTH.topo  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';

for YearGr=1992:2010
  fprintf('River Runoff Year: %i, Topography: %i\n',YearGr,ntopo);

  mday=[31;28;31;30;31;30;31;31;30;31;30;31];
  if mod(YearGr,4)==0,
    mday(2)=29;
  end


% HYCOM:
  flriva=sprintf('%srivers_%s.a',PTH.data,TV);
  flrivb=sprintf('%srivers_%s.b',PTH.data,TV);
%  fltopo=sprintf('%sdepth_%s_%2.2i.nc',PTH.topo,ntopo);
  fltopo=sprintf('%sdepth_%s_%s.nc',PTH.topo,regn,TV);
%HYCOM with Greenland  
  flrivGa=sprintf('%srivers_%s_Greenland_%4.4i.a',PTH.data,TV,YearGr); 
  flrivGb=sprintf('%srivers_%s_Greenland_%4.4i.b',PTH.data,TV,YearGr);

  faold = fopen(flriva,'r','ieee-be');
  fbold = fopen(flrivb,'r');
  fanew = fopen(flrivGa,'w');
  fbnew = fopen(flrivGb,'wt');

  % Greenland:
  fgrgr=sprintf('%sGreenland_grid.mat',PTH.river);
  fgrrv=sprintf('%sGreenland_runoff_monthly-%4.4i.mat',PTH.river,YearGr);
  % Greenland runoff and grid:
  load(fgrrv);
  load(fgrgr);
  Xgr=GRgrd.LN;
  Ygr=GRgrd.LT;
  clear GRgrd

  % Get HYCOM topo and grid:
  HH  = nc_varget(fltopo,'Bathymetry');
  alat = nc_varget(fltopo,'Latitude');
  elon = nc_varget(fltopo,'Longitude');
  [m,n]=size(HH);
  IDM=n;
  JDM=m;
  IJDM=IDM*JDM;
  npad=4096-mod(IJDM,4096);
  toto=ones(npad,1);

  icst=find(HH<-1 & HH>=-150);
  Xcst=elon(icst);
  Ycst=alat(icst);

  % Write heading:
  for nl=1:5
    aa=fgetl(fbold);
    disp(aa);

    if nl==4
      aa=sprintf('Greenland Runoff, Bamber, Year %4.4i',YearGr);
    end

    fprintf(fbnew,[aa,'\n']);
  end


  % Both HYCOM and Green. rivers are monthly
  for k=1:12

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
    maxHR=max(max(A));  % max Riv, m/s, should be Yenisey ~600 km3/yr

  % Get Greenland:
    Rgr=GR(k).runoff;
    Rgr=Rgr*1e9/(3600*24*mday(k));  % km3/mo -> m3/s  

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
      D=distance_spheric_coord(ycs,xcs,alat,elon);
      [jm,im]=find(D==min(min(D)));
  %   Area of the grid cell:
      dx=distance_spheric_coord(alat(jm,im),elon(jm,im),...
			       alat(jm,im+1),elon(jm,im+1));
      dy=distance_spheric_coord(alat(jm,im),elon(jm,im),...
			       alat(jm+1,im),elon(jm+1,im));
      Acell=dx*dy;
      if Acell < 1e-3; 
	fprintf('ERR: Area = 0 %d\n',Acell);
	keyboard
      end
      
      v0=v0m3/Acell;  %m3/s-> m/s - river flux, as precip in HYCOM
  % Find 36 grid points around - to make it
  % comparable with ARCc0.08 (9pnts):
      ncl=1;
      Ntot=36;
      HRcell(ip).I(1)     = im;
      HRcell(ip).J(1)     = jm;
      HRcell(ip).X(1)     = elon(jm,im);
      HRcell(ip).Y(1)     = alat(jm,im);
      HRcell(ip).v0(1)    = v0/Ntot;
      HRcell(ip).Acell(1) = Acell;
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
	    HRcell(ip).X(ncl)=elon(jc0,ic0);
	    HRcell(ip).Y(ncl)=alat(jc0,ic0);
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
    vtot3  = 0;
    vmax   = 0;
    iGrmax = 0;
    for ip=1:nr
      Acell=HRcell(ip).Acell(ncl);
      v0m3=sum(HRcell(ip).v0(1:Ntot))*Acell;
      vtot3=vtot3+v0m3;
      v0mx=max(HRcell(ip).v0(1:Ntot))*Acell;
      VolF(ip,1)=v0m3;
      vmax=max([v0mx,vmax]);
      if vmax==v0mx,
	iGrmax=ip;
      end
    end
    if ~exist('Acell_mn','var')
      Acell_mn=mean(HRcell(ip).Acell(ncl));
    end

    fprintf('Green. Runoff m3/s: orig %d, in HYCOM %d\n',vtot,vtot3);
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
      hmm(elon<-80)=nan;
      hmm(elon>40)=nan;
      hmm(alat<50)=nan;

      II=find(~isnan(Rgr));
      xmm=Xgr(II);
      ymm=Ygr(II);
      plot(xmm,ymm,'b.');  % Orig. Gr. runoff loc.
      hold on;
      plot(RgrInt.x0,RgrInt.y0,'ro'); % integerated sources
      plot(RgrInt.x0(iGrmax),RgrInt.y0(iGrmax),'gd',...
	   'MarkerSize',10,'linewidth',2); % maxrunoff
      contour(elon,alat,hmm,[0 0],'k');
      stt=sprintf('Month %i, Total runoff, %6.2f Sv',k,vtot*1e-6);
      title(stt,'fontsize',12);
      set(gca,'xlim',[-80 0],'ylim',[57 85]);
      set(gcf,'Color',[1 1 1]);

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
    fprintf(fbnew,[aa,'\n']);
    fwrite(fanew,dmm,'float32','ieee-be');
    fwrite(fanew,toto,'float32','ieee-be');

  end;  % for k - months

  fclose(fanew);
  fclose(fbnew);
  fclose(faold);
  fclose(fbold);

end;









