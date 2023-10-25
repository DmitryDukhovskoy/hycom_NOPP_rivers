function Rnew = sub_ncar_rivers2arc(HH,LT,LN,ACell,year,mo);
% Create river 2D field for 
% given year, month - itime( time index)
% Putting river sources  to the ARCc grid
% Locate rivers near the coast line
% distribute over N pnts - proportional to max runoff
% runoff is unevenly distributed - max is closest to original 
% point -source
% RVR -NCAR river array, structured, with river pnt sources
pthmat = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
%fmat   = sprintf('%sncar_rivers_Arctic_1993-2015.mat',pthmat);
fmat   = sprintf('%sncar_rivers_Arctic_1993-2016.mat',pthmat);
load(fmat); % <-- RVR

TM = RVR.TM;
itime = find(TM==datenum(year,mo,1));
Rindx = RVR.Riv_indx;
Nrv = length(Rindx);

Iocn = find(HH<-0.1);
Iland= find(HH>=0);
Rnew = HH*0;
Rnew(HH>=0) = nan;

[mm,nn]=size(HH);
[IIs,JJs]=meshgrid([1:nn],[1:mm]);

% Determine # of points to distribute runoff
% Use overall max - large rivers over more grid points
% of points, assuming river distributed in NxM
% 3x2, 5x3, etc.


for ir=1:Nrv
  rname = RVR.Riv_name(ir,1:7);
  x0=RVR.Lon(ir);
  y0=RVR.Lat(ir);
  D=distance_spheric_coord(LT,LN,y0,x0);
  I0 = find(D==min(min(D)));
  [j0,i0] = ind2sub(size(HH),I0);
  acell = ACell(I0);
  
% Find closest ocean pnt, if needed
% If pnt is sea pnt - move to the coast first
  iL=0;
  if HH(j0,i0)<0
    D=distance_spheric_coord(LT(Iland),LN(Iland),y0,x0);
    imm = find(D==min(min(D)));
    iL = Iland(imm);
    [j0,i0] = ind2sub(size(HH),iL);
  end    
  
  D=distance_spheric_coord(LT(Iocn),LN(Iocn),y0,x0);
  imm = find(D==min(min(D)));
  iC = Iocn(imm);
  [jc,ic] = ind2sub(size(HH),iC);
  
  QRv = RVR.Qflow_m3_sec(:,ir);
  Qr = QRv(itime);
  if Qr<=0,
    fprintf('ERR: Check river runoff, Qr: %8.3f',Qr);
    keyboard
  end
  
%  Npnts = ceil(log(max(QRv)));
% Approximate # of points around the river source:
%  Npnts = ceil(log(Qr));
  Qrms = Qr/acell; % m3/s->m/s
  Npnts = round(Qrms/0.5e-4); % ~max runoff in a grid cell 
  if Npnts<9, Npnts=9; end;
%  if Npnts>81, Npnts=81; end;

% Find how many points can be distributed around the source, 
% for given coastline
% move them "around" the source like in a spiral
% further away points have less runoff - to reduce 
% the gradient jump
  NPP = struct;
  ncl=1;
  NPP.di(1)=1;
  NPP.I(1)=ic;
  NPP.J(1)=jc;
  di = 0;
  clear i
  while ncl<Npnts;
    di=di+1;  % #gr pnts from i,j - search radius
    LL=2*di;
    ic0=ic+di;
    jc0=jc+di;
    ivct=-1+0*i;  % dir vector
    for kk=1:4*LL % grid pnts around ic,jc
  %	kk
      if HH(jc0,ic0)<0
        ncl=ncl+1;
        NPP.I(ncl)=ic0;
	NPP.J(ncl)=jc0;
        NPP.di(ncl)=di;
      end
    
      if ncl==Npnts, break; end;
      if mod(kk-1,LL)==0 & kk>1,
        ivct=ivct*i;
      end
      ic0=ic0+real(ivct);
      jc0=jc0+imag(ivct);
    end
  end
  fprintf('%s, Qr=%8.1f m3/s, estimated Np %i, actual %i, di=%i\n',...
	  rname, Qr, Npnts,ncl,di);

% Calculate weight for points:
  ndi = max(NPP.di(:));
  Np  = ncl; 
  if Np<0
    error('Found Np %i <9 pnts ...',Np);
  end
  
% For 8 points in 1 di layer  around the source pnt (9pnts) - 
% 0 grad., simply mean average:
% otherwise bin-average for function
% q = Qr-alf*p, where p - point number from 1 to Np
% alf = Qr/Np
% Bin-average is qB(i) = 0.5*(qb1+qb2)
  Qrm3 = [];
  if ndi==1,
    Qrm3(1:Np)=Qr/Np;
  else % do bin-averaging for di=1:ndi
    qm0 = 2*Qr/Np;
    alf = qm0/Np;
    for ndd=1:ndi
      nb1=min(find(NPP.di(:)==ndd));
      nb2=max(find(NPP.di(:)==ndd));
      nb1=nb1-1;
 %     if nb1>1, nb1=nb1-1; end;
 %     qb1 = Qr-alf*(nb1-1);
 %     qb2 = Qr-alf*(nb2-1);
 %     qB = 0.5*(qb1+qb2);
 %     Qrm3(nb1:nb2)=qB/(nb2-nb1+1); % evenly distribute within 1 di layer
 %     Qt = qm0*(nb2-nb1+1)-alf/2*(nb2^2-nb1^2+1);
      Qt = qm0*(nb2-nb1)-alf/2*(nb2^2-nb1^2);
      qav=Qt/(nb2-nb1); % average vol for bin1 = di=1
      Qrm3(nb1+1:nb2) = qav;
    end
  end
  
% check:
  dmm = sum(Qrm3);
  if abs(1-dmm/Qr)>1e-6
    fprintf('sub_ncar_rivers2arc: Volume by points not conserved ...\n');
    fprintf('Check ...\n');
    keyboard;
  end
%keyboard
  RR = Rnew;
  RR(HH<0)=0;
  RR(HH>=0)=nan;
  if (nansum(nansum(RR))~=0)
    error('Initial RR is not zero ...');
  end
  
  for kk=1:Np
    i0=NPP.I(kk);
    j0=NPP.J(kk);
    RR(j0,i0)=Qrm3(kk);
  end
  
% 9-box filter
% Note that when land is around, filtering can 
% produce more mass, as near-coast points
% will be smoothed unequally compared to the
% offcoast - points, thus have to modify 
% "box" to any closest 9 pnts 
% Spread rivers over borader region
% Locate Nf ocean points closest to the 
% river source ic,jc:
  fbox = 0; % fbox=1 standard box filter, =0 - smooth over N closest points
  if fbox==0
    ii1=min(NPP.I(:))-4;
    ii2=max(NPP.I(:))+4;
    jj1=min(NPP.J(:))-4;
    jj2=max(NPP.J(:))+4;
    Nf = (ii2-ii1+1)*(jj2-jj1+1);
    D = sqrt((IIs-ic).^2+(JJs-jc).^2);
    D(HH>0)=1e9;
    cp=0;
    while cp<Nf
      iF=find(D==min(min(D)),1);
      cp = cp+1;
      IFlt(cp,1) = iF; 
      D(iF) = 1e9;
    end
    [jF,iF] = ind2sub(size(HH),IFlt);

    Nbx = 9; % N-pnt box filtering
    for iflt=1:3
      R0=RR;
      dmm = R0(IFlt);
      for kk=1:Nf
  %      i0=NPP.I(kk);
  %      j0=NPP.J(kk);
	i0=iF(kk);
	j0=jF(kk);
	D = sqrt((iF-i0).^2+(jF-j0).^2);
	[sD,sI] = sort(D);
	sI = sI(1:Nbx);
	RR(j0,i0)=mean(dmm(sI));
      end
    end
  else % Standard 9-box filter 
    ii1=min(NPP.I(:))-5;
    ii2=max(NPP.I(:))+5;
    jj1=min(NPP.J(:))-5;
    jj2=max(NPP.J(:))+5;
    for iflt = 1:3
      R0=RR;
      for ii=ii1:ii2
	for jj=jj1:jj2;
	  if HH(jj,ii)>=0; continue; end;
	  dmm = R0(jj-1:jj+1,ii-1:ii+1);
	  dav = nanmean(nanmean(dmm));
	  RR(jj,ii) = dav;
	end
      end
    end
  end
%
% Adjust total volume flux after smoothing
  Qtot=nansum(nansum(RR));
  dQ = Qr-Qtot;
%  if dQ<0, % spatial averaging may produce extra "mass"
%    fprintf('dQ<0, %8.3d\n',dQ);
%    keyboard;
%  end
  Ib = find(RR>0);
  dQcell = dQ/length(Ib);
  RR(Ib) = RR(Ib)+dQcell;
  Ing = find(RR<0);
% Reduce surplus flux by eliminating grid cells with small runoff
% those were produced during smoothing
  Nng = length(Ing);
  while Nng>0
    spQ = abs(sum(RR(Ing)));
    RR(Ing) = 0;
    Isb = find(RR>0);
    RR(Isb) = RR(Isb)-spQ./length(Isb);
    Ing = find(RR<0);
    Nng = length(Ing);
    if isempty(Nng); break; end;
  end
% 
% Check again:
  Qtot = nansum(nansum(RR));
  if abs(1-Qtot/Qr)>1e-3
    fprintf('sub_ncar_rivers2arc: Volume by points not conserved ...\n');
    fprintf('Second Check failed...\n');
    keyboard;
  end

% Check negative runoff:
% Some "mass" may be produced during filtering
% due to uneven averaging with land around
  In = find(RR<0);
  if ~isempty(In),
    fprintf('Found %i pnts with Q<0\n',length(In));
    [jn,in] = ind2sub(size(HH),In);
    keyboard
  end
  
  
% ------------------------------
  fchck=0;
  if fchck>0
    figure(10); clf;
    contour(HH,[0 0],'k');
    hold on;
    plot(i0,j0,'b*');
    plot(ic,jc,'r*');
    
    pcolor(RR); shading flat;
    
    set(gca,'xlim',[ic-40 ic+40],...
	    'ylim',[jc-40 jc+40],...
	    'xgrid','on','ygrid','on');
  end
% ------------------------------
  

% Normalize by cell area:
  RR = RR./ACell;
  fprintf(':: River %s on ARCC grid, %i/%i, max %6.3d min %6.3d m/s\n',...
	rname,year,mo,max(max(RR)), min(min(RR(RR>0))));

  Rnew = Rnew+RR;
end  % river sources
  

Rnew(isnan(Rnew))=0; % NO LAND mask in rivers

return