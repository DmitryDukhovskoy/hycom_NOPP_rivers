function R2 = sub_riversARCc08toARCc04(HH1,R1,ACell1,X1,Y1,...
				       HH2,ACell2,X2,Y2);
% Have to conserve the total Volume flux
% i.e. hycomm river flux (m/s)*area_river_pnts = m3/s
%
% Interpolate river sources (R) from the  old grid ARCc0.08
% to the new grid ARCc0.04
% new Topo (HH2)
% old topo = HH1


IR1    = find(R1>0);
nr1    = length(IR1);
rFlxT1 = nansum(nansum(R1.*ACell1)); % total R.flx m3/s
%keyboard
% Onland river sources:
% that need to be moved

R2 = HH2*0;

fprintf(':: Found %i river grid cells in ARCc0.08\n',nr1);
if nr1==0; return; end;

[mm,nn]=size(Y2);
%keyboard
% Each grid point from ARCc0.08 where river > 0
% distribute over n points on ARCc0.04
% using Acell ratio
% some points in ARCc0.04 may overlap
% Note that river fluxes are m/s in HYCOM
% i.e. normalized by cell area
% no need to adjust for different areas
kk=1;
cc=0;
IPTot = [];
cpT = 0;
Atot = 0; % for checking, tot area of ARC04 river points
for ia=1:nr1
  ndone=ia/nr1*100;
  
  if mod(ia,50)==0
    fprintf('Distributing arc08->arc04,   river points done %4.2f%%\n',ndone);
  end
  
  ir0=IR1(ia);
  [j0,i0]=ind2sub(size(R1),ir0);
  xx1=X1(ir0);
  yy1=Y1(ir0);
  a1=ACell1(ir0);
  riv1=R1(ir0);
% find closest river points
% and distribute runoff from 1 grid point from ARCc0.08
% to rpnt in ARCc0.04
% conserving area-integrated riv. Vol flx
  d=distance_spheric_coord(Y2,X2,yy1,xx1);
  irN=find(d==min(min(d)),1);
  [jN,iN]=ind2sub(size(Y2),irN);
  a2=ACell2(irN);
  rpnt = round(a1/a2);
  nx = round(sqrt(rpnt));
  ny = round(rpnt/nx);
%  rpnt = nx*ny;
  
  pp=[1:nx]-ceil(nx/2);
  i1=iN+min(pp);
  i2=iN+max(pp);
  pp=[1:ny]-ceil(ny/2);
  j1=jN+min(pp);
  j2=jN+max(pp);
  
  ji=max([1,j1]);
  j2=min([mm,j2]);
  i1=max([1,i1]);
  i2=min([nn,i2]);
  rpnt=abs(j2-j1+1)*abs(i2-i1+1); % actual # of pnts
  a2tot = sum(sum(ACell2(j1:j2,i1:i2)));
  riv2 = riv1*a1/a2tot;
  
  if rpnt == 0, 
    error('sub_riversARCc08toARCc08:  # river pnts =0'); 
  end;
%  riv2=riv1; 
  
%ke 
  sr=zeros(ny,nx)+riv2;
  sh=HH2(j1:j2,i1:i2);
% Check if the volume flux conserved
% in 1 grid point of ARCc0.08 and N pnts ARCc0.04
  rFlx1=riv1*a1; % m3/sec
  rFlx2=sum(sum(sr*a2));
  
  if abs(1-rFlx2/rFlx1)>0.01
    fprintf('River vol. flux is not conserved\n');
    fprintf('River Arc08 %8.5f, Arc04 %8.5f m3/s\n',...
	    rFlx1, rFlx2);
    keyboard
  end
% 
% Updated rivers in arc04:
  R2(j1:j2,i1:i2)=R2(j1:j2,i1:i2)+sr;
  iR2    = find(R2>0);
  nr2    = length(iR2);
  rFlxT2 = nansum(nansum(R2.*ACell2)); % total R.flx m3/s
%  fprintf('1.  # of Riv. pnts=%i\n',length(iR2));  
  
  
% Check if any river points are onland
  cp=0;
  clear IP
  for ii=i1:i2
    for jj=j1:j2
      cp=cp+1;
      IP(cp,1)=sub2ind(size(Y2),jj,ii);
    end
  end

  if max(max(sh))<0, continue; end; % all points in the ocean
%keyboard
% For check:
  dg=5;
  rr1=R2(j1-dg:j2+dg,i1-dg:i2+dg);
  ss1=HH2(j1-dg:j2+dg,i1-dg:i2+dg);

  for ipp=1:cp
    sh0=HH2(IP(ipp));
%    rtot2=rtot2+R2(IP(ipp));
    if sh0<0, 
      continue; 
    end;
    x0=X2(IP(ipp));
    y0=Y2(IP(ipp));
    D = distance_spheric_coord(Y2,X2,y0,x0);
    D(HH2>=0)=nan; % exclude land
    D(IP)=nan;     % exclude existing points
    ip0=find(D==min(min(D)));
    R2(ip0)=R2(ip0)+R2(IP(ipp)); % add onland river pnt to closest ocean pnt
    R2(IP(ipp))=0;  % runoff in land point = 0
  end
% CHeck if the river vol. flux conserved
% after moving riv. pnt to land
  iR2n    = find(R2>0);
  nr2n    = length(iR2);
  rFlxT2n = nansum(nansum(R2.*ACell2)); % total R.flx m3/s
%  fprintf('2.  # of Riv. pnts=%i\n',nr2n);
%  fprintf('  Vol. Flux 1 = %d, vol. Flux 2 = %d\n',...
%           rFlxT2, rFlxT2n);

  if abs(1-rFlxT2/rFlxT2n)>0.001
    fprintf('Check total runoff in land-relocated river\n');
    fprintf('Old %8g, New %8g\n',rFlxT2,rFlxT2n);
    keyboard;
    
    figure(1); clf;
    axes('Position',[0.08 0.55 0.4 0.4]);
    pcolor(rr1); colorbar

    axes('Position',[0.08 0.05 0.4 0.4]);
    pcolor(rr1n); colorbar
    
  end
  
end;
R2(isnan(R2))=0; % NO LAND mask in rivers
IR2 = find(R2>0 & HH2>=0);
nr2 = length(IR2);

if nr2>0,
  error(' == found %s river pnts onland for ARCc0.04',nr);
end


rFlxT2 = nansum(nansum(R2.*ACell2)); % m/s -> m3/s
dlt = abs(1-rFlxT1/rFlxT2);
if dlt>0.001
  fprintf(' == Mismatch: total River runoff m3/s in final data\n');
  fprintf('Old %8g, New %8g\n',Rtot1,Rtot2);
  fprintf('Mismatch = %10.6d\n',dlt);
  
  iR2   = find(R2>0);
  nr2   = length(IR2);
  adjR  = (rFlxT1-rFlxT2)./nr2;
  R2(iR2) = R2(iR2)+adjR./ACell2(iR2);
  rFlxT2 = nansum(nansum(R2.*ACell2));
  dltN  = abs(1-rFlxT1/rFlxT2);
  fprintf('Adjusted, correction 1 river source m3/s: %10.6d\n',adjR);
  fprintf('Adjusted, Mismatch = %10.6d\n',dltN);
  if dltN>dlt
    error('Adjusted mismatch > old: Adj %10.6d  old %10.6d',...
	  dltN,dlt);
  end  
%  keyboard;
end

% Compare total volume flux m3/s
rFlxT2 = sum(sum(R2.*ACell2));
dltV = abs(1-rFlxT1/rFlxT2);

if dltV>0.001
  fprintf('Mismatch: total Riv.Flx m3/s\n');
  fprintf('Old %8g, New %8g\n',rFlxT1, rFlxT2);
  fprintf('Mismatch = %10.6d\n',dltV);
  keyboard
end

%fprintf(':: Rivers for ARCc0.04 created ...\n');
fprintf(' =========   Done   ==========  \n');

return