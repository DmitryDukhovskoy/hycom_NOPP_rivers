function Fnew = sub_interp_fld32to41(vname,finaOLD,finbOLD,...
			      HH,ZZn,ZZo,nsigma,Lx_new, Lx_old);
% Interpolate 3D fields from 32 -> 41 layers
% Interpolation is "mass" weighted, where mass is the layer 
% thickness
% ZM*, ZZ* - mean depths and interface depths of vertical layers
% Zdeep,dPdeep - dpeth/thicknesses of the new layers in deep ocean region
%                assuming they are fixed in time
%     do not apply those in the shallow region where sigma-layers
%     are generated
% Lx_new & Lx_old  - layer # in 41-layer grid where target densities
%        start to match the 32 layer grid = tdens(Lx_old)
%  -------------------------------------------------- 
%  For each interface in new t.dens. find
%  nearest upper and bottom interfaces in old t.dens.
%
%   Old v. layers             New v. layers
%   -------------
%       Z(kold) *           ----------   
%                            Zn(knew)   *
%   -------------
%
%     Z(kold+1) *           ----------   
%                            Zn(knew+1) *
%   -------------

hg=2^100;
[F,n,m,l] = read_hycom(finaOLD,finbOLD,vname);% read in 32-layer field
F(F>0.0001*hg)=nan;
IOcean=find(HH<0);
np=length(IOcean);

[nz,mm,nn]=size(ZZn);
NLold=l;
NLnew=nz-1;

Fnew = zeros(NLnew,m,n)*nan;
fprintf('Start interpolation %s, Npoints=%i\n',vname,np);

po=-1;
for ik=1:np
  pn=ik/np*100;
  if floor(pn)>po
    fprintf('>> %s:  Interp Done %6.3f%%\n',vname,pn);
    po=floor(pn);
  end
  
  [j,i] = ind2sub(size(HH),IOcean(ik));
%j=60;
%i=403;
%i=171;
%j=466;
i=184;
j=464;

% 32-layer fields:  
  h0  = HH(j,i);
%  dp0 = squeeze(DP_old(:,j,i));
  fold = squeeze(F(:,j,i));
  Zold=ZZo(:,j,i);
  
% 41-layer fields:
  Znew=ZZn(:,j,i);
%  sgmDP=Znew(nsigma+1);
%  sgmH  = abs(1-abs(sgmDP)./abs(h0)); % sigma-coord. points
% In theory, interpolation is needed
% only in the new layers added above Lx_new
% and all layers below should be identical
% to the old 32-layer grid (below Lx_old)
% however, in several locations there was
% some mismatch, so have to interpolate
% over all layers
  fnew = ones(NLnew,1)*fold(end);
  fnew(Lx_new:end)=fold(Lx_old:end-1);
%keyboard  
  for k=1:NLnew
    dHn=abs(Znew(k+1)-Znew(k)); % L. thickness
    if dHn<1e-3, continue; end
    z1_new=Znew(k);
    z2_new=Znew(k+1);
    z2_new=floor(abs(z2_new*100))*(-0.01); % trunc. errors bottom depth
    
% Find corresponding layer surfaces in old grid
% that contain new surfaces
    k1=max(find(Zold>=z1_new));  % top HYCOM surface
    k2=max(find(Zold>z2_new)); % btm HYCOM surface
    
    z1_old=Zold(k1);
    z2_old=Zold(k2);
% How many complete old z-levels inside the new layer:
    nZ = k2-k1-1;

    chck=0;
    if chck==1
      figure(10); clf;
%      xx=ones(nZ,1);
      hold on;
      for kz=1:NLold
        plot([0.98 1.02],[Zold(kz) Zold(kz)],'b.-');
      end
      xx=ones(NLnew,1);
      for kz=1:NLnew
        plot([0.98 1.02],[Znew(kz) Znew(kz)],'r--');
      end
      plot([0.9 1.1],[z1_old z1_old],'b-','linewidth',1.5);
      plot([0.9 1.1],[z2_old z2_old],'b-','linewidth',1.5);
      plot([0.9 1.1],[z1_new z1_new],'m:','linewidth',1.5); % HYCOM sigma0 interface
      plot([0.9 1.1],[z2_new z2_new],'m:','linewidth',1.5);
    end
    
% interpolate into a new layer
    if k1>NLold, 
      fprintf('k1>NLold - error in layer index\n');
      keyboard; 
    end;

    if nZ<0 % new layer is completely within the old layer:
      fnew(k)=fold(k1);
    else 
% top and bottom of new layer are in different layers of the old grid
% nZ = can be 0 or >0 complete "old" layers within the new layer
% Delta thickness of partial layers      
      dZS1=abs(z1_new-Zold(k1+1));
      dZS2=abs(z2_new-Zold(k2));
      pU1=fold(k1)*dZS1;
      if k2<=NLold
        pU2=fold(k2)*dZS2;
      else % hit bottom
	pU2=0;
      end
      
      
      dZfull=0;
      Ufull=0;
      cc=0;
      for kf=k1+1:k2-1
	cc=cc+1;
	dZ1=abs(Zold(kf)-Zold(kf+1));
	dZfull=dZfull+dZ1;
	Ufull=Ufull+fold(kf)*dZ1;
      end
% Check total thickness of the layer:
      dZtot=dZS1+dZfull+dZS2;
      err=abs(abs(dZtot/(z2_new-z1_new))-1);
      if err>1e-2
	error('*** Check total layer thickness over %i nZ = %8.2m\n',...
	      nZ,dZtot);
      end
% Interpolated value in the new layer:
      fnew(k)=(pU1+Ufull+pU2)./dZtot;
    end;  % if nZ<0
%      keyboard
  end;    % for k=1:Lx_new-1   v. layers

%
% Check interpolation:
  dZold=diff(Zold);
%  ibtm=min(find(Zold<=h0+1e-6))-1;
%  if isempty(ibtm), ibtm=NLold; end;
  Utot_old=sum(fold.*abs(dZold));
  
  dZnew=diff(Znew);
  Utot_new=sum(fnew.*abs(dZnew));
  
  if abs(1-Utot_new/Utot_old)>0.1 & ...
	abs(Utot_new-Utot_old)>1e-2,  % small values - truncation error
    fprintf('ERR: Total value is not conserved after interpolation\n');
    keyboard;
  end
  
keyboard
  Fnew(:,j,i)=fnew; 
  
end % ocean points


f_chck=0;
if f_chck>0
%  io=8000;
  II=IOcean(io);
  [jj,ii]=ind2sub(size(HH),II);
%  jj = m-1;
%  ii = 465;
  Zold = squeeze(ZZo(:,jj,ii));
%  Zmold= squeeze(ZMo(:,jj,ii));
  fold = squeeze(F(:,jj,ii)); % U on Z
  fnew = squeeze(Fnew(:,jj,ii));
  hb=HH(jj,ii);

  figure(10); clf;
  
  axes('Position',[0.08 0.06 0.35 0.85]);
  hold on;
%  plot(fold,Zmold,'g-'); % Profile, old v.layers
  % Plot U in Z-level interpolated
  % into mid-depth grid layers
  for k=1:NLold-1
    u=fold(k);
    plot([u u],[Zold(k) Zold(k+1)],'r-');
  end
  % Plot profile on new v.grid
  for k=1:NLnew
    u=fnew(k);
    plot([u u],[Znew(k) Znew(k+1)],'b--');
  end
  plot([min(fold) max(fold)],[hb hb],'k--');
  title('U, old ZZ(r), New ZZ(b)' );
  set(gca,'ylim',[1.05*hb 0],...
	  'ygrid','on')
end




return