function ZZn=sub_adjust_vlayers(ZZn1,ZZo,Lx_new,Lx_old,HH,nsigma,dPdeep);
% First: update new vertical grid
% above Layer Lx_new - fixed depths, as given in 
%       dPdeep -1D depth array of minimum layer thickness, deep
% below Lx_new - layers identically match the 32-layer grid
%         for layer >=Lx_old
% Target densities of the layers Lx_new and Lx_old 
% are the same
% In some locations, Lx_old (=14) happens to be
% shallower than Lx_new(=24) (although should match)
% In this case, stick with the new grid depth
% the old grid will be assumed at the same depth
% during field interpolation
%
% the bottom-most layer combines 2 bottom layers of 32 grid
% Bathymetry for both vertical grids should match!
% 
% ZZo and ZZn1 - negative depths
%
% ZZn1 - depths from any 41-lr field
%    e.g. climatology or 
%    any output file
% It is used for shelf regions
% to get sigma layers that are fixed in time
% z-levels in the upper ocean above Lx_new
% are from blkdat.input dp0k
%
IOcean=find(HH<0);
np=length(IOcean);

% Check sign convention
if min(min(min(ZZn1)))>=0
  ZZn1=-abs(ZZn1);
end
if min(min(min(ZZo)))>=0
  ZZo=-abs(ZZo);
end

% # of layers, old, new:
[nzz,m,n]=size(ZZo);
nlo=nzz-1;
nln=size(ZZn1,1)-1;

ZZn=ZZn1*0;
%ZMn=zeros(nln,m,n);

% Get z-layer thicknesses from blkdat
%pthb='/home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/';
%fnm=sprintf('%sblkdat.input_ARCc0.08_41lev',pthb);
%VL=read_layers_blkdat(fnm);
%dp0k=VL.dp0k;

po=-1;
ZZn=ZZn1;
for ik=1:np
  pn=ik/np*100;
  if floor(pn)>po
    fprintf('Adjusting layers: Done %6.3f%%\n',pn);
    po=floor(pn);
  end
  [j,i] = ind2sub(size(HH),IOcean(ik));
% Test point: deep ocean 
%  j=15;
%  i=370;
%j=15;
%i=370;
%i=667;
%j=2518;
% shelf:
%   j=15;
%   i=250;
%i=21;
%j=2518;
%i=236;
%j=3;
%i=237;
%j=5;
%j=20;
%i=800;
%  keyboard
  
% 32-layer fields:  
  h0  = HH(j,i);
%  dp0 = squeeze(DP_old(:,j,i));
  Zold=ZZo(:,j,i);
  
% 41-layer fields:
% DP and thicknesses need to be changed
% in the deep ocean to make z-level fixed
% for layers < Lx_new
% layers >=Lx_new should match the old grid
% (layers >=Lx_old)
  Znew=ZZn1(:,j,i);
  zn=ZZn1(:,j,i); % shelf: sigma layers(nsigma) - fixed in time
  
  sgmDP = Znew(nsigma+1);
  sgmH = abs(1-abs(sgmDP)./abs(h0)); % sigma-coord. points
%keyboard
  if sgmH>1e-2 % deep ocean, not sigma coord.
    clear zn
    zn(1,1)=0;
    dHb=abs(Znew-h0);
    Ib=max(find(dHb>1e-2));
    if Ib<Lx_new % bottom shallower than Lx_new where two grids match
      zn(2:Ib+1,1)=-cumsum(dPdeep(1:Ib));
      zn(Ib+1:nln+1)=h0;
      k1=nln+1;
      n1=nlo+1;
      zXold=h0-100;
      zXnew=h0;
    else 
      zn(2:Lx_new,1)=-cumsum(dPdeep(1:Lx_new-1)); % intrfaces above Lx_new
% Bottom of Lx_old layer (Lx_old+1) should be collocated with
% Lx_new+1, same depth
% Sometimes they do not match
% Check that bottom interface Lx_old is not shallower than layer Lx_new
% if it is, then ignore Lx_old and match layers (>Lx_old) in both 
% grid only for those Old Layers (32 lr) that are deeper or same depth
% as the Lx_new
% then Layers in the new grid are added as z-level grid until
% the corresponding layer in the old grid becomes
% deeper, after that match the rest of the layers
      k1=Lx_new; % interface counter (not layers)
      n1=Lx_old;
      zXold=Zold(Lx_old+1); % depth of the bottom interface Lx_old layer
%      zXnew=zn(Lx_new)-dPdeep(Lx_new+1);<--???dPdeep(Lx_new) ???
      zXnew=zn(Lx_new)-dPdeep(Lx_new);
      zXnew=max([zXnew,h0]);
    end
    
% zXold should = zXnew, if zXold shallower zXnew
% This means that the targ. dens sigma41(Lx_new)< density at this depth
% Start adding z-levels using z-level spacing dp0k (dPdeep)
% specified in the blkdat until the new (41) hycom layer
% becomes shallower than the old (32) hycom layer
    f_chck=0;
    if f_chck>0
     figure(21); clf;
      hold on;
      for kz=1:n1
	plot([1 1.5],[Zold(kz) Zold(kz)],'b.-');
      end
      plot([1 1.5],[Zold(kz) Zold(kz)],'b-','linewidth',2);
      
      for kz=1:k1
        plot([1.5 2],[zn(kz) zn(kz)],'r--');
      end
      plot([1.5 2],[zn(kz) zn(kz)],'r--','linewidth',2);
    end;
    
    while round(zXold)>round(zXnew)
      k1=k1+1;
      n1=n1+1;
      if k1>42; break; end;
      zn(k1)=zn(k1-1)-dPdeep(k1-1);
      if abs(zn(k1)-h0)<0.01 | zn(k1)<h0; 
	zn(k1:nln+1)=h0; % hit bottom
	k1=nln+1;
	n1=nlo+1;
	break;
      end; % bottom
      zXold=Zold(n1);
      zXnew=zn(k1);

      if f_chck>0
	plot([1 1.5],[zXold zXold],'b-');
	plot([1.5 2],[zXnew zXnew],'r--');
      end
    end  % while loop
    
%    keyboard
    if k1<=nln
      zn(k1+1:nln+1)=Zold(n1+1:end-1); % deep layers unchanged
    end
    
    zn(end)=Zold(end); % match the bottom
    Znew=zn;
    
  end; % deep ocean
%
  if ~exist('zn','var'); zn=Znew; end;
% Plot old & new interface levels:
  f_chck=0; % plot layers
  if f_chck>0
    figure(20); clf;
    hold on;
    for kz=1:nlo+1
      plot([0.98 1.02],[Zold(kz) Zold(kz)],'b.-');
      if kz == Lx_old+1,
	plot([0.98 1.02],[Zold(kz) Zold(kz)],'k.-');
      end
    end

    for kz=1:nln+1
      plot([0.98 1.02],[Znew(kz) Znew(kz)],'r--');
      if kz == Lx_new+1,
	 plot([0.98 1.02],[Znew(kz) Znew(kz)],'g--');
      end

    end
%    
% Plot new levels:
    for kz=1:nln+1
      plot([1.03 1.1],[zn(kz) zn(kz)],'r--');
      if kz == Lx_new+1,
	 plot([1.03 1.1],[zn(kz) zn(kz)],'g--');
      end
    end
    
    text(1.04,0.3,'Adjusted 41 layers');
    
    set(gca,'xlim',[0.96 1.1]);
    title('Old (solid) and New (dashed) v. grid layer interfaces');
    
    keyboard
  end	

%  zm=zn(1:end-1)-0.5*abs(diff(zn));
  ZZn(:,j,i)=zn;
%  ZMn(:,j,i)=zm;
  
end;  % ik - ocean points


return;