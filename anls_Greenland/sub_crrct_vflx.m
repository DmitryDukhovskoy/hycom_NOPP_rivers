function VCT = sub_crrct_vflx(GC,U,V,dH,DX,DY,HH,S,T,Tref,hflg);
% Volume flux is not conserved
% in HYCOM mean files
% calculate net vol flux correction 
% (+Vct*dh*dz - correction for flux in 1 grid cell) 
% across the contour
% assuming it should be ~0
% disbalance distribute over
% all grid points
% hflg = 1 - weight transport misbalance on rate of 
%            negative heat flux
% Apply correction to the regions of
% intense heat loss in deep layers - obvious bias
% due to unrealistically high transport
%      = 0 - distribute evenly
%
%dH(isnan(dH))=0;
Hs = GC.Hbottom;
[nlr,mm,nn] = size(U);
II = GC.cntr_Iindx;
JJ = GC.cntr_Jindx;
np = length(II);
ZZ = zeros(nlr,np);
Vflx = zeros(1,np)*nan; % volume flux for checking
Asgm = zeros(1,np);
Vct = 0; 
Vct0= 0;
VCT = ZZ*0;

cc = 0;
VF=1e9;
VF0 = VF;
add = 1e-8;
while abs(VF)>1e-1 | VF<0;
%Asgm2 = zeros(1,np);
  cc=cc+1;
  if cc>50, error('Endless loop, VF =%8.6f\n',VF); end

  clear Vflx
  
  FSGM = sub_vhFlx_sgm(GC,U,V,S,T,dH,DX,DY,HH,VCT,Tref);
  

  Vflx = FSGM.Vflx_Dintgr_m3_s;
  Hflx = FSGM.Hflx_grcell_W;
  ZZ   = FSGM.ZZ;
  Asgm = FSGM.Asgm_m2;
  
  VF   = nansum(Vflx);
  if abs(VF)<1e-1 & VF>=0, break; end;
%
% To avoid oscillating sollutions: 
  if sign(VF)~=sign(VF0) & cc>2, 
    VctN=0.5*(Vct+Vct0);
  else
    VctN=0;
  end
  VF0=VF;
  
  
  %Asct = nansum(Asgm); % lateral surface area, m2
  Asct = nansum(nansum(Asgm)); % lateral surface area, m2
%keyboard

  % Correction of u (m/s) or flux per 1 m2 of area
%  if sign(VF)~=sign(VF0), 
%    add=add*0.5;
%  end
%  VF0=VF;
%
if abs(VctN)>0
  Vct = VctN;
else
  Vct  = Vct-(VF+sign(VF)*0.1)/Asct;
end
Vct0 = Vct;

VCT  = Hflx*0;
if hflg ==1
  Hflx(1:20,:)=0;
  Ing = find(Hflx<0);
  Mng = median(Hflx(Ing));
  Ing = find(Hflx<Mng);
  Asgm(Asgm==0)=nan;
%  HFm2= Hflx./Asgm;
  wgt = zeros(size(Hflx));
  wgt(Ing) = abs(Hflx(Ing));
  nsm = sum(sum(wgt));
  wgt = wgt/nsm;
  
  VCT = Vct*Asct.*wgt./Asgm; % m/s
  
%  Vct = Vct-VF/Asct-sign(VF)*add;
else
  VCT = VCT+Vct;
end

  
  fprintf('cc=%i, VF=%10.8d Vct=%10.8d, max|Vct|=%10.8d\n\n',...
	  cc,VF, Vct, max(max(abs(VCT))));
%keyboard

end


%keyboard
%v12=1.12186749e-01
%vf12=533.86613576, vf34=489.98154984, vflx=3.17296393e+05, VF=6.78946018e+04
%cc=15, VF=3.11165911e-02 Vct=3.35958640e-03
return