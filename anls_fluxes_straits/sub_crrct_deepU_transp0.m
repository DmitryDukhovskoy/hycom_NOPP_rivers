function Vadj = sub_crrct_deepU_transp0();
% Volume flux is not conserved
% in HYCOM mean files
% bogus deep currents
% up to 2.5 m/s
% adjust deep flow only
% correction is weighted by 
% the speed of the deep flow
%
% Balance the integrated vol flux = 0
%
% Input: un - 2D normal velocity component along contour
%        
%
Asrf=[];
Asrf=DZ.*DX;
Atot=nansum(nansum(Asrf));

vf0=un.*Asrf;
Vflx=nansum(nansum(vf0));
hf1=nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
Hflx1=nansum(hf1);
hf2=nansum(un.*Cp.*rhow.*(tn-Tref2).*DX.*DZ);
Hflx2=nansum(hf2);
fprintf('Adjust flux Fluxes before adj\n');
fprintf('=== VFlux=%6.1f Sv, HFlux1=%10.8d, HFlux2=%10.8d\n',...
	Vflx*1e-6, Hflx1, Hflx2);

% Only near-bottom layers:
nbt=5;
[nlr,npt]=size(DX);
Msk=DX*nan;
for ii=1:npt
  dz=DZ(:,ii);
  imn=max(find(dz>1e-10));
  if isempty(imn), continue; end;
  ib=imn-nbt;
  Msk(ib:ib+nbt,ii)=1; % bottom layers
  Msk(1:ib-1,ii)=-1;
end

Ibt=find(Msk==1);
Itp=find(Msk==-1);
Itot=find(~isnan(Msk));

ubt=un(Ibt);
utp=un(Itp);
vfbt=vf0(Ibt);
vftp=vf0(Itp);
dv=sum(vftp)+sum(vfbt); % 

% distribute excess flux over bottom layers:
% weighted by U btm
% Weight adjusted volume by 
% vol flux at bottom layers
% to give higher correction for
% strong U
Abt=Asrf(Ibt); % surface area of bottom cells
wt=vfbt./sum(vfbt);
dVflx=-dv*wt;
dU=dVflx./Abt; % vel correction m/s

du0=-dv/Abt;  % correction without Uflx weighting

ubt_adj=ubt+dU;
Uadj=un;
Uadj(Ibt)=ubt_adj;

% Double check:
Vfa=Uadj.*Asrf;
Vflx=nansum(nansum(Vfa));
hf1=nansum(Uadj.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
Hflx1=nansum(hf1);
hf2=nansum(Uadj.*Cp.*rhow.*(tn-Tref2).*DX.*DZ);
Hflx2=nansum(hf2);
fprintf('Adjusted: \n');
fprintf('=== VFlux=%6.1f Sv, HFlux1=%10.8d, HFlux2=%10.8d\n',...
	Vflx*1e-6, Hflx1, Hflx2);


return