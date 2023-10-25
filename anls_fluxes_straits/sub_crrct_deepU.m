function Vadj = sub_crrct_deepU(un,tn,Tref1,Tref2,DX,DZ,Cp,rhow,vfmin,vfmax);
% Volume flux is not conserved
% in HYCOM mean files
% bogus deep currents
% up to 2.5 m/s
% adjust deep flow only
% correction is weighted by 
% the speed of the deep flow
% Adjust only bogus currents
%
% The integrated vol flux may not necessarily balanced = 0
%
% Input: un - 2D normal velocity component along contour
%        
%
Acell=[];
Acell=DZ.*DX;
Atot=nansum(nansum(Acell));

vf0=un.*Acell;
Vflx=nansum(nansum(vf0));
HF1=un.*Cp.*rhow.*(tn-Tref1).*Acell;
HF2=un.*Cp.*rhow.*(tn-Tref2).*Acell;
hf1=nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
Hflx1=nansum(hf1);
hf2=nansum(un.*Cp.*rhow.*(tn-Tref2).*DX.*DZ);
Hflx2=nansum(hf2);
fprintf('Adjust flux Fluxes before adj\n');
fprintf('===> VFlux=%6.1f Sv, HFlux1=%10.8d, HFlux2=%10.8d\n',...
	Vflx*1e-6, Hflx1, Hflx2);

% Only near-bottom layers:
nbt=5;
[nlr,npt]=size(DX);
Msk=DX*nan;
for ii=1:npt
  dz=DZ(:,ii);
  imn=min(find(dz<1e-10));
  if imn==1, continue; end;
  ib=imn-nbt;
  Msk(ib:ib+nbt-1,ii)=1; % bottom layers
  Msk(1:ib-1,ii)=-1;
  un(imn:end,ii)=nan; % get rid of "bubles" in the middles of nans
end

Ibt=find(Msk==1);
Itp=find(Msk==-1);
Itot=find(~isnan(Msk));

ubt=un(Ibt);
utp=un(Itp);
vfbt=vf0(Ibt);
vftp=vf0(Itp);
hfbt=HF2(Ibt);
hftp=HF2(Itp);
dv=sum(vftp)+sum(vfbt); % 

ubtmax=max(abs(ubt));

if (dv>vfmin | dv<vfmax) & ubtmax<0.25
  fprintf('No flux adjustment needed \n');
  Vadj=un;
  return;
end


% Find excessive fluxes
% Volume
if dv<0
  upp=prctile(ubt,10);
  uff=prctile(hfbt,10);
elseif dv>0
  upp=prctile(ubt,90);
  uff=prctile(HF2,90);
end  
dmm=un;
dmm(Msk==-1)=nan;
fmm=HF2;
fmm(Msk==-1)=nan;
Ixx=find(dmm<upp & fmm<=uff);

ubt0=un(Ixx);
vfbt0=vf0(Ixx);
hfbt0=HF2(Ixx);
dmm=vf0;
dmm(Ixx)=nan;
dv=sum(vfbt0)+nansum(nansum(dmm)); % 

%keyboard

% distribute excess flux over bottom layers:
% weighted by U btm
% Weight adjusted volume by 
% vol flux at bottom layers
% to give higher correction for
% strong U
Abt=Acell(Ixx); % surface area of adjusted bottom cells
%wt=vfbt0./sum(vfbt0); % weight based on vol flux
wt=hfbt0./sum(hfbt0);  % weight based on heat flux
dVflx=-dv*wt;
dU=dVflx./Abt; % vel correction m/s


ubt_adj=ubt0+dU;
Vadj=un;
Vadj(Ixx)=ubt_adj;

% Double check:
Vfa=Vadj.*Acell;
Vflx=nansum(nansum(Vfa));
hf1=nansum(Vadj.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
Hflx1=nansum(hf1);
hf2=nansum(Vadj.*Cp.*rhow.*(tn-Tref2).*DX.*DZ);
Hflx2=nansum(hf2);
fprintf('Adjusted Fluxes\n');
fprintf('<<= VFlux=%6.1f Sv, HFlux1=%10.8d, HFlux2=%10.8d\n',...
	Vflx*1e-6, Hflx1, Hflx2);


return