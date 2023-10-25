function iFmn = sub_indx_flx(Fmn,dp,cff,IIs,JJs,SCT);
% find indices of fluxes Fmn
% for points along the section 
% for plotting
npp = length(IIs);
Endp = SCT.End_Pnts;
iFmn=[];

%dp=10;
ip1 = min(find(~isnan(Fmn)));
%ngrd = 100; % max grid points to show flux
%cff = nanmax(abs(Fmn))/ngrd;
clc = 0;
for ipp=ip1:dp:npp-dp;
  x1=IIs(ipp);
  x2=IIs(ipp+dp);
  y1=JJs(ipp);
  y2=JJs(ipp+dp);
  if x1<Endp(2,1); % assuming section of 2 segments
    xS=Endp(1,1);
    xE=Endp(2,1);
    yS=Endp(1,2);
    yE=Endp(2,2);
  end
% Direction of the big segment
  aa  = (yE-yS)/(xE-xS);
  c   = yS-aa*xS;
  alf = atan2(aa,1); % angle of the section line
%
% Normal line:
% Note positive direction of 
% the flux: positive into Canada Basin
% so that if Flx>0, follow negative normal
  anm = -1/aa;
  cnm = yS-anm*xS;
  alf = atan2(anm,1); % angle of the normal line
  
  flx = nanmean(Fmn(ipp:ipp+dp))/cff;
%  if isnan(flx); continue; end;
  
  xps = 0.5*(x1+x2)-flx*cos(alf); % positive flux means negative normal dir. 
  yps = 0.5*(y1+y2)-flx*sin(alf);

  clc = clc+1;
  iFmn(clc,1)=xps;
  iFmn(clc,2)=yps;
%keyboard  
end

return