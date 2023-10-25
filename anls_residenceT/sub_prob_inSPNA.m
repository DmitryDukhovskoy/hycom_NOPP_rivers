% Calculate probability of part. being inside SPNA at least for N years after release
function PRB = sub_prob_inSPNA(PP,SLR,IGR,Iyr);

fprintf('Calculating probability ...\n');

TM  = PP(1).TM; 
TMd = TM-TM(1);
TMy = TMd./365.24;

XP = PP(1).XPcmb;
YP = PP(1).YPcmb;
ZL = PP(1).ZLcmb;

clear INdp

nlr   = length(SLR);
NDout = 90;
NDays = [];
for ilr=1:nlr
  lr   = SLR(ilr);
  IZ   = find(ZL==lr);
  npp  = length(IZ);
  dTMy = diff(TMy);
  dTMy(end+1)=dTMy(end);

  XPz=XP(IZ,:);
  YPz=YP(IZ,:);
  [a1,a2]=size(XPz);
  for ipp=1:npp
    X0=XPz(ipp,:);
    Y0=YPz(ipp,:);
    II=inpolygon(X0,Y0,IGR(:,1),IGR(:,2));
    dtm=dTMy.*II';

    i0=min(find(II==0));
    if isempty(i0), i0=a2; end;
    while i0<a2
      dmm=II;
      dmm(1:i0)=0;
      i2=min(find(dmm>0));
      ddays=sum(dTMy(i0:i2-1))*365;
      if isempty(i2) | ddays>NDout; break; end;
      dmm=II;
      dmm(1:i2)=1;
      i0=min(find(dmm==0));
      if isempty(i0), i0=a2; end;
    end
    NDays(ilr,ipp) = sum(dtm(1:i0));

  end
end

% Calculate probability of partlce staying within the SPNA for at least N years
% Probability of transit time (exclude points that never left the domain)
%Iyr = [1:27];
nnr = length(Iyr);
NPyr=[];
TRyr=[];
PRb =[];  % Prob T=N years
Dmax = sum(dTMy); % max possible transit time for particle, otherwise it staied in the domain
for jj=1:nnr
  yr0=Iyr(jj);
  for ilr=1:nlr
    aa=NDays(ilr,:);
    IP = find(aa>=yr0);
    NPyr(ilr,jj)=length(IP);
    IT = find(aa>=yr0 & aa<=Dmax);
    TRyr(ilr,jj)=length(IT);
    Ipr = find(aa>=yr0-0.5 & aa<yr0+0.5);
    PRb(ilr,jj)=length(Ipr);
  end
end

%keyboard

% Probability:
PRyr = NPyr./npp;
TRyr = TRyr./npp;
PRb  = PRb./npp;

PRB.PRyr = PRyr;
PRB.TRyr = TRyr;
PRB.TMy  = TMy;
PRB.TMd  = TMd;
PRB.PRb  = PRb;

% Calculate Mean Age and 95% 
% Althought the ages are not normally distributed,
% Estimate mean age and 95% CI:
% Central Limit thm, for large n
% mean(x) ~N(0,1)
%
npart=size(NDays,2);
fprintf('Initial # of particles is found = %i\n',npart);
%keyboard

% Mean and 95CI
alf=0.05;
for ilr=1:nlr
  npage=NDays(ilr,:);
  xmn   = mean(npage);
  sgm   = std(npage);
  zalf2 = abs(norminv(alf/2,0,1));
  ci_low = xmn-sgm/sqrt(npart)*zalf2;
  ci_up  = xmn+sgm/sqrt(npart)*zalf2;
  PRB.Age_mn(ilr) = xmn;
  PRB.CI95(ilr,:) = [ci_low,ci_up];
end


PRB.NPart_yrs=NDays;






return
