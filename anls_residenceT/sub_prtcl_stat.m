% Statistcs particles
function PR = sub_prtcl_stat(PP,HH,nlr,SLR,xlim1,xlim2,ylim1,ylim2);

% Count statistics
% of particles in the grid boxes
IL=find(HH>=0);

[mm,nn]=size(HH);

% Combined for all layers/years:
XP = PP(1).XPcmb;
YP = PP(1).YPcmb;
ZL = PP(1).ZLcmb;

fprintf('Calculating probability of Lagr partc. for layers ...\n');
for ilr=1:nlr
  fprintf('Layer %i\n',ilr);

  SMM = zeros(mm,nn);
  SMM(IL)=nan;

  lr=SLR(ilr);
  IZ=find(ZL==lr);

  xmm=XP(IZ,:);
  ymm=YP(IZ,:);
  [a1,a2]=size(xmm);

  for iz=1:a1
    x0=round(xmm(iz,:));
    y0=round(ymm(iz,:));
    x0=x0(:);
    y0=y0(:);
    In=find(~isnan(x0));
    x0=x0(In);
    y0=y0(In);
    II=sub2ind(size(HH),y0,x0);
    SMM(II)=SMM(II)+1;
  end

  sgmx=3;
  sgmy=sgmx;
  npnts=3*sgmx;
  ilm1=xlim1-10;
  ilm2=xlim2+10;
  jlm1=ylim1-10;
  jlm2=ylim2+10;
  SMM = sub_gauss_filter(SMM,sgmx,sgmy,npnts,ilm1,ilm2,jlm1,jlm2);
% Scale to have overall integral = 1
  SMM = SMM/(nansum(nansum(SMM)));

  PR(ilr).SMM=SMM;
end



return
