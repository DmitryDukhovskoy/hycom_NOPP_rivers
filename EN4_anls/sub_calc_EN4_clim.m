function sub_calc_EN4_clim(YC1,YC2,fclim,pthdat,SCT);
% Calc long-term  monthly climatology 
% to calculate S anomalies
% from EN4 data
fprintf('Calculating EN4 climatology Years: %i - %i\n',YC1,YC2);

ns=length(SCT);

cc=0;
for ii=YC1:YC2
  for im=1:12
    cc=cc+1;
    YRPLT(cc,1)=ii;
    YRPLT(cc,2)=im;
  end
end

nrc=cc;

SMEAN = struct;
SMEAN = SCT;
SMEAN(1).Code='sub_cal_EN4_clim.m';
SMEAN(1).Info='S climatology for sections, EN4';
SMEAN(1).Years = [YC1,YC2];


cc=0;
for ik=1:nrc
  yr=YRPLT(ik,1);
  im=YRPLT(ik,2);
  fnm = sprintf('%sEN.4.1.1.f.analysis.g10.%4.4i%2.2i.nc',pthdat,yr,im);
  
  fprintf('S Clim: Reading %i/%i \n',yr,im);

  S = double(nc_varget(fnm,'salinity'));
  cc=cc+1;
  
  S = squeeze(S(1,:,:,:));


  for is=1:ns
    IJ=SCT(is).IJ;
    or=SCT(is).Orient;
    i1=IJ(1,1);
    j1=IJ(1,2);
    i2=IJ(2,1);
    j2=IJ(2,2);
    if or==1;
      ss=squeeze(S(:,j1:j2,i1));
    else
      ss=squeeze(S(:,j1,i1:i2));
    end
    
    if cc==1
      [a1,a2]=size(ss);
      dmm=zeros(12,a1,a2);
      SMEAN(is).S=dmm;
      SMEAN(is).cntr(1:12)=0;
    end
    
    smn=squeeze(SMEAN(is).S(im,:,:));
    smn=smn+ss;
    SMEAN(is).S(im,:,:)=smn;
    cntr=SMEAN(is).cntr(im)+1;
    SMEAN(is).cntr(im)=cntr;
  end
end

Nyrs=YC2-YC1+1;
if Nyrs~=cntr,
  fprintf('Check year counter: %i\n',Nyrs);
  keyboard;
end


for is=1:ns
%  for im=1:12
  smn=SMEAN(is).S;
%  cc=SMEAN(is).cc;
  S=smn./Nyrs;
  SMEAN(is).S=S;
end

fprintf('Saving climat S %s\n',fclim);
save(fclim,'SMEAN');

return