function [YRp,VFWp]=sub_read_platov;
% FWC in the domain from G> Platov
pthin='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_mat/';
flnm=sprintf('%sGrFW_SbpGyre_Platov.dat',pthin);

%fid=fopen(flnm,'r');
A=load(flnm);
na=size(A,1);

for it=1:na
  dtt=A(it,1);
  yr=floor(dtt);
  fr=dtt-yr;
  dJ1=datenum(yr,1,1);
  ndyr=datenum(yr,12,31)-dJ1+1;
  yrday=floor(ndyr*fr)+1;
  dnmb=dJ1+yrday-1;
  TM(it,1)=dnmb;
end
DV=datevec(TM);


% Monthly means
y1=DV(1,1);
y2=DV(end,1);
cc=0;
for yr=y1:y2
  for im=1:12
  cc=cc+1;
  I=find(DV(:,1)==yr & DV(:,2)==im);
  dmm=A(I(end),2);
  TMy(cc,1)=TM(I(end));
  VFWp(im,yr-y1+1)=dmm;
  end
end
YRp=[y1:y2];

return
