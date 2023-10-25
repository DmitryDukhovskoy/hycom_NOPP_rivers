% Platov fluxes
% GFWA flux km3/day
% year day Davis Denmark


pthnm='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';
fnm=sprintf('%sdavis_denmark.txt',pthnm);

%fid=fopen(fnm,'r');
%frewind(fid);

A = load(fnm);

nrc=length(A);
for ii=1:nrc
  yr=A(ii,1);
  dd=A(ii,2);
  TM(ii,1)=datenum(yr,1,1)+dd-1;
end;

DV=[];
DM=[];
cc=0;
for iyr=2012:2015
  I=find(A(:,1)==iyr);
  cc=cc+1;
  DV(cc,1)=mean(A(I,3));
  DM(cc,1)=mean(A(I,4));
end



