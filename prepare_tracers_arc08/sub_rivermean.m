function [Rmean,JR] = sub_rivermean(riv_fid,GRMSK);
%
% Average river runoff for the given year
% Find seasonally disappearing rivers (freeze up)
% Rmean - annual mean for all river sources
% JR - sources that freeze up in winter

[JDM,IDM]=size(GRMSK);
IJDM = IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

disp('Caclulating annual mean ...');
frewind(riv_fid);
Rmean= zeros(JDM,IDM);
NP   = zeros(JDM,IDM);

for k=1:12
  A  = fread(riv_fid,IJDM,'float32'); % read 2D field
  dm1= fread(riv_fid,npad,'float32');  % Padding = size(toto)
  if size(dm1) ~= size(toto)
    error('Padding in HYCOM file ???');
  end
  I  = find(A>1e10);
  A(I)= NaN;
  A=reshape(A,IDM,JDM)';
  A=A.*GRMSK;

  Rmean=Rmean+A;
  I=find(A>0);
  NP(I)=NP(I)+1;
  ntot(k)=length(I); % # of rivers by month
end

Rmean=Rmean/k;
dmm=Rmean;
dmm(NP==k)=nan;
dmm(dmm==0)=nan;
JR=find(~isnan(dmm));
%keyboard

return;