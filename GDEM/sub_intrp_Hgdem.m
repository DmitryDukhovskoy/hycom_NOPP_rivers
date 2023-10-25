function Hgdem = sub_intrp_Hgdem(LAT,LON,HH,elon,alat);
% interpolate HYCOM topo
% into GDEM grid
% Quick nearest neighbor interpolation:
%
addpath /usr/people/ddmitry/codes/MyMatlab

elonM=elon;
I=elon>180;
elonM(I)=elonM(I)-360;
[XX,YY]=meshgrid(elonM,alat);

s = matlabpool('size');
if s==0
  matlabpool open 12;
end

[m,n]=size(XX);
nI=n*m;

Hgdem=zeros(m,n);
fprintf('Depth interpolation ... \n');
%cc=0;
%keyboard
parfor ii=1:nI
%  cc(ii)=1;
  x0=XX(ii);
  y0=YY(ii);
  D=distance_spheric_coord(y0,x0,LAT,LON);
  K0=find(D==min(min(D)),1);
  Hgdem(ii)=HH(K0);
%  dmm=sum(cc);
  if mod(ii,1000)==0
%    nfl=dmm/nI;
    fprintf('Interpolation i=%i \n',ii);
  end
  
end
 
matlabpool close;


return