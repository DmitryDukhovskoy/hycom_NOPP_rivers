function Hen4 = sub_intrp_Hen4(LAT,LON,HH,elon,alat);
% interpolate HYCOM topo
% into EN4 grid for northern lat>50N
% Quick nearest neighbor interpolation:
% HH, LAT, LON - HYCOM topo, grid
% elon, alat - EN4 grid
addpath /usr/people/ddmitry/codes/MyMatlab

elonM=elon;
I=elon>180;
elonM(I)=elonM(I)-360;
[XX,YY]=meshgrid(elonM,alat);

%s = matlabpool('size');
%if s==0
%  matlabpool open 12;
%end

[m,n]=size(XX);
nI=n*m;

Hen4=zeros(m,n);
fprintf('Depth interpolation ... \n');
%cc=0;
%keyboard
for ii=1:nI
  if mod(ii,5000)==0
    fprintf('Topo interpolation, %6.2f%% done ...\n',ii/nI*100);
  end
  
%  cc(ii)=1;
  x0=XX(ii);
  y0=YY(ii);
  D=distance_spheric_coord(y0,x0,LAT,LON);
  K0=find(D==min(min(D)),1);
  Hen4(ii)=HH(K0);
%  dmm=sum(cc);
  
end
 
%matlabpool close;


return
