function H2P = sub_intrp_indx_piomas2hycom(s_mat);
% For inverse-distance interpolation
% of PIOMAS -> HYCOM ice fields:
% If needed (s_mat>0)
% Derive interpolation indices and weights only for ocean points
% within the area wher sea ice may be present
%
PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.data = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/PIOMAS_ice_data/';
PTH.rest = '/nexsan/people/ddmitry/hycom/ARCc0.08/103/ice_restart/';
PTH.indx = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/PIOMAS_ice_data/';

fmat=sprintf('%sremap_piomas2hycom.mat',PTH.indx);
fmat=sprintf('%sremap_piomas2arc008.mat',PTH.indx);

%keyboard
if s_mat==0
  fprintf('Loading Interp. Indices &Weights %s\n',fmat);
  load(fmat);
  return;
end



fgrds = sprintf('%sgrid.dat',PTH.data); % grid for scalar fields
nxl=360;
nyl=120;
nxy=nxl*nyl;

% read lon/lat scalar fields
dmm = load(fgrds);
[a1,a2]=size(dmm);
nrw=nxy/a2;

LONp = dmm(1:nrw,:);
LONp = reshape(LONp',[nxy,1]);
LONp = reshape(LONp,[nxl,nyl])';
[mp,np]=size(LONp);

LATp = dmm(nrw+1:end,:);
LATp = reshape(LATp',[nxy,1]);
LATp = reshape(LATp,[nxl,nyl])';

% HYCOM:
fltopo=sprintf('%sdepth_ARCc0.08_11.nc',PTH.topo);
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LATh  = alat;
LONh  = elon;
[m,n]= size(HH);
[mm,nn]= size(HH);

% Ice gridpoints,Mask:
dmm=HH;
dmm(1:200,:)=100;
dmm(1:350,630:end)=100;
dmm(1:250,1:500)=100;
dmm(1:500,950:end)=100;
dmm(1:800,1050:end)=100;
dmm(2400:end,:)=100;
dmm(2300:end,1:1000)=100;
dmm(1800:end,1:450)=100;

IMSK=find(dmm<0);
Np=length(IMSK);
nb=5; % # of nearest point ot search

H2P.HYCOM_IceLinIndx=IMSK;
H2P.HYCOM_size=[mm,nn];
H2P.PIOMAS_size=[mp,np];

i1=1;
if s_mat==2
  fprintf('Loading %s\n',fmat);
  load(fmat);
  dmm=H2P.Weights;
  i1=length(dmm)+1;
  fprintf('Found records=%i\n',i1-1);
  fprintf('Continue from i=%i\n',i1);
end


for ii=i1:Np
  
  if mod(ii,100)==0
    fprintf('ii=%i, out of %i, %6.2f\n',ii,Np,ii/Np*100);
  end
  
  li=IMSK(ii);
  xh=LONh(li);
  yh=LATh(li);
  
  d=distance_spheric_coord(yh,xh,LATp,LONp);
  for k=1:nb
    jn=find(d==min(min(d)),1);    
    IP(k,1)=jn;
    IP(k,2)=1./d(jn);
    d(jn)=1e20;
  end
  
  H2P.Weights(ii,1:nb)=IP(:,2)./sum(IP(:,2));
  H2P.PIOMAS_LinIndx(ii,1:nb)=IP(:,1);
  
  if mod(ii,1000)==0,
    fprintf('Saving %s\n',fmat);
    save(fmat,'H2P');
  end
  
end;

fprintf('Saving %s\n',fmat);
save(fmat,'H2P');





return
