% Thickness and velocity diffusivity
% 2D fields
% No land mask, so no need to change
%forfun.f:c --- veldf2 is diffusion velocity for laplacian  background diffusion
%forfun.f:c --- veldf4 is diffusion velocity for biharmonic background diffusion
%c --- 'thkdf2' = diffusion velocity (m/s) for Laplacian  thickness diffusion
%c --- 'thkdf4' = diffusion velocity (m/s) for biharmonic thickness diffusion 
%c ---             (negative to input spacially varying diffusion velocity)   
% diffusivity = thkdf2*dx or thkdf4*dx^3 

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.04';
E = '010';
%ntopo1=09;
%ntopo2=11;
ntopo2=17;
TV = sprintf('%2.2iDD',ntopo2);

ptharc    = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
pthin     = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/010/'; % relax 32lyrs
pthout    = sprintf('/Net/mars/ddmitry/hycom/%s/relax/%s/41layers_T%s/',R,E,TV);
pthtopo   = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo04 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);

% Old topo:
fltopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo);
HHo   = nc_varget(fltopo,'Bathymetry');
LATo  = nc_varget(fltopo,'Latitude');
%LONo = nc_varget(fltopo,'Longitude');


% New topo
fltopo_new=sprintf('%sdepth_%s_%s.nc',pthtopo04,R,TV);
HH   = nc_varget(fltopo_new,'Bathymetry');
LAT = nc_varget(fltopo_new,'Latitude');
%LON = nc_varget(fltopo_new,'Longitude');
%[mm,nn]= size(HH);

%
fina = sprintf('%sthkdf4.a',pthin);
finb = sprintf('%sthkdf4.b',pthin);
fouta = sprintf('%sthkdf4_T%s.a',pthout,TV);
foutb = sprintf('%sthkdf4_T%s.b',pthout,TV);

IDo = 1600;
JDo = 2520;
IJDo = IDo*JDo;

thkdf = 0.0100; 
fida = fopen(fina,'r','ieee-be');
%fidb = fopen(finb,'r');
%fclose(fidb);
 
IDM=3200;
JDM=5040;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
%kold=32;
knew=41;

fprintf('%s domain, ID=%i JD=%i\n',R,IDM,JDM);

% thkdf4
% Write *a file:
fida = fopen(fouta,'w');
A = ones(JDM,IDM)*thkdf;
A = reshape(A',IJDM,1);
fwrite(fida,A,'float32','ieee-be');
fwrite(fida,toto,'float32','ieee-be');  % padding at the end
fclose(fida);

% Write *b file:
fidb = fopen(foutb,'wt');
astr = sprintf('thkdf4: range =    %6.4f    %6.4f\n',min(A),max(A));
fprintf(fidb,'%s',astr);
fclose(fidb);
fprintf('Written files: %s\n',fida);
fprintf('Written files: %s\n',fidb);


% -----------
% veldf2
%
% veldf2 - function of latitude
% with max in the N. Pole
% there is also a patch in Hudson Bay
% of max diff.
% -------------
read_old = 0;
if read_old>0
  fina_in = sprintf('%sveldf2.a',pthin);
  finb_in = sprintf('%sveldf2.b',pthin);
  fida_in = fopen(fina_in,'r','ieee-be');
  Ao = fread(fida_in,IJDo,'float32','ieee-be');
  Ao = reshape(Ao,IDo,JDo)';
  fclose(fida_in);
  
  Y = Ao(1:end,600:1200);
  X = LATo(1:end,600:1200);
  [a1,a2]=size(X);
  Y=reshape(Y,a1*a2,1);
  X=reshape(X,a1*a2,1);
  X=[ones(a1*a2,1),X];
  B = regress(Y,X);
  alf0 = B(1);
  alf1 = B(2);
  
% Hudson Bay:
  Y = Ao(320:900,80:250);
  X = LATo(320:900,80:250);
  [a1,a2]=size(X);
  Y=reshape(Y,a1*a2,1);
  X=reshape(X,a1*a2,1);
  X=[ones(a1*a2,1),X];
  B = regress(Y,X);
  bt0 = B(1);
  bt1 = B(2);

else
  alf0 = 0.000353678614189007;
  alf1 = 6.38479510885239e-05;
  bt0  = -0.000150671584057181;
  bt1  = 0.000107099157535391;
end

% COnstruct veldf2 field
A = alf0+alf1*LAT;
% Hudson Bay:
i1 = 105;
i2 = 540;
j1 = 567;
j2 = 1907;
X = LAT(j1:j2,i1:i2);
Y = bt0+bt1*X;
A(j1:j2,i1:i2)=Y; % path in Hudson Bay

% Filter edges:
fprintf('Filtering Hudson Bay patch, Gaussian Fltr...\n');
AF=A;
sgmx=4;
sgmy=4;
nn=5*sgmx;
for i=i1-nn:i2+nn
  for j=j1-nn:j2+nn
    if i>i1+nn & i<i2-nn & ...
       j>j1+nn & j<j2-nn 
      continue  % skip interior
    end
    
    x0=i;
    y0=j;
    [X,Y]=meshgrid([i-nn:i+nn],[j-nn:j+nn]);
    ww=1/(sgmx*sgmy*2*pi)*exp(-((X-x0).^2/(2*sgmx^2)+(Y-y0).^2/(2*sgmy^2)));
    yy=A(j-nn:j+nn,i-nn:i+nn);
    dmm=sum(sum(ww.*yy));
    AF(j,i)=dmm;
  end
end

A = AF;

f_plt=0;
if f_plt>0
  figure(1); clf;
  pcolor(A); shading flat
  axis('equal')
  hold; contour(HH,[0 0],'k');  
  caxis([0 7e-3]);
  colorbar
  title('ARCc0.04, veldf2','Fontsize',18)
end

fout_vel2a = sprintf('%sveldf2_T%s.a',pthout,TV);
fout_vel2b = sprintf('%sveldf2_T%s.b',pthout,TV);
% Write *a file:
fida = fopen(fout_vel2a,'w');
A = reshape(A',IJDM,1);
fwrite(fida,A,'float32','ieee-be');
fwrite(fida,toto,'float32','ieee-be');  % padding at the end
fclose(fida);

%2.86733941E-03 
% Write *b file:
fidb = fopen(fout_vel2b,'wt');
astr = sprintf('veldf2: range =   %14.8d  %14.8d\n',min(A),max(A));
fprintf(fidb,'%s',astr);
fclose(fidb);
fprintf('Written files: %s\n',fout_vel2a);
fprintf('Written files: %s\n',fout_vel2b);



% ------------------------
% veldf4 
% in Hudson Bay - veldf=0.01, over the rest = 0.02
% ------------------------
i1 = 105;
i2 = 540;
j1 = 567;
j2 = 1907;

A = ones(JDM,IDM)*0.02;
A(j1:j2,i1:i2) = 0.01;
% Gaussian Filter
fprintf(' veldf4, Filtering Hudson Bay patch, Gaussian Fltr...\n');
AF=A;
sgmx=4;
sgmy=4;
nn=5*sgmx;
for i=i1-nn:i2+nn
  for j=j1-nn:j2+nn
    if i>i1+nn & i<i2-nn & ...
       j>j1+nn & j<j2-nn 
      continue  % skip interior
    end
    
    x0=i;
    y0=j;
    [X,Y]=meshgrid([i-nn:i+nn],[j-nn:j+nn]);
    ww=1/(sgmx*sgmy*2*pi)*exp(-((X-x0).^2/(2*sgmx^2)+(Y-y0).^2/(2*sgmy^2)));
    yy=A(j-nn:j+nn,i-nn:i+nn);
    dmm=sum(sum(ww.*yy));
    AF(j,i)=dmm;
  end
end

A = AF;

% Vel4df
fout_vel4a = sprintf('%sveldf4_T%s.a',pthout,TV);
fout_vel4b = sprintf('%sveldf4_T%s.b',pthout,TV);
% Write *a file:
fida = fopen(fout_vel4a,'w');
A = reshape(A',IJDM,1);
fwrite(fida,A,'float32','ieee-be');
fwrite(fida,toto,'float32','ieee-be');  % padding at the end
fclose(fida);

% Write *b file:
fidb = fopen(fout_vel4b,'wt');
astr = sprintf('veldf4: range =    %6.4f    %6.4f\n',min(A),max(A));
fprintf(fidb,'%s',astr);
fclose(fidb);
fprintf('Written files: %s\n',fout_vel4a);
fprintf('Written files: %s\n',fout_vel4b);






