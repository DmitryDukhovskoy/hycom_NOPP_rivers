% Thickness and velocity diffusivity
% 2D fields
% No land mask, so no need to change
%forfun.f:c --- veldf2 is diffusion velocity for laplacian  background diffusion
%forfun.f:c --- veldf4 is diffusion velocity for biharmonic background diffusion
%c --- 'thkdf2' = diffusion velocity (m/s) for Laplacian  thickness diffusion
%c --- 'thkdf4' = diffusion velocity (m/s) for biharmonic thickness diffusion 
%c ---             (negative to input spacially varying diffusion velocity)   
% diffusivity = thkdf2*dx or thkdf4*dx^3 
%
%veldf4 has already been increased to 0.02  in this region:
%
%narwhal04 3611> head *df4.b
%==> thkdf4.b <==
%thkdf4: range =    0.0100    0.0100
%
%==> veldf4.b <==
%veldf4: range =    0.0100    0.0200
%
%narwhal04 3613> hycom_range veldf4.a 3200 5040 1976 629
%min, max =   9.99999419E-03  1.99999996E-02    (1976, 629) =   1.99999996E-02
%
% You could try increasing veldf4 in the same patch to 0.04 (say) .  Note that there is a CFL limit on veldf4, so don't increase it too much.
%
% In both GLBc0.04 and ATLc0.02 we have increased thkdf4 as high as 0.07, and in ATLc0.02 this might be in the region where you are blowing up.
%
% If you do increase a df4 field, it is best to repeat the previous month to give it time to work.

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.04';
E = '023';
%ntopo1=09;
%ntopo2=11;
ntopo2=17;
TV = sprintf('%2.2iDD',ntopo2);

%ptharc    = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
%pthin     = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/010/'; % relax 32lyrs
pthout = '/Net/mars/ddmitry/hycom/ARCc0.04/relax/023/';
%pthtopo   = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo04 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);


% New topo
fltopo_new=sprintf('%sdepth_%s_%s.nc',pthtopo04,R,TV);
HH   = nc_varget(fltopo_new,'Bathymetry');
LAT = nc_varget(fltopo_new,'Latitude');
%LON = nc_varget(fltopo_new,'Longitude');
%[mm,nn]= size(HH);

IDM=3200;
JDM=5040;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
%kold=32;
knew=41;

fprintf('%s domain, ID=%i JD=%i\n',R,IDM,JDM);

% ------------------------
% veldf4 
%  = 0.04 near British Island = 0.06
% ------------------------
i1 = 1816;
i2 = 2220;
j1 = 225;
j2 = 880;


A = ones(JDM,IDM)*0.04;
A(j1:j2,i1:i2) = 0.06;
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
fout_vel4a = sprintf('%sveldf4_0406_T%s.a',pthout,TV);
fout_vel4b = sprintf('%sveldf4_0406_T%s.b',pthout,TV);
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






