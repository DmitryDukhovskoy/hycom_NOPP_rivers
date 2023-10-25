% Plot relative Vorticity/f of the bottom layer and SSH contours
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = '110';
TV = 11;
rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 

sfig = 0;
zz0   = -100;  %  depth of calculation
%txtb = 'plot_deepU_ssh.m';
btxt = 'plot_vort008.m';

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/fig_2D/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Fcor = 2*omg*sind(LAT);
Lmsk = HH*0;
Lmsk(HH<zz0) = 1;

cnt = 0;
VRT = zeros(mm,nn);
Umn = zeros(mm,nn);
Vmn = zeros(mm,nn);
for yr=2006:2006
  fmat = sprintf('%smean_vort008_%i.mat',pthmat,yr);
  for iday=1:365
    tic;
    dnmb=datenum(yr,1,1)+iday-1;
%dnmb = datenum(2005,8,21); % date to plot
    DV =datevec(dnmb);
    yr = DV(1);
%iday = dnmb-datenum(yr,1,1)+1;
    pthbin = sprintf('/nexsan/archive/ARCc0.08_%s/data/%4.4i/',expt,yr);

    fina = sprintf('%s%s_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%s_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file') | ~exist(finb,'file');
      fprintf('Not found *a or *b: %s\n',fina);
      fprintf('                     %s\n',finb);
    %  continue;
    end

    %cnc=cnc+1;
    %TM(cnc,1)=dnmb;

    fprintf('Getting data expt %s: %4.4i_%2.2i_%2.2i: %s\n',expt,DV(1:3),fina);

    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    %  F=F(:,jnc1:jnc2,inc1:inc2);
%    U1=squeeze(F(1,:,:));
    F(F>huge)=0;
    UN=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
%    V1=squeeze(F(1,:,:));
    F(F>huge)=0;
    VN=F;

    if ~exist('ZZ','var')
      [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
      ZZ(isnan(ZZ))=100;
      ZM(isnan(ZM))=100;

      Iz=find(HH<zz0);
      for k=1:nlr
	zav(k) = nanmean(ZM(k,Iz));
      end
      dzav = abs(zav-zz0);
      Iz0 = find(dzav==min(dzav));
    end
    
    % Calculate dv/dx
    fprintf('Calculating dv/dx & du/dy\n');
    dVdX = HH*nan;
    dUdY = HH*nan;
    for j=1:m
      dvp = squeeze(VN(Iz0,j,2:end));
      dv0 = squeeze(VN(Iz0,j,1:end-1));
      dx  = DX(j,1:end-1);
      dv0=dv0(:)';
      dvp=dvp(:)';
      dx=dx(:)';
      dVdX(j,1:end-1) = (dvp-dv0)./dx;
    end
%    dVdX(:,n)=nan;
    for i=1:n
      dup = squeeze(UN(Iz0,2:end,i));
      du0 = squeeze(UN(Iz0,1:end-1,i));
      dy  = DY(1:end-1,i);
      dup=dup(:);
      du0=du0(:);
      dUdY(1:end-1,i) = (dup-du0)./dy;
    end
%    dUdY(m,:)=nan;

    ZtF = (dVdX-dUdY)./Fcor;
    ZtF(HH>zz0)=nan;

    cnt=cnt+1;
    VRT = VRT+ZtF;
    Umn = Umn+squeeze(UN(Iz0,:,:));
    Vmn = Vmn+squeeze(VN(Iz0,:,:));
    
    fprintf('1 record: %6.4f min\n\n',toc/60);
  end
end

VRT = VRT./cnt;
Umn = Umn./cnt;
Vmn = Vmn./cnt;

MVRT.Depth    = zz0;
MVRT.Nrec     = cnt;
MVRT.Vrt_Fcor = VRT;
MVRT.U_mean   = Umn;
MVRT.V_mean   = Vmn;

fprintf('Saving %s\n',fmat);
save(fmat,'MVRT');


