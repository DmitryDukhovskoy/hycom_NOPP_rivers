% Find relation between amount of FW and tracer values
% for rivers and Berins strait
% to adjust all values
% Use regression models
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.08';
expt = 110;  

s_mat = 1;  %=0 - do not save; =1 - save; =2 - load
plr   = 1;  % layer to analyze
s_fig = 0;
s_reg = 0; %=0 - archm files are whole domain
         %=1 - archm files subsampled into nAtl. Region,
	 %      check fortran/hycom/extract_subdomain.F90
	 
nTr  = 1;   % tracer 
TrMn = 1e-10; % Threshold value to plot
rg   = 9806;
Sref = 40;  % reference S

fprintf('Tracer #: %i, Threhold value: %8.5d\n',nTr,TrMn);


pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);
fmat    = sprintf('%strac_fw_Tr%i_Lr%i',pthmat,nTr,plr);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

switch(nTr)
 case(1), % Greenland
   IJBX(1,1:2) = [565,780]; 
   IJBX(2,1:2) = [620,780]; 
   IJBX(3,1:2) = [620,675]; 
   IJBX(4,1:2) = [565,675]; 
 case(2), % Mackenzie R.
   IJBX(1,1:2) = [401,1655]; 
   IJBX(2,1:2) = [425,1655]; 
   IJBX(3,1:2) = [425,1630]; 
   IJBX(4,1:2) = [401,1630]; 
 case(4), % Western Eurasian Rivers
   IJBX(1,1:2) = [1385,1375]; 
   IJBX(2,1:2) = [1385,1070]; 
   IJBX(3,1:2) = [1580,1070]; 
   IJBX(4,1:2) = [1580,1375]; 
 case(5), % Bering Str.
   IJBX(1,1:2) = [622,1915]; 
   IJBX(2,1:2) = [684,1915]; 
   IJBX(3,1:2) = [684,1902]; 
   IJBX(4,1:2) = [622,1902]; 
end

YRPLT=[];
cc=0;
for iyr=2005:2005
  for idd=1:365
%for iyr=2008:2008
%  for idd=1:7:365
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

[X,Y] = meshgrid((1:nn),(1:mm));
dmm = inpolygon(X,Y,IJBX(:,1),IJBX(:,2));
IN = find(dmm==1 & HH<0);

if s_mat<2 
  cnc=0;
  clear TM FWC TRC
  for ip=1:np
    yr=YRPLT(ip,1);
    iday=YRPLT(ip,2);

    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)

  % Mean files - all variables are collocated
    if s_reg==0
      fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
      finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
    elseif s_reg==1
      fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12nAtl.a',pthbin,expt,yr,iday);
      finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12nAtl.b',pthbin,expt,yr,iday);
    end

    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end

    cnc=cnc+1;
    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    TM(cnc,1) = dnmb;
  %fprintf('Plotting: Tracer #%i V.Layer=%i\n\n',nTr,plr);
    fprintf('Tr# %i, sfig=%i, V.Lr=%i, %4.4i_%2.2i_%2.2i: %s\n',...
	    nTr,s_fig,plr,DV(1:3),fina);
  % in archm u_vel=utot; total velocity
  %  keyboard
    [UN,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',plr);
    [VN,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',plr);
    UN(UN>1e6)=nan;
    UN=squeeze(UN); 
    VN(VN>1e6)=nan;
    VN=squeeze(VN); 

   % Layer thickness:
    [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
    F=F./rg;
    F(F>1e10)=0;
    F(F<1e-2)=0;
    dP=squeeze(F); 

  % Tracer  
    [F,n,m,l] = read_hycom(fina,finb,'tracer','r_layer',plr,'r_tracer',nTr);
    F(F>1e6) = 0;
    TrN=squeeze(F);

  % Salinity
    [F,n,m,l] = read_hycom(fina,finb,'salin','r_layer',plr);
    F(F>1e6) = 0;
    S=squeeze(F);

    FWC(cnc,1) = sum((Sref-S(IN))./Sref.*dP(IN).*Acell(IN))*1e-9; % km3 FW content
    TRC(cnc,1) = sum(TrN(IN).*dP(IN).*Acell(IN))*1e-3; % tonn, Tracer content, mass

  %  keyboard

  end  % time loop

  FTR.TM       = TM;
  FTR.FWC_km3  = FWC;
  FTR.TRC_tonn = TRC;
  FTR.nTr      = nTr;
  FTR.Layer    = plr;
  FTR.Area_m2  = sum(Acell(IN));
  if s_mat == 1
    fprintf('Saving %s\n',fmat);
    save(fmat,'FTR');
  end
end
% 

fprintf('Loading %s\n',fmat);
load(fmat);

plr = FTR.Layer;
nTr = FTR.nTr;
% Normalize
AA = FTR.Area_m2;
dmm = FTR.FWC_km3*1e9; % km3->m3
FWC = dmm/AA; % m of water in 1 layer

dmm = FTR.TRC_tonn*1e3/1e3; % tonn -> m3, assuming rho=1000
TRC = dmm/AA;

nt = length(TRC);
X = [(ones(nt,1)),TRC];
[B,bint,r,rint,stats] = regress(FWC,X);
rF = X*B; 

% Annual mean:
cff = mean(FWC)/mean(TRC)

figure(1); clf;
axes('Position',[0.08 0.55 0.8 0.35]);
plot(FWC); hold;
plot(rF,'r');
title('FWC, m/m2, and regr fit (red)');

axes('Position',[0.08 0.08 0.8 0.35]);
plot(TRC);
title('Tracer Content, m/m2');

fprintf('Regression fit to FWC, Tracer: %i, Layer %i\n',nTr,plr);
fprintf('    FWC(m/m2) = %6.4f + %6.4f*Tracer(m/m2)\n',B);
fprintf('    Annual mean coeff: FWC = cff*TRC (m/m2) %6.2f\n',cff);

