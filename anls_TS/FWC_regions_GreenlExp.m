% Plot difference of 
% Extract S fields and calc mean for several regions
% N. Atlantic
% over specified days for given year 
% Monthly mean S is averaged within the layers
% and exactly within the specified depth layers
% for 2 experiments (with & without Greenland runoff)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

yr1=1993;
yr2=2016;


LRS = load('LRS.dat');
nlrs= length(LRS)-1; % skip whole depth  
       
regn = 'ARCc0.08';
expt = 110; % no Greenland runoff  
%expt = 112;  % Greenland runoff

s_fig = 0;
s_mat = 0;

% Experiments:
% 110 - no Greenland runoff
% 112 - Greenland runoff

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% N Atl regions
%BX = sub_define_boxes(HH,LON,LAT,1);
[XX,YY] = meshgrid((1:nn),(1:mm));

% Subpolar Gyre:
% Labr - Irm - eastern & central N. Atlantic
IGR = [  430         729
         463         677
         571         666
         810         621
         855         530
         930         501
        1031         445
        1099         410
        1192         403
        1228         321
        1218         185
        1080         145
         444         145
         416         246
         367         492
         374         591
         379         729];

INP = inpolygon(XX,YY,IGR(:,1),IGR(:,2));
INGr= find(INP==1 & HH<0);
II = find(HH<0);

%fss = sprintf('%sFWC_GreenlExp.mat',pthmat);
%fss = sprintf('%sFWC_regions_GreenlExp.mat',pthmat);
fss = sprintf('%sFWC_SubpolarG_GreenlExp.mat',pthmat);

if s_mat==1
  ic=0;
  for iyr=yr1:yr2;
    YR=iyr;
    cc=0;
    for mo=1:12
      ic=ic+1;
      cc=cc+1;
      expt1=110;
      pthmat1  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt1);
      pthm1 = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt1);
      fmat1 = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthm1,expt1,YR,mo);
      fprintf('Loading %s\n',fmat1);
      load(fmat1);

    % Calculate S content
    % Assuming S is in mg/m3 = 1e-3 kg/m3
% FWC = (Sold*Vol-Snew*Vol)/Snew    
      for ilv=1:4
	dz  = abs(LRS(ilv,2)-LRS(ilv,1));
	S   = meanS(ilv).Savrg;
	sms =  S*1e-3*dz.*Acell;
	if cc==1; S1(ilv).Smass=0; end;
	S1(ilv).Smass = S1(ilv).Smass+sms;
	S1(ilv).S = S;
      end  


      expt2=112;
      pthmat2  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt2);
      pthm2 = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt2);
      fmat2 = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthm2,expt2,YR,mo);
      fprintf('Loading %s\n',fmat2);
      load(fmat2);

      for ilv=1:4
	dz  = abs(LRS(ilv,2)-LRS(ilv,1));
	S   = meanS(ilv).Savrg;
	sms =  S*1e-3*dz.*Acell;
	if cc==1; S2(ilv).Smass=0; end;
	S2(ilv).Smass = S2(ilv).Smass+sms;
	S2(ilv).S = S;
      end  

% Estimatee FWC=(S1-S2)*V/S1
% Monthly FWC from monthly mean S
      VfwA=0;
      Vfw=0;
      for ilv=1:4
	dz = abs(LRS(ilv,2)-LRS(ilv,1));
	s1 = S1(ilv).S;
	s2 = S2(ilv).S; 
        dz = abs(LRS(ilv,2)-LRS(ilv,1));
        Vol= dz*Acell;
        dV = (s1-s2).*Vol./s1;
        VfwA= VfwA+nansum(dV(INGr));
        I=find(dV<0);
        dV(I)=0;
        Vfw= Vfw+nansum(dV(INGr));
      end
      VfwM(ic,1)=Vfw;
      VfwMA(ic,1)=VfwA;
      
    end; % month

%keyboard    
    if cc>1
      for ilv=1:nlrs
	dmm=S1(ilv).Smass;
	S1(ilv).Smass=dmm/cc;

	dmm=S2(ilv).Smass;
	S2(ilv).Smass=dmm/cc;
      end
    end

    nyr=iyr-yr1+1;
    fprintf('==== mFWC=%8.2f/%8.2f km3\n\n',...
	    mean(VfwMA(ic-11:ic))*1e-9, mean(VfwM(ic-11:ic))*1e-9);
%keyboard

  end; % year
  
  FWC.dFWC_500m_negS = VfwM; % only dS<0 
  FWC.dFWC_500m = VfwMA;
  
  fprintf('Saving %s\n',fss);
  save(fss,'FWC');
else
  fprintf('Loading %s\n',fss);
  load(fss);
end

% Calculate FWC change
dz=500;
Vol=nansum(nansum(Acell(INGr)))*dz;

clear dFWC
dFWC = FWC.dFWC_500m_negS*1e-9; % m3-> km3

% Plot overall FWC change 
yrs=[1993:2016];
yrm=[1993:1/12:2016.99];
figure(1); clf; 
axes('Position',[0.1 0.5 0.85 0.42]);
bar(yrm,dFWC,'FaceColor',[0.1 0.6 0.9]);
hold on;
for ii=yrs(1):yrs(end)
  ny=ii-yrs(1)+1;
  ic1=(ny-1)*12+1;
  ic2=ny*12;
  dfwm=mean(dFWC(ic1:ic2));
  plot([ii ii+0.99],[dfwm dfwm],'r-');
  plot([ii ii],[0 dfwm],'r:')
  plot([ii+1 ii+1],[0 dfwm],'r:')
end



set(gca,'tickdir','out',...
	'xtick',[1993:2016],...
	'Fontsize',14);
title('dlt FWC, km^3, 0-500m, Subpolar Gyre');
btx='FWC_regions_GreenlExp.m';
bottom_text(btx,'pwd',1,'Position',[0.03 0.3 0.4 0.05]);

	

f_br=0;
if f_br==1
% Integrate over depths:
  AA = squeeze(sum(dSM,2));
  
figure(2); clf;
axes('Position',[0.09 0.4 0.84 0.42]);
bar(AA');
set(gca,'xticklabel',{'Lbr','Irm','Icl','Grn','Bff'});
%hl=legend('50','150','300','500');
%set(hl,'Position',[0.71 0.27 0.06 0.12]);
ttl=sprintf('Salt Mass difference, expt112/110, upper 500 m');
title(ttl);

  
btx='FWC_GreenlExp.m';
bottom_text(btx,'pwd',1,'Position',[0.03 0.3 0.4 0.05]);
end


