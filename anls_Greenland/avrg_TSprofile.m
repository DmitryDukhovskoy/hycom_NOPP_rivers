% Average over all years 
% monthly mean TS profiles
% along the contour
%
% Plot contours
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

%pfld = 'salin'; % field to extract
%pfld = 'temp';
YR1 = 1993;
YR2 = 2016;
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
regn = 'ARCc0.08';
s_mat = 0; % =0 - extract data, no save; =1 extract & save, =2 - load saved

btx = 'avrg_TSprofile.m';


fprintf('Averaging T&S, expt %3.3i, mean %i-%i\n',expt,YR1,YR2);

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';
fout=sprintf('%sarc08_%3.3i_monthlyClimTS_GreenlCntr.mat',pthmat,expt);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section


if s_mat==1
  cc=0;
  for iyr=YR1:YR2
    fmat = sprintf('%s%3.3i_Greenl_TScontour_%i.mat',...
		     pthmat,expt,iyr);
    if ~exist(fmat,'file');
      fprintf('Not found %s\n',fmat);
      continue;
    end

    fprintf('Loading %s\n',fmat);
    load(fmat);
    cc=cc+1;
    for imo=1:12
      DZ=TSGR(imo).DZ;
      DZ(isnan(DZ))=0;
      T=TSGR(imo).T;
      S=TSGR(imo).S;
% Extend all the way to end for avoiding nan in averaging
      [m,n]=size(S);
      for ik=1:n
	inn=max(find(~isnan(S(:,ik))));
	if inn==m, continue; end;
	S(inn+1:end,ik)=S(inn,ik);
	T(inn+1:end,ik)=T(inn,ik);
%
% Fill in mid-column nans caused 
% by 0-thickness iospycn. layers
	I=find(isnan(S(:,ik)));
	if ~isempty(I) & Hs(ik)<0;
	  fprintf('  Fixing Middle-column nan,imo=%2.2i ik=%i\n',imo,ik);
	  for ip=1:length(I)
	    ii0=I(ip);
	    S(ii0,ik)=S(ii0-1,ik);
	    T(ii0,ik)=T(ii0-1,ik);
	  end
%	  keyboard;
	end
	
      end

      
      
      if cc==1
	TSCNTR(imo).DZ=DZ;
	TSCNTR(imo).T=T.*DZ;
	TSCNTR(imo).S=S.*DZ;
      else
	TSCNTR(imo).DZ=TSCNTR(imo).DZ+DZ;
	TSCNTR(imo).T=TSCNTR(imo).T+T.*DZ;
	TSCNTR(imo).S=TSCNTR(imo).S+S.*DZ;
      end
    end

  end

%keyboard
  % Finish Averaging
  % For plotting, add zurface layer at zz=0
  for imo=1:12
    DZ = TSCNTR(imo).DZ;
    i0 = find(DZ<1e-6);
    DZ(i0) = nan;
    [m,n]=size(DZ);

    dmm = TSCNTR(imo).T./DZ;
    dmm=[dmm(1,:);dmm];
    T=dmm;

    dmm = TSCNTR(imo).S./DZ;
    dmm=[dmm(1,:);dmm];
    S=dmm;
%
% Extend all the way down
    for ik=1:n
      inn=max(find(~isnan(S(:,ik))));
      if inn==m, continue; end;
      S(inn+1:end,ik)=S(inn,ik);
      T(inn+1:end,ik)=T(inn,ik);
    end
    TSCNTR(imo).T = T;
    TSCNTR(imo).S = S;
    

    dmm = TSCNTR(imo).DZ;
    TSCNTR(imo).DZ = dmm./cc;
    dz=TSCNTR(imo).DZ;
    i0=find(isnan(dz));
    dz(i0)=0;
    zz=-cumsum(dz);
    zz=[zeros(1,n);zz];

%    [m,n]=size(S);
%    for ik=1:n
%      zb=zz(:,ik);
%%      dHb=abs(zb-Hs(ik));
%      dzz=abs(diff(zb));
%      Ib=find(dzz<1e-3)+1;
%      zz(Ib,ik)=nan;
%    end
    
    TSCNTR(imo).ZZ =zz;
  end

  TSCNTR(1).Info=sprintf('T/S Monthly Vertical Distr, Gr. Contour %s-%s',YR1,YR2);
  TSCNTR(1).Hbottom=GC.Hbottom;
  TSCNTR(1).Distance_km=GC.Distance_m*1e-3; 

  fprintf('Saving monthly clim.: %s\n',fout);
  save(fout,'TSCNTR');

else 
  fprintf('Loading %s\',fout);
  load(fout);
end

% Plot Sections:
imo=1;

dx = TSCNTR(1).Distance_km; % 
Hs = TSCNTR(1).Hbottom;
%plot(dx,Hs); % plot Bottom profile along contour


T  = TSCNTR(imo).T;
S  = TSCNTR(imo).S;
ZZ = TSCNTR(imo).ZZ;

nlr = size(ZZ,1);
ar = ZZ*0;
[Dst,dmb] = meshgrid(dx,[1:nlr]);

tstr=sprintf('arc008_%3.3i, T, Mo=%i, %i-%i',expt,imo,YR1,YR2);
% For plotting add surface layer      
%ZZp = [ZZ(1,:);ZZ];
%ZZp(1,:)=0;
%Tp  = [T(1,:); T];
%Dstp = [Dst(1,:);Dst];
c1=-2;
c2=10;
fnmb = 1;
pfld='temp';
sub_plot_TS_Zcntr(T,ZZ,Dst,fnmb,tstr,Hs,c1,c2,pfld);
btx='avrg_TSprofile.m';  
bottom_text(btx,'pwd',1,'position',[0.08 0.2 0.4 0.03]);


c1=31;
c2=35;
fnmb = 2;
tstr=sprintf('arc008_%3.3i, S, Mo=%i, %i-%i',expt,imo,YR1,YR2);
pfld='salin';
sub_plot_TS_Zcntr(S,ZZ,Dst,fnmb,tstr,Hs,c1,c2,pfld);

bottom_text(btx,'pwd',1,'position',[0.08 0.2 0.4 0.03]);

