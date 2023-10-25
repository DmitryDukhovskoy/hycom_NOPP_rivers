% Calculate volume transport across straits/sections
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;

regn = 'ARCc0.08';
expt = '110';
TV   = 11;  % topo version
%segm = 'BeringS';
 segm = 'FramS';
%segm = 'BarOp';
%segm = 'DavisS';
YR1  = 1993;
YR2  = 2016;

%Zmn = 100; % average over the top Zmn m
rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%s/fig_sections/';
pthfgG  = '/Net/ocean/ddmitry/hycom/ARCc0.08/110/fig_GrSct/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);

fmat = sprintf('%s%s_Vflux_%s_%4.4i_%4.4i.mat',...
	       pthmat,expt,segm,YR1,YR2);

fprintf('Loading %s\n',fmat);
load(fmat);
nrc = length(FLXV); 

d1=TM(1);
dv=datevec(d1);
DV = datevec(TM);
yr1=dv(1);
dj1=datenum(yr1,1,1);

nrc1=nrc;
DV1=DV;

dn  = 365/12;
dTM = TM-dj1;
dM  = dTM/dn+1;
Vmn = nanmean(FLXV);

% Group by years
cc = 0;
yr_old = 0;
Iyr = [];
for j=1:nrc
  yr = DV(j,1);
  if yr~=yr_old
    cc = cc+1;
    yr_old = yr;
    Iyr(cc,1) = j;
  end
end
nyrs=cc;

% Calculate annual means
nday_yr = floor(nrc/nyrs);
Fyr=[];
Fstdv=[];
F2d = [];
for jj=1:nyrs-1
  j1 = Iyr(jj);
  j2 = Iyr(jj+1);
  A = FLXV(j1:j2);
  Fyr(jj,1) = nanmean(A);
  Fstdv(jj,1) = std(A);
  F2d(1:nday_yr,jj) = A(1:nday_yr);
end
A = FLXV(Iyr(end):end);
Fyr(end+1,1) = nanmean(A);
Fstdv(end+1,1) = std(A);
F2d(1:nday_yr,jj+1) = A(1:nday_yr);

xyr = (TM-TM(1))/365.25+DV(1,1);

mF = mean(FLXV);
p10 = prctile(FLXV,10);
p90 = prctile(FLXV,90);
jYR = DV(Iyr,1);

figure(1); clf;
axes('Position',[0.08 0.55 0.87 0.35]);
%plot(dTM,FLXV,'linewidth',2);
plot(xyr,FLXV,'linewidth',2); 
hold on;
%plot([xyr(1) xyr(end)],[0 0],'k-');
plot([xyr(1) xyr(end)],[mF mF],'r-');
plot([xyr(1) xyr(end)],[p10 p10],'r--');
plot([xyr(1) xyr(end)],[p90 p90],'r--');
set(gca,'xlim',[xyr(1) ceil(xyr(end))],...
	'tickdir','out',...
	'xtick',[xyr(1):1:2020],...
	'xgrid','on',...
	'ygrid','on');

stt=sprintf('%s-%s, Vol Flux (Sv), %s, %i-%i, Vmn=%4.1f Sv',...
	    regn,expt,segm,YR1,YR2,Vmn);
title(stt);
txb='plot_vol_transp.m';

axes('Position',[0.08 0.08 0.87 0.35]);
boxplot(F2d,jYR);
stt=sprintf('%s-%s, Vol.Flux, %s, %i-%i, Vmn=%4.1f Sv',...
	    regn,expt,segm,YR1,YR2,Vmn);
title(stt);
set(gca,'tickdir','out',...
	'ylim',[min(FLXV)-0.2 max(FLXV)+0.2],...
	'yminortick','on')

bottom_text(txb,'pwd',1);


% ---------------------------------
% Monthly climatology
% ---------------------------------
for im=1:12
  Im=find(DV1(:,2)==im);
  dmm=FLXV(Im);
  prcL(im)=prctile(dmm,25);
  prcU(im)=prctile(dmm,75);
  Fmn(im)=nanmean(dmm);
end

mnth=[1.5:12.5];
clr=[0 0.4 0.7];
dw=0.4;
dh=0.18;
POS=[0.07 0.80 dw dh; ...
     0.07 0.55 dw dh; ...
     0.07 0.30 dw dh; ...
     0.07 0.05 dw dh];

% Save for plotting
%ftmp=sprintf('%shycom_vol.mat',pthmat);
%save(ftmp,'Fmn','prcL','prcU','mnF','stdv');


figure(8); clf;
axes('position',POS(2,:));
plot(mnth,Fmn,'Linewidth',2.5,'Color',clr);
hold on;
for ik=1:12
  plot([ik+0.5 ik+0.5],[prcL(ik) prcU(ik)],'k--','Color',clr);
end

mnF=nanmean(FLXV);
stdv=nanstd(FLXV);
stxt=sprintf('%4.1f+/-%4.1f',mnF,stdv);
text(5,0.9*min(prcL),stxt);

set(gca,'tickdir','out',...
	'xlim',[1 12.8],...
	'xtick',[1:12],...
	'ylim',[1.05*min(prcL) max([0, max(prcL)])]);
stlL=sprintf('arc08-%s, Fram Str Vol Flux, Sv',expt);
title(stlL,'Interpreter','none');

bottom_text(txb,'pwd',1,'pos',[0.04 0.5 0.4 0.04]);

if s_fig==1
  fgnm=sprintf('%sarc08_%s_VolFlux_FramStr',...
	       pthfgG,expt);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end


% -------------------------
% Vol transports northward/southward
% ---------------------------------
% UV sections in straits
fmat = sprintf('%s%s_UV_straits_%4.4i_%4.4iv2.mat',...
	       pthmat,expt,YR1,YR2);
fprintf('Loading %s\n',fmat);
load(fmat);  % UV
nuv=length(UV);

for ik=1:nuv
  nm=UV(ik).Name;
  if strncmp(nm,'FramS',4),
    break;
  end
end

% Calculate positive flux (Nortward or Eastward)
TM=UV(ik).TM;
uv=UV(ik).UV_normal;
ZZ=UV(ik).ZZ;
dx=UV(ik).Dist;
nrc=length(TM);
DV=datevec(TM);
xyr = (TM-TM(1))/365.25+DV(1,1);

nlv=41;
[DX,dmm]=meshgrid(dx,(1:nlv));

fprintf('Calculating Northward/southward transports ...\n');
for it=1:nrc
  zz=squeeze(ZZ(it,:,:));
  dz=diff(zz,1,1);
  agrd=abs(dz.*DX);
  v=squeeze(uv(it,:,:));
  Ip = find(v>0);
  In = find(v<0);
  Fp(it) = nansum(v(Ip).*agrd(Ip));
  Fn(it) = nansum(v(In).*agrd(In));
end
Fp = Fp*1e-6; %Sv
Fn = Fn*1e-6; 

mFp = mean(Fp);
mFn = mean(Fn);


% Annual means
% bar diagrams
% Group by years
cc = 0;
yr_old = 0;
Iyr = [];
for j=1:nrc
  yr = DV(j,1);
  if yr~=yr_old
    cc = cc+1;
    yr_old = yr;
    Iyr(cc,1) = j;
  end
end
nyrs=cc;

NDyr = diff(Iyr);
NDyr(cc+1) = nrc-Iyr(end); 

FPyr=[]; % positive fluxes
FP10=[];
FP90=[];
FNyr=[]; % negative fluxes
FP10=[];
FP90=[];
for jj=1:nyrs
  j1 = Iyr(jj);
  if jj==nyrs, 
    j2 = nrc; 
  else
    j2 = Iyr(jj+1);
  end;
%  fprintf('jj=%i, j1=%i, j2=%i\n',jj,j1,j2);
  A = Fp(j1:j2);
  FPyr(jj,1) = nanmean(A);
  FP10(jj,1) = prctile(A,10);
  FP90(jj,1) = prctile(A,90);
  A = Fn(j1:j2);
  FNyr(jj,1) = nanmean(A);
  FN10(jj,1) = prctile(A,10);
  FN90(jj,1) = prctile(A,90);
end
jYR = DV(Iyr,1);



figure(2); clf;
axes('Position',[0.08 0.55 0.87 0.35]);
plot(xyr,Fp,'Color',[1 0.5 0]);
hold on;
plot([xyr(1) xyr(end)],[mFp mFp],'-','Color',[1 0. 0.9]);
%plot([xyr(1) xyr(end)],[p10 p10],'r--');
%plot([xyr(1) xyr(end)],[p90 p90],'r--');
plot(xyr,Fn,'-','Color',[0 0.5 0.8]);
plot([xyr(1) xyr(end)],[mFn mFn],'-','Color',[0 1 1]);

set(gca,'xlim',[xyr(1) ceil(xyr(end))],...
	'tickdir','out',...
	'ylim',[-32 32],...
	'ytick',[-30:10:30],...
	'xtick',[xyr(1):1:2020],...
	'xgrid','on',...
	'ygrid','on');

stt=sprintf('%s-%s, +/- Vol.Fluxes, %s, %i-%i, V+=%4.1f, V-=%4.1fSv',...
	    regn,expt,segm,YR1,YR2,mFp,mFn);
title(stt);

axes('Position',[0.08 0.09 0.87 0.35]);
hold on;
bb=bar(jYR,FPyr,0.9);
set(bb,'Facecolor',[1 0.5 0],'EdgeColor','none');
bb2=bar(jYR,FNyr,0.9);
set(bb2,'Facecolor',[0 0.5 0.8],'EdgeColor','none');
for jj=1:nyrs,
  y0=jYR(jj);
  p1=FP10(jj);
  p2=FP90(jj);
  plot([y0 y0],[p1 p2],'k');
  p1=FN10(jj);
  p2=FN90(jj);
  plot([y0 y0],[p1 p2],'k');
end
set(gca,'xlim',[jYR(1)-0.5 jYR(end)+0.5],...
	'tickdir','out',...
	'ylim',[-25 25],...
	'ytick',[-30:5:30],...
	'xtick',[jYR(1):1:2020]);

title('Annual & 10-90%');


txb='plot_vol_transp.m';
bottom_text(txb,'pwd',1);

% Cumulative vol flux across Fram/
ik=2;
it=483;
V=squeeze(UV(ik).UV_normal(it,:,:));
ZZ=UV(ik).ZZ;
dx=UV(ik).Dist;
zz=squeeze(ZZ(it,:,:));
dz=diff(zz,1,1);
agrd=abs(dz.*DX);
VFlx=V.*agrd;
a=nansum(VFlx);
vf=cumsum(a);

