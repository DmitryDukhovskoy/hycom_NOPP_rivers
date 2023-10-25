% Calculate EOF of the monthly mean field at Fram Strait
% use 3z interpolated fields
% see interp2z_UTS_Fram.m
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear

expt = 112;
TV   = 11;
sctnm= 'Fram';

YR1=2006;
YR2=2011;

rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
btx='anls_UTS_Fram.m';

btx='calc_eof_Fram.m';

% combine fields:
cc=0;
VV=struct;
for yr=YR1:YR2
  fmat_out = sprintf('%sarc08-%3.3i_UTSonZ_daily_%s_%4.4i.mat',...
                 pthmat,expt,sctnm,yr);
  fprintf('Loading %s\n',fmat_out);
  load(fmat_out);
% SCTZ 
% Compute monthly means
% For EOF, use only segments that go east (j=const)
  if ~exist('Jpr','var');
    JJ=SCTZ.J;
    dJ=diff(JJ);
    Jpr=find(dJ==0);
  end

  U=SCTZ.NVel;
  TM=SCTZ.TM;
  ndays=size(U,1);
  DV=datevec(TM);
% Monthly means:
  for im=1:12
    I=find(DV(:,2)==im);
    aa=U(I,:,:);
    aa=squeeze(nanmean(aa,1));
    cc=cc+1;
    VV(cc).V=aa(:,Jpr);
    TMe(cc)=datenum(yr,im,15);
  end
% Weekly means:  
%  di=7;
%  for ii=1:di:ndays-di
%    aa=U(ii:ii+di-1,:,:);
%    aa=squeeze(nanmean(aa,1));
%    cc=cc+1;
%    VV(cc).V=aa(:,Jpr);
%    TMe(cc)=TM(ii+round(di/2));
%  end 
end

% Prepare a matrix F:
% rows are "maps" (spatial points)
% columns are time series
% take only points that are not NaN
ll=length(VV);
F=[];
fprintf('Forming time-space data matrix ...\n');

for ii=1:ll
  if mod(ii,20)==0,
    fprintf('%5.1f%% completed ...\n',ii/ll*100);
  end
  A=VV(ii).V;
  [a1,a2]=size(A);
  A=reshape(A,a1*a2,1);
  I=find(~isnan(A));
  F(ii,:)=A(I);
end

% detrend time series:
F=detrend(F,1);

fprintf('Calculating SVD ...\n');
[Up,Sp,Vp] = svd(F);

% Squared diagonal values of Sp (e/values of SVD) are the e/values of covariance matrix F'*F
% e/vectors are in Vp
Lmb = diag(Sp).^2;

% Reconstruct fields:
Eof=struct;
for iof=1:4
  dmm=Vp(:,iof);
  A=VV(1).V;
  [a1,a2]=size(A);
  A=reshape(A,a1*a2,1);
  I=find(~isnan(A));
  J=find(isnan(A));
  E=A;
  E(I)=dmm;
  E=reshape(E,a1,a2);
  Eof(iof).E=E;

% Pr. Components:
  pc=F*Vp(:,iof);
  Eof(iof).PC=pc;
  Eof(iof).perc=Lmb(iof)/sum(Lmb)*100;
end;


fmat_out=sprintf('%sarc008_%3.3i_%s_eof_%i-%i.mat',pthmat,expt,sctnm,YR1,YR2);
fprintf('Saving mat file %s ...\n',fmat_out);
save(fmat_out,'Eof','Lmb','TMe');

%/home/ddmitry/codes/anls_mtlb_utils/hycom_gom04/eof_YucV.m

fprintf('All done \n');







