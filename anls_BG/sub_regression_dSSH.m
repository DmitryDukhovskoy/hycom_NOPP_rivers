function sub_regression_dSSH(fmat);
% See if there is a good relationship
% btw wind stress curl and BG intensity
%
fprintf('Loading %s\n',fmat);
load(fmat);

cf=1e9;
TM = LAPLE.TM;
dE = LAPLE.d2E_BG;

DV=datevec(TM);
nc=length(TM);
TT=(TM-TM(1))/365.25+DV(1,1);

% Annual means:
% Year is Jan - Dec
cp=0;
clear dEY
for iyr=DV(1,1):DV(end,1),
  cp=cp+1;
  IY=find(DV(:,1)==iyr);
  dEY(cp,1)=mean(dE(IY));
end


% Winter d2E:
IS=find(DV(:,2)>4 & DV(:,2)<10);
IW=find(DV(:,2)<=4 | DV(:,2)>9);

clear TY dEW dES
cp=0;
for iyr=DV(1,1):DV(end,1),
% Winter: Oct - March (previous/current year)
%  IYW=find(DV(:,1)==iyr & (DV(:,2)<4 | DV(:,2)>9));
  IYW=find((DV(:,1)==iyr-1 & DV(:,2)>9) |...
	   (DV(:,1)==iyr & DV(:,2)<4));
  dmm=mean(dE(IYW));
  vmm=mean(dE(IYW));
  cp=cp+1;
  dEW(cp,1)=dmm; % winter curl by years
  TY(cp,1)=iyr;
  
% Summer
  IYS=find(DV(:,1)==iyr & (DV(:,2)>4 & DV(:,2)<10));
  smm=mean(dE(IYS));
  dES(cp,1)=smm;
end



% Wind curl - see anls_atms/calc_wndcrl.m
pthmat  ='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/data_mat/';
fcrl = sprintf('%swnd_curl_Arctic_cfsr.mat',pthmat);
fprintf('Loading %s\n',fcrl);
load(fcrl);

TMc = WCRL.TM;
TYc = WCRL.Wind_Yrs;
CY  = WCRL.Wind_Annual; % annual curl
CW  = WCRL.Wind_Curl_wnt;
CS  = WCRL.Wind_Curl_smr;


nc=length(TY);
% Autoregressive multivariate model
clear X
cc=1;
X=ones(nc-1,1);

%x=dEY(1:nc-1);% previous year dSSH
%%x=x-mean(x);
%cc=cc+1;
%X(:,cc)=x; 

x=CW(2:nc); % current winter dSSH
%x=x-mean(x);
cc=cc+1;
X(:,cc)=x;

x=CW(1:nc-1); % previous winter
%x=x-mean(x);
cc=cc+1;
X(:,cc)=x;

x=CS(2:nc); % current summer
%x=x-mean(x);
cc=cc+1;
X(:,cc)=x;

x=CS(1:nc-1); % previous summer
%x=x-mean(x);
cc=cc+1;
X(:,cc)=x;

y = dEY(2:nc);
%y=y-mean(y);
Y=y;
[B,Bint,R,Rint,stat] = regress(Y,X);
yf = X*B;

figure(1); clf;
axes('Position',[0.08 0.42 0.83 0.5]);
plot(TY(2:nc),Y);
hold;
plot(TY(2:nc),yf);
pl=legend('d2E','Yfit');
set(gca,'tickdir','out',...
	'xlim',[1994 2015],...
	'xtick',[1994:2015],...
	'ygrid','on',...
	'xgrid','on');

stl=sprintf('Multivar. a/regr., Lapl.SSH, R^2=%3.2f',stat(1));
title(stl);
btx='sub_regression_dSSH.m';
bottom_text(btx,'pwd',1,'Position',[0.02 0.3 0.8 0.2]);



return
