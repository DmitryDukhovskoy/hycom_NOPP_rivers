% get GrSh particles
%
function PP = sub_GrSh(fnmout);

fprintf('Loading %s\n',fnmout);
load(fnmout);

nll=length(PP);
ntm=0;
for iip=1:nll
  tm=PP(iip).TM;
  nt=length(tm);
  ntm=max([ntm,nt]);
  if ntm==nt
    Li=iip;
  end
end
PP(1).ntm=ntm;

TM0=PP(Li).TM;
% Exclude repeated dates:
dTm=diff(TM0);
I=find(dTm>0);
I=[I;ntm];
TM0=TM0(I);
ntm=length(TM0);

% Fill in missing dates in the experiments 
% so that all experiments have same # of days (ntm)
clear nnm
nll=length(PP);
for iip=1:nll
  Xp=PP(iip).Xp;
  Yp=PP(iip).Yp;
  Tm=PP(iip).TM;
  zl=PP(iip).layer;
  [a1,a2]=size(Xp);

  Xn=[];
  Yn=[];
  for it=1:ntm
    t0=TM0(it);
    jt=find(Tm==t0,1);

    if ~isempty(jt),
      Xn(:,it)=Xp(:,jt);
      X0=Xp(:,jt);
      Yn(:,it)=Yp(:,jt);
      Y0=Yp(:,jt);
    else
      Xn(:,it)=X0;
      Yn(:,it)=Y0;
    end
  end

  PP(iip).Xp=Xn;
  PP(iip).Yp=Yn;
  PP(iip).ZL=ones(a1,1)*zl;
end


% Combine all experiments      
% Rows - all particles from nll runs for all depth layers
% Columns - Time, should be = TM0
XP=[];
YP=[];
ZL=[];
TM=TM0;

for iip=1:nll
  Xp=PP(iip).Xp;
  Yp=PP(iip).Yp;
  zl=PP(iip).ZL;
  XP=[XP;Xp];
  YP=[YP;Yp];
  ZL=[ZL;zl];
end;

PP(1).XPcmb=XP;
PP(1).YPcmb=YP;
PP(1).ZLcmb=ZL;
PP(1).TM=TM;

return
