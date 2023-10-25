  function [ip1,icnt,UN,VN] = sub_start_mat(PRTCL,YRPLT,expt,pthbin);
%
% Start Lagr. atlantic particle tracking from last record
% in mat file

dmm=PRTCL.TRACK;
icnt=length(dmm);
TM = PRTCL.TRACK(icnt).TM;
dv = datevec(TM(end));
dnmbYR = datenum(YRPLT(:,1),1,1)+YRPLT(:,2)-1;
ip1 = find(dnmbYR==TM);
ixy = min(find(YRPLT(:,1)==dv(1)));
yr   = YRPLT(ip1,1);
yrF  = dv(1);
iday = YRPLT(ip1,2);

dnmb=datenum(yrF,1,1)+iday-1;
dnmbF=dnmbYR(ip1);
if dnmb~=dnmbF,
  error('Check saved mat file time stamp and YRPLT');
end


dnmbUV=datenum(yr,1,1)+iday-1;
DV=datevec(dnmbUV);

fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

fprintf('Get old U and V: %4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

if ~exist(fina,'file');
  fprintf('Not found: %s\n\n',fina);
end

[F,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',1);
F(F>1e6)=nan;
UN=squeeze(F);
[F,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',1);
F(F>1e6)=nan;
VN=squeeze(F);


ip1=ip1+1;

yr   = YRPLT(ip1,1);
iday = YRPLT(ip1,2);

dnmbN=datenum(yr,1,1)+iday-1;
fprintf('RESTART: Last saved record: ip=%i\n',ip1-1);
fprintf('RESTART: Last saved record: %s\n',datestr(dnmb));
fprintf('RESTART: Will start from: record %i, %s\n\n',...
	ip1,datestr(dnmbN));


return


