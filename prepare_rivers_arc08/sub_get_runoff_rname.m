function [TM,QRv] = sub_get_runoff_rname(rnm);
% Get runoff for specific NCAR rivers - proccessed data (filled gaps, etc)
pthmat = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
%fmat   = sprintf('%sncar_rivers_Arctic_1993-2015.mat',pthmat);
fmat   = sprintf('%sncar_rivers_Arctic_1993-2016.mat',pthmat);
load(fmat); % <-- RVR

lnn = length(rnm);

TM = RVR.TM;
Rindx = RVR.Riv_indx;
Nrv = length(Rindx);
QRv= [];
for ir=1:Nrv
  rname = deblank(RVR.Riv_name(ir,:));
  if ~strncmp(rnm,rname,lnn); continue; end;
  
  QRv = RVR.Qflow_m3_sec(:,ir);
  if strncmp(rnm,'Mackenzie',9)
    QRv = 3*QRv; % in NCAR river data I splitted total runof in 3 branches
  elseif strncmp(rnm,'Lena',9)
    QRv = 3*QRv; % in NCAR river data I splitted total runof in 3 branches
  end
end;

if isempty(QRv)
  fprintf('%s not found \n',rnm);
  QRv = nan;
end


return