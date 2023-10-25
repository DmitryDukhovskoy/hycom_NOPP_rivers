function [ic,jc,Qr] = sub_get_loc_rivername(rnm,LT,LN,HH,year,mo);
% Find river location from NCAR rivers - proccessed data (filled gaps, etc)
% Also get runoff
pthmat = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
%fmat   = sprintf('%sncar_rivers_Arctic_1993-2015.mat',pthmat);
fmat   = sprintf('%sncar_rivers_Arctic_1993-2016.mat',pthmat);
load(fmat); % <-- RVR

lnn = length(rnm);
Iocn = find(HH<-0.1);

TM = RVR.TM;
itime = find(TM==datenum(year,mo,1));
Rindx = RVR.Riv_indx;
Nrv = length(Rindx);
Qr= [];
iC = 0;
jC = 0;
for ir=1:Nrv
  rname = deblank(RVR.Riv_name(ir,:));
  if ~strncmp(rnm,rname,lnn); continue; end;
  x0=RVR.Lon(ir);
  y0=RVR.Lat(ir);
  D=distance_spheric_coord(LT,LN,y0,x0);
  I0 = find(D==min(min(D)));
  [j0,i0] = ind2sub(size(HH),I0);

% Find closest ocean pnt, if needed
% If pnt is sea pnt - move to the coast first
  iL=0;
  if HH(j0,i0)<0
    D=distance_spheric_coord(LT(Iland),LN(Iland),y0,x0);
    imm = find(D==min(min(D)));
    iL = Iland(imm);
    [j0,i0] = ind2sub(size(HH),iL);
  end    
  
  D=distance_spheric_coord(LT(Iocn),LN(Iocn),y0,x0);
  imm = find(D==min(min(D)));
  iC = Iocn(imm);
  [jc,ic] = ind2sub(size(HH),iC);
  
  QRv = RVR.Qflow_m3_sec(:,ir);
  Qr = QRv(itime);
  
  if strncmp(rnm,'Mackenzie',9)
    Qr = 3*Qr; % in NCAR river data I splitted total runof in 3 branches
  elseif strncmp(rnm,'Lena',9)
    Qr = 3*Qr; % in NCAR river data I splitted total runof in 3 branches
  end
end;

if isempty(Qr)
  fprintf('%s not found \n',rnm);
  Qr = nan;
end


return