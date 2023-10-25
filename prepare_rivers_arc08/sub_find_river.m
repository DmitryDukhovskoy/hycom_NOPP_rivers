function RR = sub_find_river(rname,rivname,Flow,LAT,LON);
% Find river rname 
% in NCAR data - unproccessed data
%
a1=length(LAT);
lnm = length(rname);
RR = [];
for ik=1:a1
  rnm = deblank(rivname(ik,:));
  if strncmp(rnm,rname,lnm)
    Q = Flow(:,ik);
    RR.Riv_name = rnm;
    RR.Indx = ik;
    RR.Lon = LON(ik);
    RR.Lat = LAT(ik);
    RR.Q = Q;
    break;
  end
  
end

if isempty(RR),
  fprintf('%s not found\n');
end


return
