% Derive Gr Runoff
% and save
function sub_get_GrRunoff(friv);

[TMg,Fgr]=sub_read_Greenland_v3; % Greenland total FWF, km3/mo
% Calculate annual anomalies
% and cumulative up to date
DVg = datevec(TMg);
ngr = length(TMg);
nyr =ngr/12;
dmm = reshape(Fgr,[12,nyr]);

Fyr = sum(dmm); % km3/yr
Ygr = [DVg(1,1):DVg(end,1)];
ii=find(Ygr==1990);
Fmn = mean(Fyr(1:ii));
ism = find(Ygr==dv0(1,1));
cFWF = cumsum(Fyr-Fmn);
fwf0 = cFWF(ism); % km3


% Greenland runoff with step-function increase of surplus runoff
% increase = mean (1993-2016) = 209 km3/yr
Fmn = 818.3;
dFmn = 209;
Fstep = Fyr;
iyr = find(Ygr==1993);
Fstep(iyr:end) = Fmn+dFmn; % Greenland runoff with step-function increase rate start on 1993
cFWF_step = cumsum(Fstep-Fmn);


save(frv,'cFWF','Fyr','Ygr','Fstep','cFWF_step');



return

