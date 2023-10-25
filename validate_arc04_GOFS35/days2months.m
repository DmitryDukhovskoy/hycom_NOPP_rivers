% Transfer mat days to months
% for plotting
function MM = days2months(TM);
DV = datevec(TM);
icc = 0;
yold = DV(1,1);
mold = -1;
%keyboard
for ill=1:length(TM)
  yr = DV(ill,1);
  mo = DV(ill,2);
  dm = DV(ill,3);

  if yr~=yold
    icc=icc+1;
    yold = yr;
  end

  dmm = datenum(yr,mo,1)+32;
  dv2 = datevec(dmm);
  dmm = datenum(dv2(1),dv2(2),1);
  ndays = dmm - datenum(yr,mo,1);

  MM(ill,1) = icc*12+mo+dm/ndays;

end
  

return	
