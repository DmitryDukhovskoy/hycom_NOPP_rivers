function ichck = sub_check_fld(Aold, Anew, iiS);
%Check if saved fields from consecutive 
% runs match at                                         
% last/first overlapping points
nrc_old = size(Aold,1);
ss = size(Anew);
dim = length(ss);
if iiS ~= nrc_old
  fprintf('Indcies do not match for checking flds, skipped...\n');
  ichk=2;
  return
end
if dim == 3
  a1=Aold(nrc_old,1,1);
  a2=Anew(1,1,1);
elseif dim == 2
  a1=Aold(nrc_old,1);
  a2=Anew(1,1);
else
  error('CHeck dimension: sub_check_fld\n');
end

if abs(a1-a2)<1e-7,
  ichck = 0; % good
else
  ichck = 1; % bad
end


return