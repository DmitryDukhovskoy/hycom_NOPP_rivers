function A = sub_iceA(fld,irc,ICE,A);

stl = sprintf('dmm = ICE(irc).%s;',fld);
eval(stl);
stl = sprintf('A.%s_mn(irc,1) = median(dmm);',fld);
eval(stl)
stl = sprintf('A.%s_iqr(irc,1) = prctile(dmm,25);',fld);
eval(stl)
stl = sprintf('A.%s_iqr(irc,2) = prctile(dmm,75);',fld);
eval(stl);

return
