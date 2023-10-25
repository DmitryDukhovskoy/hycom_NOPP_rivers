function zLdp = mean_ZM_41lrs(plr);
% Mean layer depths for deep ocean > 1000m
%
%
s_get_dz=0;

zLdp = 0;
if s_get_dz>0
  fprintf('Getting DZ ...\n');
  [ZM,ZZ] = sub_zz_zm(fina, finb,HH);
%  zL = squeeze(ZM(plr,:,:));
  Idp = find(HH<-1000);
  Ish = find(HH<-10 & HH>-200);
%  zLdp = nanmean(nanmean(zL(Idp)));
%  zLsh = nanmean(nanmean(zL(Ish)));
%  fprintf('Zldp=%6.1f, Zlsh=%6.1f\n',zLdp,zLsh);

  for kt=1:41
    zL = squeeze(ZM(kt,:,:));
    zLdp = nanmean(nanmean(zL(Idp)));
    zLsh = nanmean(nanmean(zL(Ish)));
    fprintf('%6.1f, %6.1f;\n',zLdp,zLsh);
  end      
end

A = [  -0.5,   -0.5;
  -1.9,   -1.8;
  -4.4,   -3.8;
  -8.4,   -6.8;
 -13.2,  -10.3;
 -18.6,  -14.3;
 -24.9,  -18.9;
 -32.3,  -24.2;
 -40.3,  -29.9;
 -48.3,  -35.6;
 -56.3,  -41.3;
 -64.3,  -47.0;
 -72.3,  -52.7;
 -80.3,  -58.4;
 -88.3,  -62.9;
 -96.3,  -66.1;
-104.3,  -69.0;
-112.3,  -71.5;
-120.3,  -73.8;
-128.3,  -75.9;
-136.3,  -77.7;
-144.5,  -79.3;
-153.6,  -80.6;
-166.6,  -81.8;
-191.1,  -83.1;
-237.4,  -84.0;
-311.8,  -84.3;
-405.2,  -84.3;
-514.6,  -84.3;
-630.6,  -84.3;
-752.8,  -84.3;
-888.0,  -84.3;
-1034.7,  -84.3;
-1256.2,  -84.3;
-1588.1,  -84.3;
-2031.1,  -84.3;
-2520.0,  -84.3;
-2872.8,  -84.3;
-3059.4,  -84.3;
-3203.5,  -84.3;
-3364.5,  -84.3;];

zLdp = A(plr,1);

return