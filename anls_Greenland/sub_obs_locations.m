function POBS = sub_obs_locations;
%
%OK, Dmitry.
%Now you will see, why it was taking so long - just finished 
% putting the numbers together, so have not checked everything 
% there, but if there is a problem we must be able to spot it right away.
%The reason it was taking so long 
% (for the four regions with many points comprising their boundaries) 
% is obvious - I do not use circles or rectangles or preset regions 
% for area selection. Normally have certain criteria based on geography, 
% topography and even data density and monitoring line positions. 
% These are all programmed - a combination of conditions to be met, 
% including topography. So, what I did was plotting all points selected 
% in each domain and going around digitizing the boundaries 
% of the displayed data clusters by hand.
%Note there are point series and data syntheses. 
% In all data syntheses I breaks the signals into seasonal, 
% high-frequency and interannual variability and also clean 
% noise at several steps. We may not need all these details, 
% but may what is the sufficient minimum.
%I could have made polygon files in any format making 
% it easy for you to visualize (like bln or anything else), 
% but did not know what you prefer for that.  
%Cheers,
%Igor
i=1;
POBS(i).Name = 'Rockall Trough';
POBS(i).Lat   = 56.75;
POBS(i).Lon   = -11.0;

i=i+1;
POBS(i).Name = 'Iceland Basin';
POBS(i).Lat   = 60;
POBS(i).Lon   = -20;

i=i+1;
POBS(i).Name = 'North Irminger Sea';
POBS(i).Lat   = 64.2;
POBS(i).Lon   = -28;

i=i+1;
POBS(i).Name = 'Fylla Section';
POBS(i).Lat   = 63.88;
POBS(i).Lon   = -53.37;

i=i+1;
POBS(i).Name = 'SW Iceland Shelf';
POBS(i).Lat   = 63.00;
POBS(i).Lon   = -21.47;

i=i+1;
POBS(i).Name = 'OWS Mike';
POBS(i).Lat   = 66;
POBS(i).Lon   = 2;

i=i+1;
POBS(i).Name = 'Faroe Shetland Channel';
POBS(i).Lat   = [61, 61.5, 63, 61];
POBS(i).Lon   = [-3, -6, -6, -3];

i=i+1;
POBS(i).Name = 'North Labrador Sea';
POBS(i).Lat   = [58.9628
60.1117  
61.3245  
62.1223  
62.5372  
62.6649  
63.2713  
63.3351  
63.1117  
63.2394  
63.1436  
63.1117  
62.8245  
62.9202  
63.2713  
62.7926  
62.6968  
61.9309  
61.484  
61.2287  
60.7181  
60.75  
60.4628  
59.25  
59.0266  
59.0585  
59.5691  
59.8564  
59.7606  
59.633  
59.4096  
59.0904  
58.9947  
58.9628  
58.9628];
    
    
POBS(i).Lon   = -1*[59.44  
60.0260  
59.9787  
59.1990  
58.3719  
57.8521  
57.2377  
56.6233  
56.2925  
55.8436  
55.3474  
55.1347  
53.3861  
53.1734  
53.1262  
52.7009  
52.2519  
50.9995  
50.6215  
49.6999  
49.5818  
49.3691  
48.7784  
44.9031  
44.9976  
49.6999  
49.8417  
50.9050  
53.2680  
53.3861  
53.2207  
53.0789  
53.1734  
59.4589  
59.44];


i=i+1;
POBS(i).Name = 'Central Labrador Sea';
POBS(i).Lat   = [ 58.0465  
 58.3019  
 58.4056  
 58.4614  
 58.5891  
 58.7487  
 58.7646  
 58.6449  
 58.629  
 58.5572  
 58.3737  
 58.3816  
 58.1184  
 57.6875  
 57.129  
 56.7779  
 56.2593  
 56.1396  
 55.6928  
 55.5013  
 55.5093  
 55.5971  
 55.7008  
 56.012  
 56.1476  
 56.1077  
 56.1715  
 56.2832  
 56.4428  
 56.5465  
 56.6024  
 56.7939  
 56.9295  
 57.2566  
 57.3285  
 57.3045  
 57.3045  
 57.4721  
 57.6636  
 57.7513  
 58.0465];

POBS(i).Lon   = -1*[ 54.5506
55.0014  
55.0354  
54.9844  
54.1678  
54.0742  
53.8871  
53.7254  
52.577  
51.9646  
51.2925  
50.8842  
50.646  
50.4078  
50.1356  
50.0165  
49.7869  
49.6508  
49.6508  
49.957  
51.216  
51.7349  
52.0156  
52.3899  
52.4835  
52.7216  
53.1044  
53.3426  
53.4107  
53.4107  
53.2235  
53.2491  
53.2491  
53.2831  
53.4022  
53.5128  
53.6233  
53.5808  
53.6489  
53.87  
54.5506];

i=i+1;
POBS(i).Name = 'Central Irminger Sea';
POBS(i).Lat   = [   58.0053  
 58.8457  
 59.3883  
 59.5798  
 60.6649  
 61.0479  
 61.6436  
 61.633  
 61.25  
 61.2074  
 61.0053  
 60.6968  
 59.7926  
 58.4096  
 58.3351  
 58.0691  
 57.9947  
 57.9734  
 57.9947  
 58.016  
 58.0053];

POBS(i).Lon   = -1*[ 40.9518  
41.0274  
40.9518  
41.0274  
40.9802  
40.706  
40.1389  
39.0331  
37.7193  
37.3129  
36.5  
36.0274  
34.9688  
34.9972  
35.0444  
35.2713  
37.0104  
38.6834  
40.7344  
40.9329  
40.9518];

f_wLab=1;
if f_wLab==1
  i=i+1;
  POBS(i).Name = 'Whole Labrador Sea';
  POBS(i).Lat   = [ 54.0798  
 54.9734  
 55.5798  
 57.9734  
 60.1436  
 61.3564  
 62.4734  
 63.2074  
 63.1436  
 62.8883  
 62.7926  
 61.484  
 60.1436  
 59.6968  
 59.4096  
 59.1223  
 58.4202  
 58.1649  
 54.0798  
 54.0798  
 54.0798]; 

  POBS(i).Lon   = -1*[52.3937
52.4173  
56.5052  
58.9863  
60.026  
59.9551  
58.4192  
57.1905  
55.371  
53.0317  
52.2755  
50.6687  
48.0931  
46.6044  
46.439  
45.1158  
44.9031  
44.0288  
43.9579  
52.3464  
52.3937];

end

i=i+1;
POBS(i).Name = 'SW Iceland Shelf';
POBS(i).Lat   = 63.00;  
POBS(i).Lon   = -21.47;  

i=i+1;
POBS(i).Name = 'NorwSea Svinoy';
POBS(i).Lat   = 63.00;  
POBS(i).Lon   = 3;  

i=i+1;
POBS(i).Name = 'NorwSea Gimsoy';
POBS(i).Lat   = 69.00;  
POBS(i).Lon   = 12;  

i=i+1;
POBS(i).Name = 'Barents Sea Fugloya';
POBS(i).Lat   = 73.00;  
POBS(i).Lon   = 20;  

i=i+1;
POBS(i).Name = 'W Iceland Sea Langanes';
POBS(i).Lat   = 67.5;  
POBS(i).Lon   = -13.5;  



return