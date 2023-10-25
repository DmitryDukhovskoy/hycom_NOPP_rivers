functions [IJr,hmin] = sub_regions_vort(regn);
% Regions for calculation area-mean vorticity

switch(regn),
 case('natl');
%  ii1=41;
%  ii2=122;
%  jj1=9;
%  jj2=113;
  IJr = [41   9
   41  113
   122 113
   41   9];
  hmin= -200;
 case('arctA');
  IJr=[ 67   212
    75   213
   130   198
   167   133
   167    87
   142    78
   122   103
    96   108
    81   113
    30   169];
  hmin = -10;
 case ('arctB');
  IJr = [    52   185
    69   189
    90   187
   118   174
   132   169
   131   136
   121   110
   106   108
    87   122
    58   156
    46   170
    48   181];
  hmin=-10;
 case('arctic'); % exclude near-coastal regions
  IJr=[     61   188
    75   187
    90   189
   107   175
   124   175
   130   174
   134   167
   132   154
   130   149
   133   144
   135   139
   143   138
   145   128
   144   121
   148   113
   151   107
   152    94
   147    89
   139    85
   122   104
   103   110
    96   118
    78   130
    65   141
    59   151
    51   163
    47   173
    48   180
    56   186];
  hmin = -10;
end


return 


