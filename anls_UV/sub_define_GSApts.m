function PTS = sub_define_GSApts(npath); 
% Define set of points for
% different pathways of GSA/ Freshwater
% It typically follows the maximum |U|
% on the annual field

switch(npath);
 case(1)
  PTS.IJ=[         530         504
         518         532
         496         547
         470         550
         430         529
         417         498
         414         450
         420         391
         430         353
         450         332
         477         316
         489         272
         493         244
         506         226
         510         198
         504         167
         521         149
         553         151
         569         164
         566         182
         544         188
         533         208
         556         224
         579         230
         601         244
         632         244
         670         244
         696         254
         712         257
         737         278
         759         305
         775         323
         796         332
         822         334
         838         339
         840         362
         855         369
         873         354
         879         333
         893         346
         894         367
         918         381
         948         399
         972         427
         998         421
        1015         401
        1049         398
        1089         430
        1101         454
        1136         477
        1160         499
        1155         536
        1144         574
        1172         619
        1213         672
        1219         717
        1177         760
        1158         783
        1162         820
        1137         841
        1124         864
        1096         882
        1075         908 ];
  
end


return