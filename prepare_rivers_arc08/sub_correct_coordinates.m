function RVR = sub_correct_coordinates(RVR);
% Adjust river mouth locations 
% to get correct positioning in ARCc
dmm = RVR.Lon;
a2 = length(dmm);

for ik = 1:a2
  rnm = RVR.Riv_name(ik,:);
  rnm = deblank(rnm);
  switch(rnm),
   case('Pechora')
    RVR.Lon(ik) = 54.39;
    RVR.Lat(ik) = 68.15;
   case('Mezen')
    RVR.Lon(ik) = 44.14;
    RVR.Lat(ik) = 65.87;
   case('Severnaya Dvina')
    RVR.Lon(ik) = 40.37;
    RVR.Lat(ik) = 64.58;
   case('Onega')
    RVR.Lon(ik) = 38.036;
    RVR.Lat(ik) = 63.918;
   case('Varzuga')
    RVR.Lon(ik) = 36.83;
    RVR.Lat(ik) = 66.295;
   case('Ob')
    RVR.Lon(ik) = 72.25;
    RVR.Lat(ik) = 68.895;
   case('Nadym')
    RVR.Lon(ik) = 74.61;
    RVR.Lat(ik) = 68.91;
   case('Pur')
%    RVR.Lon(ik) = 77.665;
%    RVR.Lat(ik) = 67.417;
    RVR.Lon(ik) = 75.56;
    RVR.Lat(ik) = 71.14;
   case('Taz')
%    RVR.Lon(ik) = 78.6; % actual locations, but move a bit to match ARCc
%    RVR.Lat(ik) = 67.48; % coastal line
    RVR.Lon(ik) = 77.97;
    RVR.Lat(ik) = 71.09;
   case('Yenisey')
    RVR.Lon(ik) = 82.723;
    RVR.Lat(ik) = 71;
   case('Pyasina')
    RVR.Lon(ik) = 85.967;
    RVR.Lat(ik) = 73.862;
   case('Khatanga')
    RVR.Lon(ik) = 109.66;
    RVR.Lat(ik) = 74.195;
   case('Anabar')
    RVR.Lon(ik) = 113.56;
    RVR.Lat(ik) = 73.35;
   case('Olenek')
    RVR.Lon(ik) = 119.74;
    RVR.Lat(ik) = 72.99;
   case('Lena A');
    RVR.Lon(ik) = 123.41;
    RVR.Lat(ik) = 72.94;
   case('Lena B');
    RVR.Lon(ik) = 128.94;
    RVR.Lat(ik) = 72.65;
   case('Lena C');
    RVR.Lon(ik) = 128.695;
    RVR.Lat(ik) = 71.98;
   case('Kolyma');
    RVR.Lon(ik) = 161.44;
    RVR.Lat(ik) = 69.45;
   case('Yana');
    RVR.Lon(ik) = 136.71;
    RVR.Lat(ik) = 71.51;
   case('Indigirka');
    RVR.Lon(ik) = 150.397;
    RVR.Lat(ik) = 71.295;
   case('Palyavaam');
    RVR.Lon(ik) = 170.399;
    RVR.Lat(ik) = 68.814;
   case('Mackenzie A');
    RVR.Lon(ik) = -135.58;
    RVR.Lat(ik) = 68.97;
   case('Mackenzie B');
    RVR.Lon(ik) = -135.27;
    RVR.Lat(ik) = 69.37;
   case('Mackenzie C');
    RVR.Lon(ik) = -134.19;
    RVR.Lat(ik) = 69.25;
   case('Anderson');
    RVR.Lon(ik) = -128.49;
    RVR.Lat(ik) = 69.893;
  end
  
end;


return