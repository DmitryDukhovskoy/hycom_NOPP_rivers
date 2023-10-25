function A = sub_Glb2Arc(Ain,IJ);
% Get Arctic cap region from the GLBb GOFS3.5 grid
% IJ are corner indices
% Note that the Global grid is not exactly the same as ARCc0.04
%
% Indices where Global grid 
%iir1=1855; 
%jr1=2501;
%ir2=2091;
%jr2=1695;

%keyboard

i11=IJ(1,1);
j11=IJ(1,2);
i12=IJ(2,1);
j12=IJ(2,2);
i21=IJ(3,1);
j21=IJ(3,2);
i22=IJ(4,1);
j22=IJ(4,2);

dm1=Ain(j11:end,i11:i12);
dm2=Ain(j21:end,i22:i21);
dm2=fliplr(dm2);
dm2=flipud(dm2);
A=[dm1;dm2];


return
