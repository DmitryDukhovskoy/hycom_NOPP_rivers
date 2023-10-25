function [XT,YT] = read_tiles;
% Read tiles extracted from output log file
%
A = load('tiles.dat');
NN = size(A,1);
XT(:,1) = A(:,2);
XT(:,2) = A(:,4);
YT(:,1) = A(:,3);
YT(:,2) = A(:,5);


return
