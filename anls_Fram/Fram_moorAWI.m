function SBE=Fram_moorAWI;
% Note that mooring in the western Fram are NPI moorings
% https://data.npolar.no/dataset/9e01a801-cddf-4f2d-8ed5-b367ad73ea41
%
% Get locations X,Y of 
% Deep water moorings F1 - F16
% of AWI observations (SBE records)
% 1997 - 2010 - some moorings were added after 1997
% See website: http://cmrecords.net/quick/arctic/fram/fram.htm
% Note that the moorings were rotated 
% every year and exact positions were slightly
% different
% Here I use the initial positions on 2002
i=0;
i=i+1;
SBE(i).Name = 'F1';
SBE(i).lat_lon = [78.830, 8.637];

i=i+1;
SBE(i).Name = 'F2';
SBE(i).lat_lon = [78.760, 8.335];

i=i+1;
SBE(i).Name = 'F3';
SBE(i).lat_lon = [78.847, 7.962];

i=i+1;
SBE(i).Name = 'F4';
SBE(i).lat_lon = [78.828, 6.933];

i=i+1;
SBE(i).Name = 'F5';
SBE(i).lat_lon = [78.823, 6.45];

i=i+1;
SBE(i).Name = 'F6';
SBE(i).lat_lon = [78.828, 4.995];

i=i+1;
SBE(i).Name = 'F7';
SBE(i).lat_lon = [78.812, 4.045];

i=i+1;
SBE(i).Name = 'F8';
SBE(i).lat_lon = [78.832, 2.612];

i=i+1;
SBE(i).Name = 'F9';
SBE(i).lat_lon = [78.99,0.262];

i=i+1;
SBE(i).Name = 'F10';
SBE(i).lat_lon = [79.008, -2.042];

i=i+1;
SBE(i).Name = 'F11';
SBE(i).lat_lon = [79.015, -3.018];

i=i+1;
SBE(i).Name = 'F12';
SBE(i).lat_lon = [78.997, -4.228];

i=i+1;
SBE(i).Name = 'F13';
SBE(i).lat_lon = [78.972,-5.319];

i=i+1;
SBE(i).Name = 'F14';
SBE(i).lat_lon = [79.029, -6.582];

i=i+1;
SBE(i).Name = 'F15';
SBE(i).lat_lon = [78.833, 1.612];

i=i+1;
SBE(i).Name = 'F16';
SBE(i).lat_lon = [78.835, 0.4];

% Mooring on shelf
i=i+1;
SBE(i).Name = 'F17';
SBE(i).lat_lon = [79.029, -8];

return
