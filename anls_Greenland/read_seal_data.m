% Data of seal T data around Greenland
% Shared with me by Theresa Morison
%
% From: David Sutherland <dsuth@uoregon.edu>
%Date: Fri, May 4, 2018 at 4:17 PM
%Subject: Re: Elephant Seal Data in the Irminger Sea
%To: Theresa Morrison <t4morris@ucsd.edu>
%

%Hi Theresa,

%Sure- I’ve been asked before, so I have it in matlab format already. 
% Note these are from hooded seals (primarily, see paper for details).
%
% The variables in the seals structure should be self-explanatory. There are 4339 dives in there. 
% The only matrix is then #dives by z-levels (4339 x 141), which is the temperature (Tint), 
% where I’ve interpolated the temperatures on to that z-grid. Max_dbar tells you 
% the max depth they dove to. Some of the other variables have to do with comparing 
% bathymetry products to deepest dives (see paper if interested in that). If you want original data, etc. we can chat.
%This file is focused on the SE Greenland region, as described in that paper. 
% To use the other data, we’d have to contact the folks who deployed the tags.
% Good luck on your research- sounds interesting! And, note, that I’ll be down 
% at Scripps this Fall for a sabbatical…
% Dave Sutherland

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

regn = 'ARCc0.08';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_theresa/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

Lmsk = HH*0;
Lmsk(HH<0)=1;
cmsk=[0 0 0; 1 1 1];

fin = sprintf('%sTheresa_sealdata.mat',pthmat);
load(fin); % depth - Z
A = seals;
Z=-Z;

slt = A.LAT;
sln = A.LON;
TM  = A.END_DATE;
T   = A.Tint;

fprintf('Converting lon,lat ---> I,J hycom indices\n');
IJ = sub_XY2indx([sln,slt],LON,LAT);



figure(1); clf;
%contour(HH,[0 0],'k');
pcolor(Lmsk); shading flat;
colormap(cmsk);
caxis([0 1]);
hold on;
contour(HH,[-800 -800],'Color',[0. 0.4 0.8],'Linewidth',1.8);
contour(HH,[-6000:1000:-1000],'Color',[0.7 0.7 0.7]);

plot(IJ(:,1),IJ(:,2),'r.');

axis('equal');
set(gca,'xlim',[450 950],...
        'ylim',[300 800]);


dv1=datevec(min(TM));
dv2=datevec(max(TM));

stl=sprintf('Seal T prof, %i/%i/%i - %i/%i/%i',dv1(1:3),dv2(1:3));
title(stl);

btx = 'read_seal_data.m';
bottom_text(btx,'pwd',1);


figure(2); clf;
axes('Position',[0.1 0.08 0.5 0.82]);
ii=1;
plot(T(ii,:),Z);




