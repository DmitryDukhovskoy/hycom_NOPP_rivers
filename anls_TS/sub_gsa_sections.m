function SCT = sub_gsa_sections;
% Coordinates of GSA sections/stations
%% based on Belkin et al., 1998
% to compare with model S changes
ist=1;
SCT(ist).Name= 'FyllaBank';
SCT(ist).XY=[-52.37, 63.95;
             -59.17, 63.2];

ist=ist+1;
SCT(ist).Name = 'SealHamilton';
SCT(ist).XY = [-55.65, 53.23;
	       -52.5, 55.07];

ist=ist+1;
SCT(ist).Name = 'FlemishCap';
SCT(ist).XY = [-52.03, 47;
	       -42.0, 47.0];

ist=ist+1;
SCT(ist).Name = 'Bonavista';
SCT(ist).XY = [-52.97, 48.73;
	       -49., 50];

ist=ist+1;
SCT(ist).Name = 'OWS_B';
SCT(ist).XY = [-51.88, 55.78;
	       -51.4, 55.81];

ist=ist+1;
SCT(ist).Name = 'OWS_C';
SCT(ist).XY = [-35.5, 52.75;
	       -35.2, 52.75];

ist=ist+1;
SCT(ist).Name = 'S10';
SCT(ist).XY = [-48.12, 49.07;
	       -36.68, 52.4];

ist=ist+1;
SCT(ist).Name = 'Rockall';
SCT(ist).XY = [-13.63, 56.67;
	       -6.13,  57.58];

ist=ist+1;
SCT(ist).Name = 'NF_FaroeShetl';
SCT(ist).XY = [-6.2, 62.0;
	       -1.0, 60.93];

ist=ist+1;
SCT(ist).Name = 'OWS_Mike';
SCT(ist).XY = [-2.0, 66.0;
	       -1.8, 66.1];

ist=ist+1;
SCT(ist).Name = 'Svinoy';
SCT(ist).XY = [5.2, 62.37;
	       0.0, 64.67];

ist=ist+1;
SCT(ist).Name = 'Gimsoy';
SCT(ist).XY = [14.08, 68.4;
	       8.2, 70.4];

ist=ist+1;
SCT(ist).Name = 'BearIsl'; % Bjornoya
SCT(ist).XY = [20.0, 70.5;
	       19.17, 74.25];



return

