function RV = sub_define_rvr_clim
% Define locations for searcing Arctic rivers
% on ARCc grid with old rivers!
% Also prepare river climatology
% From AOMIP
% Specify rivers
% in terms of HYCOM ARCc grid;
% Conversion runoff: 0.031536* m3/s -> km3/yr
% Annual gauged vol. flux is 2456 km3/yr (AWI climatology)
% AOMIP river runoff data 
RV = struct;
ir = 0;
RNMS{1}='Mackenzie';
RNMS{2}='Yenisey';
RNMS{3}='Lena';
RNMS{4}='Kolyma';
RNMS{5}='Yana';
RNMS{6}='Severnaya Dvina';
RNMS{7}='Pechora';
RNMS{8}='Ob';
RNMS{9}='Olenek';




ir=ir+1;
RV(ir).Name='Mackenzie';
RV(ir).I=[400,450]; % indices of old rivers, climatology
RV(ir).J=[1600,1700];
RV(ir).Clim_m3s=[3814.35, 3605.9, 3349.0, 3384.45,...
		 13218.8, 21413.0, 17854.7, 13984.2, ...
		 11268.7, 9038.15, 4756.0, 3595.0];

ir=ir+1;
RV(ir).Name='Severnaya Dvina';
RV(ir).I=[147,1530];
RV(ir).J=[650,800];
RV(ir).Clim_m3s=[1033.45, 826.6, 724.9, 2415.1, ...
		 13839.4, 7029.7, 2943.8, 2149.4,...
		 2320.4, 2912.4, 2363.0, 1400.7];

ir=ir+1;
RV(ir).Name='Pechora';
RV(ir).I=[1490,1530];
RV(ir).J=[900,1000];
RV(ir).Clim_m3s=[959.4, 773.8, 695.4, 950.2,...
		 15502.6, 17126.4, 5534.2, 3227.8,...
		 3917.0, 4197.6, 1894.4, 1277.4];

ir=ir+1;
RV(ir).Name='Ob';
RV(ir).I=[1520,1590];
RV(ir).J=[1110,1200];
RV(ir).Clim_m3s=[4986.7, 4120.2, 3635.8, 3698.0,...
		 15122.7, 36715.9, 31694.8, 23122.7,...
		 14747.4, 11000.7, 6695.4, 5734.2];

ir=ir+1;
RV(ir).Name='Yenisey';
RV(ir).I=[1450,1600];
RV(ir).J=[1200,1260];
RV(ir).Clim_m3s=[6038.6, 6022.9, 5983.9, 6001.3, ...
		 27533.5, 77386.6, 26586.8, 17485.4,...
		 16896.1, 13969.2, 6855.9, 5839.8];

ir=ir+1;
RV(ir).Name='Pyasina';
RV(ir).I=[1400,1500];
RV(ir).J=[1300,1350];
RV(ir).Clim_m3s=[501.1, 501.1, 501.1, 501.1,...
		 501.1, 7516.7, 10022.3, 2818.8,...
		 3758.4, 4071.6, 1252.8, 689.0];

ir=ir+1;
RV(ir).Name='Olenek';
RV(ir).I=[1300,1380];
RV(ir).J=[1550,1600];
RV(ir).Clim_m3s=[7.00000, 3.20000, 2.10000, 1.65000,...
		 306.850, 7160.00, 2203.50, 846.100,...
		 1050.45, 308.900, 78.7500, 25.1000];

ir=ir+1;
RV(ir).Name='Lena';
RV(ir).I=[1150,1300];
RV(ir).J=[1600,1750];
RV(ir).Clim_m3s=[2783.04, 2136.78, 1651.78, 1350.20,...
		 6235.90, 73917.4, 39683.4, 27340.0,...
		 24126.0, 13771.7, 3502.06, 2927.86];

ir=ir+1;
RV(ir).Name='Kolyma';
RV(ir).I=[930,1000];
RV(ir).J=[1840,1870];
RV(ir).Clim_m3s=[260.0, 191.765, 190.706, 162.647,...
		 1739.24, 13996.6, 7527.12, 6011.70,...
		 4804.71, 1646.18, 436.059, 340.647];


return