% Data from
%http://www.whoi.edu/page/preview.do?pid=66578
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear



AOO=[1946	5
1947	3.33333
1948	2.5
1949	2.38095
1950	2.77778
1951	2.85714
1952	2
1953	-3.88889
1954	-2.77778
1955	-1.25
1956	-2.85714
1957	-2
1958	1.33333
1959	1.25
1960	3.33333
1961	2.59259
1962	1.42857
1963	-1.33333
1964	-1.25
1965	-1.33333
1966	-0.833333
1967	-1.66667
1968	-2
1969	-1.66667
1970	-1.85185
1971	-1.25
1972	2.22222
1973	3.33333
1974	1.66667
1975	2
1976	2.22222
1977	3.33333
1978	1.66667
1979	2
1980	-1.66667
1981	-1.11111
1982	-1.66667
1983	-1.66667
1984	-2
1985	2.38095
1986	2.22222
1987	2.85714
1988	2.22222
1989	-5
1990	-1.90476
1991	-2.38095
1992	-1.66667
1993	-0.333333
1994	-3.33333
1995	-1.66667
1996	-2.33333
1997	1.38889
1998	2.5
1999	1.78571
2000	0.833333
2001	1.51515
2002	1.5625
2003	0.925926
2004	3.88889
2005	2.77778
2006	1.85185
2007	4.58333
2008	2.91667
2009	1.92308
2010	1.85185
2011	1.96078
2012	0.833333
2013	0.952381
2014	1.4881
2015	1.19048
2016.   1.89
2017.   0.38 ];

yrs = AOO(:,1);
aoo = AOO(:,2);

ips = find(aoo>0);
ing = find(aoo<0);


figure(1); clf;
axes('Position',[0.08 0.55 0.85 0.38]);
hold on
hb=bar(yrs(ips),aoo(ips),0.98);
set(hb,'edgecolor',[0.3 0.3 0.3],...
       'FaceColor',[0 0.5 0.9]);

hb=bar(yrs(ing),aoo(ing),0.98);
set(hb,'edgecolor',[0.3 0.3 0.3],...
       'FaceColor',[0.9 0.6 0]);

set(gca,'tickdir','out',...
	'xlim',[yrs(1)-0.5 2017.5],...
	'ylim',[-5 5],...
	'xtick',[1945:5:2020],...
	'xminortick','on',...
	'ytick',[-5:5],...
	'Fontsize',16,...
	'ygrid','on',...
	'xgrid','on');
title('AOO Index');

axes('Position',[0.08 0.08 0.85 0.38]);
hold on
hb=bar(yrs(ips),aoo(ips),0.98);
set(hb,'edgecolor',[0.3 0.3 0.3],...
       'FaceColor',[0 0.5 0.9]);

hb=bar(yrs(ing),aoo(ing),0.98);
set(hb,'edgecolor',[0.3 0.3 0.3],...
       'FaceColor',[0.9 0.6 0]);

set(gca,'tickdir','out',...
	'xlim',[1993 2017.5],...
	'ylim',[-5 5],...
	'xtick',[1945:2020],...
	'xminortick','on',...
	'ytick',[-5:5],...
	'ygrid','on',...
	'xgrid','on',...
	'Fontsize',16,...
	'ydir','reverse');

set(gcf,'Position',[944 508 1596 831]);
btx = 'aoo_index.m';
bottom_text(btx,'pwd',1);


