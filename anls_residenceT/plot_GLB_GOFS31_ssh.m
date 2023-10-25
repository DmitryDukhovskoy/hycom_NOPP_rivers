% Mean SSH over the SPG region
% from Global HYCOM reanalysis GOFS3.1 (expt_53.x)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 0; % =0 - skip reading hycom files, load saved mean SSH


rg=9806;  % convert pressure to depth, m
Sref=35; % N.Atl. is too saline
TV = '11';

pthARC = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mnth_mean/';
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_ssh/';
%ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';

fmat = sprintf('%sGLB_GOFS31_annualSSH_SPG.mat',pthmat);


i1=1210;
i2=2445;
ni=i2-i1;
j1=2100;
j2=3100;
nj=j2-j1;

if s_mat==1
  iyy=0;
  for iyr=1994:2015
    iyy=iyy+1;
    pthbin=sprintf('/nexsan/hycom/GLBv0.08_53X/data/%i/',iyr);
  %  fmat = sprintf('%sGLB_GOFS31_annualSSH_SPG_%4.4i.mat',pthmat,iyr);

    SSHM=zeros([nj,ni]);
    icc=0;

    for imo=1:12
      for iday=1:7:31
	if iyr<1999
	 expt=530;
	elseif iyr==1999 & imo<4
	  expt=530;
	elseif iyr==1999 & imo>=4
	  expt=531;
	elseif iyr==2000
	  expt=531;
	elseif iyr>=2001 & iyr<2003
	  expt=532;
	elseif iyr==2003 & imo<=6
	  expt=532;
	elseif iyr==2003 & imo>6
	  expt=533;
	elseif iyr>2003 & iyr<2005
	  expt=533;
	elseif iyr==2005 & imo<=6
	  expt=533;
	elseif iyr==2005 & imo>6
	  expt=534;
	elseif iyr>2005 & iyr<2007
	  expt=534;
	elseif iyr==2007 & imo<=6
	  expt=534;
	elseif iyr>2007 & iyr<2009
	  expt=535;
	elseif iyr==2009 & imo<=6
	  expt=535;
	elseif iyr==2009 & imo>6
	  expt=536;
	elseif iyr>2009 & iyr<2011
	  expt=536;
	elseif iyr==2011 & imo<=6
	  expt=536;
	elseif iyr==2011 & imo>6
	  expt=537;
	elseif iyr>2011 & iyr<2013
	  expt=537;
	elseif iyr==2013 & imo<=6
	  expt=537;
	elseif iyr==2013 & imo>6
	  expt=538;
	elseif iyr==2014
	  expt=538;
	elseif iyr==2015
	  expt=539;
	end

	ihr=0;
	fnm=sprintf('%shycom_GLBv0.08_%3.3i_%4.4i%2.2i%2.2i12_t%3.3i.nc',...
		    pthbin,expt,iyr,imo,iday,ihr);

	if ~exist(fnm,'file');
	  while ihr<21
	    ihr=ihr+3;
	    fnm=sprintf('%shycom_GLBv0.08_%3.3i_%4.4i%2.2i%2.2i12_t%3.3i.nc',...
		    pthbin,expt,iyr,imo,iday,ihr);
	     if exist(fnm,'file'); break; end;
	  end
	end

	if ~exist(fnm,'file');
	  fprintf('Does not exist: %s\n',fnm);
	  continue;
	end

	fprintf('Reading %s\n',fnm);

	ssh=squeeze(nc_varget(fnm,'surf_el',[0,j1,i1],[1,nj,ni]));
  %      sshm=nanmean(nanmean(ssh));
  %      ssh=ssh-sshm;
	icc=icc+1;
	SSHM=SSHM+ssh;


      end   % iday
    end     % imo

    SSH(iyy).ssh_mean=SSHM./icc;
    SSH(iyy).year=iyr;
    SSH(iyy).reanalysis=fnm;

    fprintf('Saving %s\n\n',fmat);
    save(fmat,'SSH');

  end       % iyr
else
  fprintf('Loading %s\n',fmat);
  load(fmat);
end


nyr=length(SSH);
A=SSH(1).ssh_mean*0;
for ii=1:nyr
  A=A+SSH(ii).ssh_mean;
end
A=A./nyr;
sshm=nanmean(nanmean(A));
ssh=A-sshm;

clr1=colormap_blue(100);
clr2=colormap_red(100);
for ik=1:10
  clr2(ik,:)=[1 1 1];
  clr1(ik,:)=[1 1 1];
end
clr1=flipud(clr1);
cmp=[clr1;clr2];
cmp=smooth_colormap(cmp,5);


figure(1); clf;
pcolor(ssh); shading flat;
caxis([-0.8 0.8]);
colormap(cmp);

axis('equal');
set(gca,'xtick',[],...
	'ytick',[],...
	'xlim',[100 1100],...
	'ylim',[1 900],...
	'color',[0 0 0]);

pbb=colorbar;
set(pbb,'Fontsize',14);
title('Demeaned ssh, HYCOM+NCODA Reanalysis GOFS3.1, 1994-2015');

btx='plot_GLB_GOFS31_ssh.m';
bottom_text(btx,'pwd',1);





