% Combine individual files TMP_*{01,02,...}.mat -> rest_ice_*mat
% saved in remap_piomas2arc.m
% run in several-jobs mode
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

YR = 1993;
MM = 1;
s_njobs = 6;
pthmat  = '/nexsan/people/ddmitry/hycom/ARCc0.08/103/ice_restart/';

% Do not use f_mat = 3 unless there is a bug
% in remap_piomas2arc.m
f_mat = 1; % = 3 - fix a bug in VSNON and VICEN,

FLDS{1}='aicen';
FLDS{2}='vicen';
FLDS{3}='vsnon';
FLDS{4}='eicen';
FLDS{5}='esnon';

ncat = 5;
nlyr = 4; 

nf = length(FLDS);

for jj=1:nf
  fld = FLDS{jj};
  fmatout = sprintf('%srest_ice_%s_%4.4i%2.2i.mat',...
		 pthmat,fld,YR,MM);
  clear AA
  if f_mat == 3 ...
     & ~strcmp('vsnon',fld) ...
     & ~strcmp('vicen',fld)
    continue
  end
  
  for jb=1:s_njobs
    s_jobid = jb;
    fmatin = sprintf('%sTMP_%s_%4.4i%2.2i_%2.2i.mat',...
                     pthmat,fld,YR,MM,s_jobid);
    fprintf('Loading %s\n',fmatin);
    load(fmatin);
% Correct error
    if f_mat== 3 & strcmp('vsnon',fld);
      VSNON = VSNON./ncat;
      fprintf('Modifying VSNON: %s\n',fmatin);
      save(fmatin,'VSNON','iiS','iiE','INDX');
      continue
    end
    if f_mat== 3 & strcmp('vicen',fld);
      VICEN = VICEN./ncat;  % volume per category (not per each layer*categ!!!)
      fprintf('Modifying VICEN: %s\n',fmatin);
      save(fmatin,'VICEN','iiS','iiE','INDX');
      continue
    end
    
    
%    nrc=size(EICEN,1);
%    nrc_old = size(AA,1);
      
%keyboard    
    switch(lower(fld))
     case('aicen')
      if jb>1,
	ichck = sub_check_fld(AA,AICEN,iiS);
	if ichck == 1, error('Non matching old/new fields'); end;
      end
      nrc=size(AICEN,1);
      AA(iiS:iiE,:) = AICEN(1:nrc,:);
     case('vicen')
      if jb>1,
	ichck = sub_check_fld(AA,VICEN,iiS);
	if ichck == 1, error('Non matching old/new fields'); end;
      end
      nrc=size(VICEN,1);
      AA(iiS:iiE,:) = VICEN(1:nrc,:);
     case('vsnon')
      if jb>1,
	ichck = sub_check_fld(AA,VSNON,iiS);
	if ichck == 1, error('Non matching old/new fields'); end;
      end
      nrc=size(VSNON,1);
      AA(iiS:iiE,:) = VSNON(1:nrc,:);
%      keyboard
     case('eicen')
      if jb>1,
	ichck = sub_check_fld(AA,EICEN,iiS);
	if ichck == 1, error('Non matching old/new fields'); end;
      end
      nrc=size(EICEN,1);
      AA(iiS:iiE,:,:) = EICEN(1:nrc,:,:);
     case('esnon')
      if jb>1,
	ichck = sub_check_fld(AA,ESNON,iiS);
	if ichck == 1, error('Non matching old/new fields'); end;
      end
      nrc=size(ESNON,1);
      AA(iiS:iiE,:) = ESNON(1:nrc,:);
 %     keyboard
    end
  end
  
  if f_mat == 3; continue; end;
  
  switch(lower(fld))
   case('aicen')
    clear AICEN
    AICEN = AA;
    fprintf('Saving %s\n\n',fmatout);
    save(fmatout,'AICEN');
   case('vicen')
    clear VICEN
    VICEN = AA;
    fprintf('Saving %s\n\n',fmatout);
    save(fmatout,'VICEN');
   case('vsnon')
    clear VSNON
    VSNON = AA;
    fprintf('Saving %s\n\n',fmatout);
    save(fmatout,'VSNON');
   case('eicen')
    clear EICEN
    EICEN = AA;
    fprintf('Saving %s\n\n',fmatout);
    save(fmatout,'EICEN');
   case('esnon')
    clear ESNON
    ESNON = AA;
    fprintf('Saving %s\n\n',fmatout);
    save(fmatout,'ESNON');
  end
end


  
