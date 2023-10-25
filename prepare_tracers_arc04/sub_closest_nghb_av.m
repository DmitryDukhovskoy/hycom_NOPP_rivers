    function Av = sub_closest_nghb_av(A,Nav,Npts);
%
% Smooth field A that has lots of NaNs
% using only closest Npts neighbour points
% weights are equal, 
% if not - need to modify code
% X, Y - grid incides or coordinates
Ip=find(~isnan(A));
if isempty(Ip), 
  fprintf('sub_closest_nghb_av: no nonNaN values to average ...\n');
  fprintf('sub_closest_nghb_av: QUITTING no averaging is done\n');
  return
end

[m,n]=size(A);
[X,Y]=meshgrid([1:n],[1:m]);

%wght=1/Npts;

nI=length(Ip);
Av=A*0;
clear PP;
PP=struct;
for nf=1:Nav
  fprintf('Averaging: cycle %i out of %i\n',nf,Nav);
  for ii=1:nI;
%    disp(ii);
    if mod(ii,2000)==0 & nf==1
      fprintf('Averaging: %i%% done\n',...
	      round(ii/nI*100));
    end
    
    [j,i]=ind2sub(size(A),Ip(ii));
    if nf==1
      x0=X(j,i);
      y0=Y(j,i);
    %  D=distance_spheric_coord(y0,x0,Y,X);
      D=sqrt((X-x0).^2+(Y-y0).^2);
      Ds=D(Ip);
      [sDs,Is]=sort(Ds);
      Is=Is(1:Npts);
      Iav=Ip(Is);

      if length(Is)<Npts
	error('Npts > available points for averaging');
      end
      PP(ii).Iav=Iav;
    else
      Iav=PP(ii).Iav;
    end
    
    dmm=mean(A(Iav));
    if isnan(dmm)
      error('Check averaging - nans');
    end
    
    Av(j,i)=dmm;

  end
  A=Av;
end;
  
return



