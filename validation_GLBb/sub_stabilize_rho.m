function [Tz,Sz] = sub_stabilize_rho(Tz,Sz,ZMz,Zz);
% Some T/S profiles result in unstable
% density profiles, which complicates 
% MLD calculation - need to stabilize the
% density profile
% Stbilization is done by mixing
% unstable layers and producing a homogeneous 
% layer
% Mixing is performed in the T/S space
%
% D. Dukhovskoy, COAPS FSU, June 2017
%
nlv=length(Tz);

% Check density profile - should be stably-stratified:
% if not -mix the layers
Tz0=Tz;
Sz0=Sz;
Rho0 = sw_dens0(Sz,Tz);
Rho = Rho0;
ik=1;
%keyboard
while ik<=nlv-2
  if isnan(Tz(ik+1)); break; end;
%  fprintf('ik=%i\n',ik);
%  if ik==40; keyboard; end;
  dR = Rho(ik+1)-Rho(ik);
%  if abs(dR)<1e-6; dR=0; end;
  if dR<0;
%    fprintf('Unstable Rho ik=%i\n',ik);
% Check dRho above:
    ik1=ik;
 %   if ik==nlv-1, keyboard; end;
    for iup=ik:-1:2
      udR=Rho(iup)-Rho(iup-1);
      if abs(udR)<1e-6 % homog. layer
	ik1=iup-1;
      else
	break;
      end
    end
% Sum up over upper homogeneous layers 
    dz1=0;
    t1=0;
    s1=0;
    for jlr=ik1:ik
      dz=abs(Zz(jlr+1)-Zz(jlr));
      dz1=dz1+dz;
      t1=t1+dz*Tz(jlr);
      s1=s1+dz*Sz(jlr);
    end

    dz2=abs(Zz(ik+2)-Zz(ik+1));
    t2=dz2*Tz(ik+1);
    s2=dz2*Sz(ik+1);
    t0=(t1+t2)/(dz1+dz2);
    s0=(s1+s2)/(dz1+dz2);
    Tz(ik1:ik+1)=t0;
    Sz(ik1:ik+1)=s0;
    Rho = sw_dens0(Sz,Tz);
    ik=ik1-2;
    ik=max([ik,0]); % move 1 pnt above the homogenized layer
                    % to check if neg. density 
		    % jumps have appeared above
		    % after mixing layers
  end
  ik=ik+1;
  
  mchck=0;
  if mchck>0
    clf;
    subplot(1,3,1);
    plot(Rho,ZMz,'.-');
    hold on;
    plot(Rho0,ZMz,'r.-');
    set(gca,'ylim',[Zz(ik+1)-25 0]);
    subplot(1,3,2);
    plot(Tz,ZMz,'.-');
    hold on;
    plot(Tz0,ZMz,'r.-');
    set(gca,'ylim',[Zz(ik+1)-25 0]);
    subplot(1,3,3);
    plot(Sz,ZMz,'.-');
    hold on;
    plot(Sz0,ZMz,'r.-');
    set(gca,'ylim',[Zz(ik+1)-25 0]);
    keyboard
  end
end;    
%keyboard
chck=0;
if chck>0
  txb = 'sub_stabilize_rho.m';
  Rho0=sw_dens0(Sz0,Tz0);
  figure(12); clf;
  %axes('position',[0.08 0.08 0.2 0.82]);
  plot(Rho0,ZMz,'.-');
  hold on
  plot(Rho,ZMz,'r.-');
  title('Rho0');
  bottom_text(txb,'pwd',1);
  figure(13); clf;
  %axes('position',[0.08 0.08 0.2 0.82]);
  plot(Tz0,ZMz,'.-');
  hold on
  plot(Tz,ZMz,'r.-');
  title('Temperature');
  figure(14); clf;
  %axes('position',[0.08 0.08 0.2 0.82]);
  plot(Sz0,ZMz,'.-');
  hold on
  plot(Sz,ZMz,'r.-');
  title('Salinity');
  keyboard
end;

return