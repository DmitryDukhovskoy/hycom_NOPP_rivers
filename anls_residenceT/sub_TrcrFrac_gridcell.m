% Calculate tracer fruction in 
% grid cells for converting 
% Tracer conc to Gr FW volume etc
function sub_TrcrFrac_gridcell(fextr,expt,Acell,HH,ilv,hZ,IntSrf);

regn='ARCc0.08';
nTr = 1;
iyy=0;
for iyr=1993:2016
		iyy=iyy+1;
		for im=1:12
				dnmb = datenum(iyr,im,15);
				dv0  = datevec(dnmb);
				yr   = dv0(1);
				iday = dnmb-datenum(yr,1,1)+1;

				ism = find(Ygr==dv0(1,1));
				fwf0 = cFWF(ism); % km3
%
% FWFlux from step-function increase (constant Greenland anomaly after 1993)
				fwf0_stp = cFWF_step(ism);

				rr = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ,IntSrf);

% Estimate volume of Greenland surplus FW in grid cell
% = ratio of FW tracer*total FWF (integrated over time 1990-date requested)
				Vfw = fwf0*rr*1e9; % m3 Gr FWF anomaly integrated in time
				VfwS = fwf0_stp*rr*1e9; % Gr FWF anomaly from step-function increase
% Integrate vol of FW in the region
				vfw=nansum(Vfw(INGr))*1e-9; % m3->km3
				VFW(im,iyy)=vfw;

				vfws=nansum(VfwS(INGr))*1e-9;
				VFWstp(im,iyy)=vfws;

		end
end
fprintf('Saving %s\n',fextr);
save(fextr,'VFW','VFWstp');


return
