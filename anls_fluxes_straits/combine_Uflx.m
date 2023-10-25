% Combines strait fluxes for several years
function FLX = combine_Uflx(pthmat,YR1,YR2,expt,nmex,res); 

for YR = YR1:YR2
  fmatout=sprintf('%shycom%3.3i_%3.3i_%s_StraitFluxes_%4.4i.mat',...
                  pthmat,res*100,expt,nmex,YR);
  FF = sub_readUflx(fmatout); 
  if YR == YR1
    FLX = FF;
  else
    for ik=1:length(FF)
      dmm1=FLX(ik).Vol;
      dmm2=FF(ik).Vol;
      FLX(ik).Vol=[dmm1;dmm2];

      dmm1=FLX(ik).TM;
      dmm2=FF(ik).TM;
      FLX(ik).TM=[dmm1;dmm2];
    end
  end
end

return

