# !/bin/csh -vx
# Extract Arctic domain from Global
# analysis experiemnts 
# GLBa fields and rotate into ARCc grid
# annual mean fields
# Script to run remapGLB.x 
# for N years
# 
set echo

#@ yr1=1992
@ yr1=2013
@ yr2=2015

while (${yr1} <= ${yr2})
  if (`echo ${yr1} 2013 | \
       awk '{if ($2-$1 == 0) print 1; else print 0}'`) then
    setenv Ep 91.0
  else 
    setenv Ep 91.1
  endif
  setenv E `echo ${Ep} | awk '{printf("%3.3i",$1*10)}'`


    setenv D /nexsan/GLBa0.08/expt_${Ep}/data/meanstd/
    setenv FGLB ${E}_archMN.${yr1}_01_${yr1}_12
    setenv FARC ${E}_archMN_GLBb2ARCc.${yr1}_01_${yr1}_12

    rm -f PARAM_${yr1}.dat
    sed -e "s|pathglb.*|pathglb  = ${D}|" \
    -e "s|fnmGLB.*|fnmGLB   = ${FGLB}|"\
    -e "s|fnmARC.*|fnmARC   = ${FARC}|" < PARAM_archMN.dat >! PARAM_${yr1}.dat

    rm -f PARAM.dat
    ln -s PARAM_${yr1}.dat PARAM.dat

    ./remapGLB.x
    wait


  @ yr1++
end

exit 0
