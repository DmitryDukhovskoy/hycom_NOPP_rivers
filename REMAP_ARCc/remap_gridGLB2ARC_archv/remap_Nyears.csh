# !/bin/csh -vx
# Extract Arctic domain from Global
# GLBb fields and rotate into ARCc grid
# Annual Mean fields
# Script to run remapGLB.x 
# for N years
# 
set echo

#@ yr1=1992
@ yr1=1993
@ yr2=2015

while (${yr1} <= ${yr2})
  if (${yr1} <= 1994) then
    setenv E 190
  else
    setenv E 191
  endif

  setenv D /nexsan/GLBb0.08/GLBb0.08_${E}/data/meanstd/
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
