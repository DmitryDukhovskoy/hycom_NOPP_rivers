# !/bin/csh -vx
# Extract Arctic domain from Global
# GLBb fields and rotate into ARCc grid
# Monthly Mean fields
# Script to run remapGLB.x 
# for N years
# 
set echo

@ yr1=1993
#@ yr1=2008
@ yr2=2011

while (${yr1} <= ${yr2})
  if (${yr1} <= 1995) then
    setenv E 190
  else
    setenv E 191
  endif

  foreach mn ( 02 03 04 05 06 08 09 10 11 12 )
    if (`echo ${mn} 8 ${yr1} | \
          awk '{if ($2-$1 <= 0 && $3 == 1995) print 1; \
                else print 0}'`) then
        setenv E 191
    endif

    setenv D /nexsan/GLBb0.08/GLBb0.08_${E}/data/meanstd/
    setenv FGLB ${E}_archMN.${yr1}_${mn}
    setenv FARC ${E}_archMN_GLBb2ARCc.${yr1}_${mn}

    rm -f PARAM_${yr1}.dat
    sed -e "s|pathglb.*|pathglb  = ${D}|" \
    -e "s|fnmGLB.*|fnmGLB   = ${FGLB}|"\
    -e "s|fnmARC.*|fnmARC   = ${FARC}|" < PARAM_archMNm.dat >! PARAM_${yr1}.dat

    rm -f PARAM.dat
    ln -s PARAM_${yr1}.dat PARAM.dat

    ./remapGLB.x
    wait

  end  # month cycles

  @ yr1++
end

exit 0
