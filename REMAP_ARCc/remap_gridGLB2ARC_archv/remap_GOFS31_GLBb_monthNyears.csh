# !/bin/csh -vx
# Extract Arctic domain from Global
# analysis experiemnts GOFS3.1
# GLBa fields and rotate into ARCc grid
# monthyl mean fields
# Script to run remapGLB.x 
# for N years
# 
set echo

#@ yr1=1992
@ yr1=2013
@ yr2=2013

while (${yr1} <= ${yr2})
#  foreach mn ( 01 02 03 04 05 06 07 08 09 10 11 12 )
  foreach mn ( 08 )
    if (`echo ${yr1} | \
	 awk '{if ($1 < 2013) print 1; else print 0}'`) then
      setenv Ep 90.9
    else if (`echo ${yr1} | \
	 awk '{if ($1 == 2013) print 1; else print 0}'`) then
      setenv Ep 91.0
    else if (`echo ${yr1} | \
	 awk '{if ($1 < 2016 && $1 >= 2014) print 1; else print 0}'`) then
      setenv Ep 91.1
    else if (`echo ${yr1} | \
	 awk '{if ($1 >= 2016) print 1; else print 0}'`) then
      setenv Ep 91.2
    endif

    if (`echo ${yr1} ${mn}| \
          awk '{if ($1 == 2013 && $2 <= 7) print 1; \
                else print 0}'`) then
        setenv Ep 90.9
    else if (`echo ${yr1} ${mn}| \
          awk '{if ($1 == 2014 && $2 <= 3) print 1; \
                else print 0}'`) then
        setenv Ep 91.0
    else if (`echo ${yr1} ${mn}| \
          awk '{if ($1 == 2016 && $2 <= 3) print 1; \
                else print 0}'`) then
        setenv Ep 91.1
    endif

    setenv E `echo ${Ep} | awk '{printf("%3.3i",$1*10)}'`
    setenv D /nexsan/GLBa0.08/expt_${Ep}/data/meanstd/
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

  end  # months

  @ yr1++
end



exit 0
