#! /bin/csh -vx 
#
# Unzip and untar archv ARCc output
# or move untarred archv files from input to data dir
# range years: YR1 - YR2, or 1 year : YR1
# Usage: csh untar_archv.csh YR1 [YR2] [e.g., 1993 1995]
#

set echo
# E 022 - experiments with 0.04 HYCOM-CICE5 forced with CFSv2 different CICE options
# E 023 - expt with 0.04 HYCOM-CICE5 EAP dltEddington radiative param. mushy JRA55 forc
setenv E 123
setenv R ARCc0.08
setenv D0 ${cwd}
setenv D /nexsan/people/ddmitry/hycom/${R}_${E}
#setenv RX "_dltEddTmltEAP"  # name of expt, could be empty
setenv RX "_dltEddTmltEAPJRA"

if ($#argv < 1) then
  echo "At least 1 year has to be specified YR1 (e.g. 1993)"
  exit 1
else if ($#argv == 1) then
  @ YR1 = `echo $argv[1] | awk '{printf("%i", $1)}'`
  @ YR2 = `echo $argv[1] | awk '{printf("%i", $1)}'`
else if ($#argv == 2) then
  @ YR1 = `echo $argv[1] | awk '{printf("%i", $1)}'`
  @ YR2 = `echo $argv[2] | awk '{printf("%i", $1)}'` 
endif

echo "Untarring ${YR1} - ${YR2}"

while ( ${YR1} <= ${YR2} )
  cd ${D}/incoming
  setenv DT ${D}/data/${YR1}${RX}
  mkdir -pv ${DT}

  setenv Y `echo ${YR1} | awk '{printf("%03d", $1-1900)}'`
  
  ls -l ${E}_archv_${Y}*gz >& /dev/null
  if (! $status) then
    mv ${E}_archv_${Y}*gz ${DT}/.
    wait
  else
    echo "No tar files found for year $YR1 in ${D}/incoming ..."
# Check untarred *a, *b files:
#  /bin/rm ${E}_archv.${YR1}_*.txt
#  ls -l ${E}_archv.${YR1}*.[ab] >& /dev/null
#  if (! $status) then
#    foreach flm (${E}_archv.${YR1}*.[ab])
#      /bin/mv $flm ${DT}/.
#    end
#  endif
  endif

  cd ${DT}
  ls -l ${E}_archv_${Y}*gz >& /dev/null
  if (! $status) then
    echo "Tar files in ${DT}, start untarring ..."
    touch list_${Y}.txt
    ls -1 *${Y}*.gz | cat >> list_${Y}.txt
    foreach fltar (${E}_archv_${Y}*gz)
      echo "Untarring ${fltar} ..."
      tar xzvf ${fltar}
      wait

      /bin/mv ${fltar} ${D}/dump/.
    end

# Clean up:
#    /bin/mv ${E}_archv_${Y}*gz ${D}/dump/.
    /bin/rm ${E}_archv.${Y}_*.txt
   endif

  @ YR1++
end  # year loop

exit 0
