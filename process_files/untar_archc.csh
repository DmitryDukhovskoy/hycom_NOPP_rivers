#! /bin/csh -vx 
#
# Unzip and untar archm ARCc output
# or move untarred archm files from input to data dir
# range years: YR1 - YR2, or 1 year : YR1
# Usage: csh untar_archm.csh YR1 [YR2] [e.g., 1993 1995]
#

set echo

setenv E 123
setenv R ARCc0.08
setenv D0 ${cwd}
setenv RX "_dltEddTmltEAPJRA"  # name of expt, could be empty

if ($E == 110) then
  setenv D /nexsan/archive/${R}_${E}
else if ($E == 112) then
  setenv D /nexsan/hycom/${R}_${E}
else
  setenv D /nexsan/people/ddmitry/hycom/${R}_${E}
endif

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

echo "Untarring CICE ${YR1} - ${YR2}"

while ( ${YR1} <= ${YR2} )
  cd ${D}/incoming
  setenv DT ${D}/data/${YR1}_cice${RX}
  mkdir -pv ${DT}

  setenv Y `echo ${YR1} | awk '{printf("%03d", $1-1900)}'`
  
  ls -l ${E}_archc_${Y}*gz >& /dev/null
  if (! $status) then
    mv ${E}_archc_${Y}*gz ${DT}/.
    wait
  else
    echo "No tar files found for year $YR1 in ${D}/incoming ..."
  endif

  cd ${DT}
  ls -l ${E}_archc_${Y}*gz >& /dev/null
  if (! $status) then
    echo "Tar files in ${DT}, start untarring ..."
    touch list_${Y}.txt
    ls -1 *${Y}*.gz | cat >> list_${Y}.txt
    foreach fltar (${E}_archc_${Y}*gz)
      echo "Untarring ${fltar} ..."
      tar xzvf ${fltar}
      wait

      /bin/mv ${fltar} ${D}/dump/.
    end

    chmod 755 *.nc

# Clean up:
#    /bin/mv ${E}_archc_${Y}*gz ${D}/dump/.
    /bin/rm ${E}_archc.${Y}_*.txt
   endif

  @ YR1++
end  # year loop

exit 0
