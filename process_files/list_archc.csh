#! /bin/csh -vx 
#
# Create ASCII file with missing or incomplete files
# from untarred archc files CICE
# usage list_archc.csh YR1 [YR2]

#set echo

setenv E 110
setenv R ARCc0.08
setenv D0 /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/process_files
#setenv D /nexsan/archive/${R}_${E}/incoming
setenv D /nexsan/archive/${R}_${E}/dump
setenv DT /nexsan/archive/${R}_${E}/data

cd ${D}

#setenv Y1 108
#setenv Y2 108
if ($#argv < 1) then
  echo "At least 1 year has to be specified Y1 (e.g. 1993)"
  exit 1
else if ($#argv == 1) then
  setenv Y1 `echo $argv[1] | awk '{printf("%i", $1)}'`
  setenv Y2 `echo $argv[1] | awk '{printf("%i", $1)}'`
else if ($#argv == 2) then
  setenv Y1 `echo $argv[1] | awk '{printf("%i", $1)}'`
  setenv Y2 `echo $argv[2] | awk '{printf("%i", $1)}'` 
endif

if (`echo ${Y1} | awk '{if ($1>1900) print 1; else print 0}'`) then
  echo "<--- $Y1"
  setenv Y1 `echo $Y1 | awk '{printf("%03i\n", $1-1900)}'`
  setenv Y2 `echo $Y2 | awk '{printf("%03i\n", $1-1900)}'`
  echo "---> $Y1"
endif

echo "Checking missing archc files for ${Y1} - ${Y2}"


setenv Fout ${D}/../list_missing_archc${Y1}_${Y2}.txt
#if (-e ${Fout}) then
#  mv ${Fout} ${Fout}-bkp
#endif

touch ${Fout}
rm ${Fout}
touch ${Fout}

setenv LFile xxx # last saved

while (`echo ${Y1} ${Y2} | awk '{if ($1<=$2) print 1; else print 0}'`)
  setenv YR1 `echo ${Y1} | awk '{printf("%i\n", $1+1900)}'`
  echo "Year ${Y1}  ${YR1}"
  setenv DTY ${DT}/${YR1}
  if (! -d $DTY) then
    foreach AB(a b c d e f g h i j k l)
      setenv fnm ${E}_archc_${Y1}${AB}.tar.gz
      echo ${fnm} | cat >> ${Fout}
      setenv LFile ${fnm}
    end
  else
    cd ${DTY}  
    rm *archm*.txt

    @ dd = 1
    while (`echo ${dd} | awk '{if ($1<=365) print 1; else print 0}'`)

      setenv Day `echo ${dd} | awk '{printf("%3.3d", $1)}'`
      setenv MM `echo "DATES" | awk -f ${D0}/dates.awk y01=${Y1} d01=${Day}\
            | awk '{printf("%2.2i", $2)}'`
      setenv AB `echo "MM2AB" | awk -f ${D0}/dates.awk MM=${MM} | awk '{print $1}'`
      
      echo "Year ${YR1} Day ${Day} Month ${MM}, ${AB}" 
      setenv fina ${E}_archc.${YR1}_${Day}_12.a
      setenv finb ${E}_archc.${YR1}_${Day}_12.b
      setenv fnm ${E}_archc_${Y1}${AB}.tar.gz

# Check if files exist and right size
      if (-f ${fina} && -f ${finb}) then
        setenv sz `ls -l ${fina} | awk '{printf("%i", $5*1e-9)}'`
        if (`echo ${sz} | awk '{if ($1<7) print 1; else print 0}'`) then
          echo "Incomplete transfer ${fina} ${sz} Gb"
          if ( ${fnm} == ${LFile} ) then
            @ dd++
            continue  # already saved
          else
            setenv LFile ${fnm}
            echo ${fnm} | cat >> ${Fout}
          endif
        endif
      else
        echo "Files is missing: ${fina} or ${finb}"
	if ( ${fnm} == ${LFile} ) then
          @ dd++
	  continue  # already saved
	else
          setenv LFile ${fnm}
	  echo ${fnm} | cat >> ${Fout}
	endif
      endif

      @ dd++
    end
  endif

  setenv Y1 `echo ${Y1} | awk '{printf("%03d", $1+1)}'`
end

exit 0
