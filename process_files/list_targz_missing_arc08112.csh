#! /bin/csh -vx 
#
# Create ASCII file with missing or incomplete files
#

set echo

setenv E 112
setenv R ARCc0.08
setenv D0 ${cwd}
setenv D /nexsan/hycom/${R}_${E}/incoming

cd ${D}

setenv Y1 093
setenv Y2 105
setenv Fout ${D}/../list_arc008_${E}_targz.txt
touch ${Fout}
rm ${Fout}
touch ${Fout}

while (`echo ${Y1} ${Y2} | awk '{if ($1<=$2) print 1; else print 0}'`)
  echo "Year ${Y1}"
  foreach m(a b c d e f g h i j k l)
    setenv fnm ${E}_archm_${Y1}${m}.tar.gz
    if (-f ${fnm}) then
      setenv sz `ls -l ${fnm} | awk '{printf("%i", $5*1e-9)}'`
      if (`echo ${sz} | awk '{if ($1<50) print 1; else print 0}'`) then
        echo ${fnm} | cat >> ${Fout}
      endif
    else
      echo ${fnm} | cat >> ${Fout}
    endif
  end

  setenv Y1 `echo ${Y1} | awk '{printf("%03d", $1+1)}'`
end

exit 0
