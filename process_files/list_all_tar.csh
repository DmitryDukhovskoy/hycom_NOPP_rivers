#! /bin/csh -vx 
#
# Create ASCII file with tar.gz archm 
# for years YR1 YR2
#

set echo

setenv E 110
setenv R ARCc0.08
setenv D0 ${cwd}
setenv D /nexsan/archive/${R}_${E}/incoming

cd ${D}

setenv Y1 093
setenv Y2 105
setenv Fout list_missing_files.txt
touch ${Fout}
rm ${Fout}
touch ${Fout}

while (`echo ${Y1} ${Y2} | awk '{if ($1<=$2) print 1; else print 0}'`)
  echo "Year ${Y1}"
  foreach m(a b c d e f g h i j k l)
    setenv fnm ${E}_archm_${Y1}${m}.tar.gz
    echo ${fnm} | cat >> ${Fout}
  end

  setenv Y1 `echo ${Y1} | awk '{printf("%03d", $1+1)}'`
end

exit 0
