#! /bin/csh -vx
# Replace old data storage directories
# /Net/ocean, /Net/tholia, /Net/Movies0
# Run it within 1 directory at a time
#
setenv DR /home/ddmitry/XXX

if ($#argv < 1) then
  echo "Indicate directory"
  echo "Usage csh replace_net_directories.csh /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers" 
  exit 1
else
  setenv DR $argv[1]
endif

setenv DO /nexsan/people/ddmitry/Net_ocean
setenv DT /nexsan/people/ddmitry/Net_tholia
setenv DT2 /nexsan/people/ddmitry/Net_data2

cd ${DR}
pwd

foreach file ( *.m )
  echo $file
  /bin/cp $file $file-bkp
  sed -e "s|/Net/ocean/ddmitry|$DO|g" \
      -e "s|/Net/tholia/ddmitry|$DT|g" \
      -e "s|/Net/data2/ddmitry|$DT2|g" ${file}-bkp >! ${file}

end

/bin/rm *.m-bkp
exit 0
