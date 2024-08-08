#! /bin/bash 
# link ASTER DEM file to a different name start with yyyyMMDDHHMMSS


for line in `ls *dem.tif *meta.txt`
do
  dir1=$(dirname $line)
  ifile=$(basename $line)
  year1=${ifile:16:4}
  month=${ifile:12:4}
  rest0=${ifile:20} #20 to end
  rest=${rest0/'DEM.tif'/'AST_dem.tif'}
  
  oneyr=365
  #remainder=`expr $year1 % 4`   # remainder=$(($year1%4))
  #day=`echo "($doy1 + $doy2 )/2." | bc`
  ofile='WV02_'${line:0:8}"_"${line:8}
  echo $line $ofile
  ln -fs $line $ofile
#exit
done
