#! /bin/bash 
# link ASTER DEM file to a different name start with yyyyMMDDHHMMSS

if [ $# -lt 1 ]
then
  echo "usage: $0 filelist "
  exit
fi

filelist=$1

for line in ` cat $filelist `
do
  #GSM-2_2012336-2012366_0029_EIGEN_G---_0005
  #2003_016_0.txt
#AST14DMO_00303162016150419_20200616153100_18885_DEM.tif
#AST0_20160316150419_20200616153100_18885_dem.tif
# echo $line
  dir1=$(dirname $line)
  cd $dir1
  ifile=$(basename $line)
  year1=${ifile:16:4}
  month=${ifile:12:4}
  rest0=${ifile:20} #20 to end
  rest=${rest0/'DEM.tif'/'AST_dem.tif'}
  
  oneyr=365
  #remainder=`expr $year1 % 4`   # remainder=$(($year1%4))
  #day=`echo "($doy1 + $doy2 )/2." | bc`
  ofile="$year1""$month""$rest"
  echo $ifile $ofile  
  ln -fs $ifile $ofile
  cd -

done
