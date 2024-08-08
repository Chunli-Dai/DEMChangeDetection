#/bin/sh -f

if [ $# -lt 1 ]
then
  echo "usage: $0 filelist "
  exit
fi

filelist=$1
ores=$2  #8m

odir='./orthosubset/';
odir='./orthosubset2mfixgrid/';
#for line in `cat filelistpgc`; do
#for line in `cat nadirviewlist.txt`; do
for line in ` cat $filelist `; do
infile=$line
dir1=$(dirname $line)
odir=$dir1

filename=$(basename $line)
#"${firstString/Suzi/$secondString}
ofile1=$odir/${filename/.tif/_sub.tif}
ofile2=$odir/${filename/.tif/_sub_utm2m.tif}

ofile3=$odir/${filename/.tif/_utm.tif}

#gdal_translate -projwin -3113600 730300 -3102200 722000 /home/dai.56/chunliwork/sdm/datasite22/qb4202_Barry_Glacier_Pan_Mono_Ortho_Imagery_2020may08/ortho/WV03_20170203215639_10400100271A1300_17FEB03215639-P1BS-501172500040_01_P008_u16ns3413.tif WV03_20170203215639_10400100271A1300_17FEB03215639-P1BS-501172500040_01_P008_u16ns3413_sub.tif
#gdalwarp WV02_20150516_1030010041649400_103001004243B900_2m_lsf_seg4_ortho_prep.tif WV02_20150516utm_2m.tif -t_srs epsg:32606 -tr 2 2
#-te xmin ymin xmax ymax -te 432936 6772572 443098 6785034

#gdal_translate -projwin -3113600 730300 -3102200 722000 $infile $ofile1
#gdalwarp $ofile1 $ofile2 -t_srs epsg:32606 -tr 2 2 -te 432936 6772572 443098 6785034

#boundary of gdalinfo W2W2_20180611_103001007E7F1900_103001007F843600_2m_lsf_seg1_ortho.tif 
#-3429122 -3406392 2374856 2396770

gdalwarp $infile $ofile3 -t_srs epsg:3413 -tr $ores $ores -te -3429122 2374856 -3406392 2396770
#gdalwarp $infile $ofile3 -t_srs epsg:3413 -tr $ores $ores 

#rm $ofile1
#exit
done
