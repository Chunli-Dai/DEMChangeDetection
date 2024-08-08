#/bin/sh -f
# usage: ./merge.sh 1 
# 	./merge.sh 2

#option 1: use given list; 2 read list from file tilelist
echo usage: ./merge.sh 2 

option=$1 
echo option is $option

#module load pgc


#Banks
#tilelist=('35_22_2_2' '36_21_2_1' '37_24_2_1');

if [ $option == 1 ] #use given list
then
tilelist=('41_16_1_1');
elif [ $option == 2 ]
then
# get tilelist from a filelist
# ls -d ??_??_?_? > tilelistb
i=0
for line in `cat tilelistb`
do
i=$i+1;
tilelist[$i-1]=$line;
done
fi

nlist=${#tilelist[@]}
echo Total number of tiles: $nlist.

echo ${tilelist[*]}

resr=10 #
echo Output resolution $resr m.

for (( i=1; i<=$nlist; i++ ))
do
tilename=${tilelist[$i-1]} #39_17_2_2 "$yidc"_"$xidc"_"$xidsc"_"$yidsc"

count=`ls $tilename/*/*/*_jump.tif | wc -l`
#ratio=` printf %.1f "$((count *100/625 ))" ` #`echo "$count*100/625" | bc -l`
ratio=`awk "BEGIN {printf \"%.1f\n\", $count*100/625}"`
echo $tilename : $ratio % subtiles have results. 

continue

if [ $count == 0 ]
then
echo No files to merge!
continue
fi

#gdal_merge.py -o 40_18_2_2_jump.tif 40_18_2_2/*/*_jump.tif
#echo gdal_merge.py -o $tilename"_jump.tif" $tilename/*/*_jump.tif
gdal_merge.py -o $tilename"_jump.tif" -ps $resr $resr $tilename/*/*/*_jump.tif
# 49_60_2_2_01_01_bitmask.tif      49_60_2_2_01_01_jumpstd.tif
#49_60_2_2_01_01_eventtimeT1.tif  49_60_2_2_01_01_jump.tif
#49_60_2_2_01_01_eventtimeT2.tif
gdal_merge.py -o $tilename"_jumpstd.tif" -ps $resr $resr $tilename/*/*/*_jumpstd.tif
gdal_merge.py -o $tilename"_eventtimeT1.tif" -ps $resr $resr $tilename/*/*/*_eventtimeT1.tif
gdal_merge.py -o $tilename"_eventtimeT2.tif" -ps $resr $resr $tilename/*/*/*_eventtimeT2.tif
gdal_merge.py -o $tilename"_bitmask.tif" -ps $resr $resr $tilename/*/*/*_bitmask.tif
gdal_merge.py -o $tilename"_nov.tif" -ps $resr $resr $tilename/*/*/*_nov.tif
#gdal_merge.py -o $tilename"_nov.tif" $tilename/*/*_nov.tif

done
