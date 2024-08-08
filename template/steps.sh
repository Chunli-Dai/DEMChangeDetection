#/bin/sh -f

tilelist=('arcticdem_01_iceland' 'arcticdem_02_greenland_southeast' 'arcticdem_03_greenland_southwest' 'arcticdem_04_greenland_central' 'arcticdem_05_greenland_northeast' 'arcticdem_06_greenland_northwest' 'arcticdem_07_canada_ellesmere' 'arcticdem_08_canada_baffin' 'arcticdem_09_canada_victoria' 'arcticdem_10_canada_north_mainland' 'arcticdem_11_canada_north_hudson' 'arcticdem_12_canada_south_nwt' 'arcticdem_14_svalbard' 'arcticdem_15_russia_novaya_zemlya' 'arcticdem_18_russia_cherskly' 'arcticdem_19_russia_magadanskaya' 'arcticdem_20_russia_kamchatka' 'arcticdem_21_russia_yakutiya_east' 'arcticdem_22_russia_central_east' 'arcticdem_23_russia_yakutiya_west' 'arcticdem_24_russia_central_west' 'arcticdem_25_russia_norilsk' 'arcticdem_26_russia_petersburg' 'arcticdem_27_russia_murmansk' 'arcticdem_28_scandinavia' 'arcticdem_29_russia_franz_josef' 'arcticdem_30_russia_siberian_islands' 'arcticdem_31_alaska_south' 'arcticdem_34_alaska_north')

nlist=${#tilelist[@]}
echo Total number of tiles: $nlist.

echo ${tilelist[*]}

for (( i=14; i<=$nlist; i++ ))
do
dir=${tilelist[$i-1]}
echo $i $dir

#dir=arcticdem_12_canada_south_nwt
cd $dir
cp /home/chunlidai/blue/apps/landslide/template/* .   # warning: need to change the directory
sed -i 's/arcticdem_08_canada_baffin/'$dir'/g' constant.m #change directory
grep $dir ademtiles_join-ademreg.csv > df
sed -i 's/'$dir',/''/g' df 
mv df tilelist  #get tilelist
#Get mat0.mat;delete old  mat0.mat striplist.dat
sbatch job.pbs #run Tilemain_nov.m to update the arcticdem_nov.tif and get mat0.mat.

run ./compile.sh to compile all matlab codes.
cp /home/chunlidai/blue/apps/landslide/template/run_Tilemain.sh . #add LD_LIBRARY_PATH to run_Tilemain.sh 

#bundle all jobs
nohup ./run_change_group_pgc_par.sh > outrun1 &
#nohup ./run_change_group_pgc_bundle.sh > out2_run1 &

cd ../
#exit
done

