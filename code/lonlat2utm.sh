#/bin/sh -f

for line in `ls n*v3.tif`; do
echo $line  20000211000000_SRTM_utm_$line 
#gdalwarp 20000211000000_srtm.tif 20000211000000_SRTM_utm.tif -t_srs epsg:3413 -tr 30 30
gdalwarp $line  20000211000000_SRTM_utm_$line -t_srs epsg:32645 -tr 30 30
done


