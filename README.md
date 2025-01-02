How to run the code:
Step 1: copy template/* to your work directory.
Step 2: Download example data (1GB) to your work directory and unzip it: https://drive.google.com/file/d/1G5VHwl334AbmMIsBwHYQciEgOHrpZh9N/view?usp=sharing 
Step 3: Update paramters in constant.m and prepare files tilelist or aoi.txt (provided for this example):
        In constant.m, make sure stripdir and codedir is correct.
Step 4: run nohup ./run_change_group_pgc_par.sh > outrun1 &

Data/Software requirement:
1\ Matlab license, and matlab installed on linux system as a module,
2\ gdal installed as a module,
3\ SLURM job management,
4\ (Optional) setsm software (https://github.com/setsmdeveloper/SETSM) may be installed and executable in the command line if you choose setsm for DEM coregistration. 

#######################
Output files: Elevation change data for an area size of 2 km by 2km.
It is suggested to run two scenarios. One is using summer (July August) data only, the other is using data from all seasons.

The result includes:
41_16_1_1_01_01_listused.txt: a list of ArcticDEM strip file names used for change estimation.
41_16_1_1_01_01_jump.tif: elevation change in meters.
41_16_1_1_01_01_jumpstd.tif: elevation change uncertainty in meters.
41_16_1_1_01_01_eventtimeT1.tif: date of the closest measurement before detected change. The date format is YYYYMMDD.
41_16_1_1_01_01_eventtimeT2.tif: date of the closest measurement after change. 
41_16_1_1_01_01_nov.tif: The number of overlapping DEMs used for each pixel.
41_16_1_1_01_01_bitmask.tif: Integer flag of whether the estimated change has high confidence level. 1 is good, and 0 is bad.

They are Geotiff files, which can be visualized in QGIS (open software). You can also read the files using Matlab code (https://github.com/ihowat/setsm_postprocessing/blob/master/readGeotiff.m).

The time data now includes the T1 (date of the closest measurement before change) and T2 (date of the closest measurement after change), and the date format is YYYYMMDD.

#######
References:


