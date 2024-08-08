How to run the code:
Step 1: copy Tilemain.m, constant.m, job_group.slurm, and run_change_group_slurm.sh to your work directory.
Step 2: run ./run_change_group_slurm.sh

Parameters need to be changed:
1\ In constant.m, change stripdir and codedir.
e.g. stripdir='/fs/project/howat.4/EarthDEM/region*/strips_unf/2m/';
         where stripdir is the directory of EarthDEM data.
2\ In run_change_group.sh
Modify the xid and yid to the desired tile that you'd like to run. e.g.
if you'd like to run code for tile 41_16_1_1, let yid=41;xid=16;xids=1;yids=1. 

Data/Software preparation:
1\ setsm software should be installed and executable in the command line. 

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


