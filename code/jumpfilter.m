%Given elevation change; filter bad points and then fill voids.
%Try to make sure that data within the lava flow area has no voids.

addpath(genpath(['/home/dai.56/arcticdemapp/landslide/code1/']));

jump=readGeotiff('givensite_jump.tif');
jumpstd=readGeotiff('givensite_jumpstd.tif');

%remove points with std >10 
M=jump.z==0|jumpstd.z>=10; %(A*(A != 0 & ~( B > 10)))
jump.z(M)=nan; % bad/void to nan

projstr='polar stereo north';
OutName=['givensite_jumpfilter1.tif'];
writeGeotiff(OutName,jump.x,jump.y,double(jump.z),4,nan,projstr); 

%use gdal to fill voids up to 15 pixels in the neighbor
str='gdal_fillnodata.py -md 15 givensite_jumpfilter1.tif givensite_jumpfilter1filled.tif';
[status, cmdout]=system([str])

