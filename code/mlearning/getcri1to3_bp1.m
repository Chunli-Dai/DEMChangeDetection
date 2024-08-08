constant
codedir='/home/dai.56/arcticdemapp/landslide/code1/';  %Directory of codes.

addpath(genpath([codedir]));

load test1.mat

%[co]=system('gdal_merge.py -o 41_16_1_1_10m_dem.tif slumps_dem/*.tif -ps 10 10 -n -9999 -a_nodata nan');
%dem=readGeotiff('41_16_1_1_10m_dem.tif'); %2m
dem=readGeotiff('41_16_1_1_4m_dem.tif'); %2m
demr.z=interp2(dem.x,dem.y,double(dem.z),jump.x,jump.y','nearest',nan);
demr.x=jump.x;demr.y=jump.y;

if 0
[maskcri1to3]=dhelev(jump,demr,[],S2,Mn,Mp);
save maskcri1to3.mat maskcri1to3 jump -v7.3
end

dem=demr;

buff=3000;
x0y0=[-2498790, 3984];
x0y0=[-2495448, 9967];
rangi=[x0y0(1)-buff x0y0(1)+buff x0y0(2)-buff x0y0(2)+buff ];
Mn_sub=cropmatrix(Mn,rangi,2);
Mp_sub=cropmatrix(Mp,rangi,2);
jump_sub=cropmatrix(jump,rangi,2);
dem_sub=cropmatrix(dem,rangi,2);
%load test2.mat
[maskcri1to3_sub]=dhelev(jump_sub,dem_sub,[],S2,Mn_sub,Mp_sub);
%save maskcri1to3_sub_u.mat maskcri1to3_sub -v7.3
save test2_u.mat jump_sub dem_sub S2 Mn_sub Mp_sub maskcri1to3_sub -v7.3

