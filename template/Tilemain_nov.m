% main program for getting coastline for each ArcticDEM tile
% Requirements: gdal software
% %%%% inputs needed

% collect the nov for all area /home/chunli/chunliwork/work/landslide/arcticdem_nov.tif

constant 

macdir='/Users/chunlidai/surge/';
macdir=[];

currentdir=pwd;
%addpath(genpath(currentdir));
%addpath(genpath([macdir,'/data/chunli/landslide/']));
%addpath(genpath([macdir,'/data/chunli/scripts/']));
%addpath(genpath([macdir,'/home/dai.56/arcticdemapp/landslide/code1/']));
addpath(genpath([codedir]));

shpname='./GSHHS/GSHHS_f_L1.shp';% a priori coastline shapefile

%General directory that contains the tile DEM files, such as /elev/dem/setsm/ArcticDEM/mosaic/v2.0/
%tiledir='/Users/chunlidai/surge/data/chunli/coastline/';%ArcticDEM mosaic tile directory. 
%tiledir=[macdir,'/fs/byo/howat-data3/ArcticDEMmosaics/']; %
%stripdir='/*/ArcticDEM/region*/strips/2m/';
%stripdir='/fs/byo/howat-data2/ArcticDEM/region*/strips/2m/';
%stripdir='/fs/project/howat.4/EarthDEM/region*/strips_unf/2m/';
%stripdir=currentdir;

% %%%% control parameters
width=2e3; %buffer width of the a priori coastline, e.g., 2km.

if ~exist('mat0.mat','file') %readding boundary; time 1 hour for 90751 strip files and 10130 mono xml files

% %%% Preparation: get the list of strip files and boundries
filename='boundaries_regall_strip.dat'; %'boundaries_reg31.dat';
filename='striplist.dat';
if ~exist(filename,'file')
   %str=sprintf('find  %s -name ''*mdf.txt'' > %s',deblank(stripdir),filename);
   str=sprintf('find  %s -name ''*meta.txt'' > %s',deblank(stripdir),filename);
  [status, cmdout]=system(str);
end
fprintf ('\n Step 0: geting the boundary for all files in the region.\n')
%READ INPUT PARAMETERS; getting the boundaries for all files
% filename='boundaries_reg31.dat';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
range=zeros(n,4);XYbg=cell(n,1);
for i=1:n
   ifile=[macdir,fgetl(fid)];
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
   satname=f{i}(1:4);

   % get the boundary from xml file
   [XYbi,rangei]=imagebd(ifile);
   range(i,1:4)=rangei;XYbg{i}=XYbi;
end

save mat0.mat -v7.3

else 
load mat0.mat
end

dx=100e3;x0=-4000e3;y0=-4000e3;
dxs=dx/2/10; %5km
dxs=dx/50; %2km
xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;

%/home/chunli/chunliwork/work/landslide/arcticdem_nov.tif
%ofile='/home/chunli/chunliwork/work/landslide/arcticdem_nov.tif';
ofile='../arcticdem_nov.tif';
if exist(ofile,'file')
nov=readGeotiff(ofile);
nx=length(nov.x);ny=length(nov.y);
else %40m resolution
% integer type, 0 for nan.
resrc=400;
nov.x=x0:resrc:xe;
nov.y=ye:(-resrc):y0;
nx=length(nov.x);ny=length(nov.y);
nov.z=uint16(zeros(ny,nx));
end

resrc=mean(diff(nov.x));
%M=logical(size(nov.z));

novt=nov; % for this region
novt.z=uint16(zeros(ny,nx));

for i=1:n
        XYbi=XYbg{i};
        Xb=XYbi(:,1);Yb=XYbi(:,2);

        idx=round((Xb-nov.x(1))/resrc)+1;
        idy=round((Yb-nov.y(1))/(-resrc))+1;
        Mb = poly2mask(idx,idy, ny,nx); % build polygon mask       
	novt.z=novt.z + uint16(Mb);
end

%update nov with the new novt when nov is zero and novt is larger.
M=novt.z>nov.z;
nov.z(M)=novt.z(M);

projstr='polar stereo north';
writeGeotiff(ofile,nov.x,nov.y,uint16(nov.z),12,255,projstr)


