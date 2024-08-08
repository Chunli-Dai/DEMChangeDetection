
constant
codedir
addpath(genpath([codedir]));
%addpath(genpath('/home/dai.56/arcticdemapp/landslide/code1/mlearning/'));

pointfiles='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Toni.gmt';
imagefiles={'/Users/chunlidai/Downloads/banksnew/36_21_2_1_jump.tif'};
% [shapefiles]=changemap2polygon(pointfiles,imagefiles);

imagefiles={'/Users/chunlidai/Downloads/banksnew/merge_jump.tif'};

imagefilesall={'37_18_2_2_jump.tif' '37_19_1_2_jump.tif' '41_16_2_1_jump.tif' '41_17_1_1_jump.tif' '42_63_2_1_jump.tif' '43_59_1_2_jump.tif' '48_59_2_2_jump.tif' '55_46_2_1_jump.tif' '55_46_2_2_jump.tif' '60_42_2_2_jump.tif' '63_44_1_2_jump.tif' '63_44_2_2_jump.tif' '63_45_1_1_jump.tif' '63_45_1_2_jump.tif' '64_42_2_2_jump.tif' '64_44_1_1_jump.tif' '64_44_2_1_jump.tif' '64_44_2_2_jump.tif' '65_44_2_1_jump.tif' '65_44_2_2_jump.tif' '65_45_1_2_jump.tif' '66_44_2_1_jump.tif' '66_44_2_2_jump.tif' '67_45_1_1_jump.tif'};

%ls ../arcticdem_18_russia_cherskly/*_jump.tif > jumplist.txt
filename='jumplist.txt';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
for i=1:n
   ifile=[fgetl(fid)];
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   imagefilesall{i}=[demdir,'/',name,ext];
end

%pointfiles={'merged_thaw_slumps.shp','Slumps_2017_18_SCAR_FINAL.shp'}; %raw ground truth 
pointfiles='RTS_NESiberia.gmt'; %points one string file; polyons can be multiple files.
n=length(imagefilesall);
fprintf(['\n Total number of tiles to work on ',num2str(n) ,' .\n'])

shp1g=[];
shapefiles=['change2polygon_tilei.shp'];
shapefiles_del=strrep(shapefiles,'.shp','.*');
str1=['rm -f ',shapefiles_del];
[status, cmdout]=system(str1);
for i=1:n
    %imagefiles={['/Users/chunlidai/Downloads/banksnew/',imagefilesall{i}]};
    imagefiles={['./',imagefilesall{i}]};
    fprintf(['\n Working on i: ',num2str(i),' ',imagefiles{1}])
    tic
    [shapefiles]=changemap2polygon(pointfiles,imagefiles);
    if exist(shapefiles,'file')
        shp1{i}=shaperead(shapefiles);
        shp1g=[shp1g;shp1{i}];

	%delete shapefiles to avoid repeat loading
	shapefiles_del=strrep(shapefiles,'.shp','.*');
	str1=['rm -f ',shapefiles_del];
	[status, cmdout]=system(str1);
    end
    toc
end
%merge shapefiles
%ofile1='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Toni_merge_polygons.shp';
%ofile1='change2polygons_ingmar.shp';
ofile1='change2polygons_russia.shp';
shapewrite(shp1g, ofile1);
