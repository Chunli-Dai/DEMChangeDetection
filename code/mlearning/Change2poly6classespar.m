
constant
%addpath(genpath(['/home/chunlidai/blue/apps/landslide/code1/']));
% %addpath(genpath('/home/dai.56/arcticdemapp/landslide/code1/mlearning/'));

%pointfiles='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Toni.gmt';
% [shapefiles]=changemap2polygon(pointfiles,imagefiles);

%ls ../arcticdem_18_russia_cherskly/*_jump.tif > jumplist.txt
%regiondir='/blue/chunlidai/chunlidai/landslide/arcticdem_10_canada_north_mainland/';
[demdir,name,ext] =fileparts([strtrim(regiondir)]);
regid=name(11:12);
indir=name; %'arcticdem_10_canada_victoria';
ofile1=['change2polyg_',regid,'.shp']; %'change2polyg_09.shp';
if ~exist('jumplist.txt','file')
%str1=['ls ../../',strtrim(indir),'/mergetiles/*_jump.tif > jumplist.txt '];
str1=['ls ',strtrim(regiondir),'/mergetiles/*_jump.tif > jumplist.txt '];
[status, cmdout]=system(str1);
end

imagefilesall={};
filename='jumplist.txt';
pointfiles='';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
for i=1:n
   ifile=[fgetl(fid)];
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   imagefilesall{i}=[demdir,'/',name,ext];
end

n=length(imagefilesall);
fprintf(['\n Total number of tiles to work on ',num2str(n) ,' .\n'])

shp1g=[];
shapefiles=['change2polygon_tilei.shp'];
shapefiles_del=strrep(shapefiles,'.shp','.*');
str1=['rm -f ',shapefiles_del];
[status, cmdout]=system(str1);

shp1=cell(n,1);

sz = getenv('SLURM_NTASKS');
sz=str2num(sz);
fprintf(['\n ',num2str(sz),' worker(s) allocated in job.slurm.\n'])
poolobj=parpool(sz);
%poolobj=parpool(180); %34
%poolobj=parpool(64); %34
parfor i=1:n
    %imagefiles={['/Users/chunlidai/Downloads/banksnew/',imagefilesall{i}]};
    %imagefiles={['./',imagefilesall{i}]};
    imagefiles={[imagefilesall{i}]};
    fprintf(['\n Working on i: ',num2str(i),' ',imagefiles{1}])
    [shapefiles]=changemap2polygon(pointfiles,imagefiles);
     shp1{i}=shapefiles;

end
delete(poolobj)

save test1.mat -v7.3

for i=1:n
shp1g=[shp1g;shp1{i}];
end

shpc=struct2cell(shp1g)';
areaj=[shpc{:,6}]';

%merge shapefiles
%ofile1='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Toni_merge_polygons.shp';
%ofile1='change2polygons_ingmar.shp';
%ofile1='change2polyg_09.shp';
shapewrite(shp1g, ofile1);

%6 classes,
% 6\ [500,000 m2 <= a ) landslide; very very large
% 5\ [100,000 <= a < 500,000) very large
% 4\ [50,000 <= a < 100,000) large
% 3\ [10,000 <= a < 50,000) median
% 2\ [5,000 <= a < 10,000)   small
% 1\ [1,000 <= a < 5,000)   very small
ofilevs=strrep(ofile1,'.shp','_vs.shp');
ofiles=strrep(ofile1,'.shp','_s.shp');
ofilem =strrep(ofile1,'.shp','_m.shp');
ofilel =strrep(ofile1,'.shp','_l.shp');
ofilevl=strrep(ofile1,'.shp','_vl.shp');
ofilela=strrep(ofile1,'.shp','_landslide.shp');

nj=length(shp1g);
%shpvs=[];shps=[];shpm=[];shpl=[];shpvl=[];shpla=[];
M1=areaj>=2000&areaj<5000;
shpvs=[shp1g(M1)];
M1=areaj >=5000 & areaj<10e3;
shps=[shp1g(M1)];
M1=areaj >=10000 & areaj<50e3;
shpm=[shp1g(M1)];
M1=areaj >=50000 & areaj<100e3;
shpl=[shp1g(M1)];
M1=areaj >=100000 & areaj<500e3;
shpvl=[shp1g(M1)];
M1=areaj >=500e3;
shpla=[shp1g(M1)];

shapewrite(shpvs, ofilevs);
shapewrite(shps, ofiles);
shapewrite(shpm, ofilem);
shapewrite(shpl, ofilel);
shapewrite(shpvl, ofilevl);
shapewrite(shpla, ofilela);

nla=length(shpla); nvl=length(shpvl);nl=length(shpl);nm=length(shpm);ns=length(shps);nvs=length(shpvs);
nall=nla+nvl+nl+nm+ns+nvs;
fprintf(['There are a total of ',num2str(nla),' super large clusters.\n'])
fprintf(['There are a total of ',num2str(nvl),' very large clusters.\n'])
fprintf(['There are a total of ',num2str(nl),' large clusters.\n'])
fprintf(['There are a total of ',num2str(nm),' median clusters.\n'])
fprintf(['There are a total of ',num2str(ns),' small clusters.\n'])
fprintf(['There are a total of ',num2str(nvs),' very small clusters.\n'])
fprintf(['There are a total of ',num2str(nall),' all types clusters.\n'])

