
addpath(genpath('/home/chunlidai/blue/apps/landslide/code1/'));

if 1
%Merge multiple training images and json files
indir={'./slumpsIngmar/','../banks/pgcnew/slumpsbanks_band1val/','../ashleytuk/slumps/'};
indir={'./slumps/','../ashleytuk/slumps/'};
odir='./slumpsmerged/';
indir={'./slumps/','../ashleytuk/slumps/','../banks/pgcnew/slumpsbanks_band1val/'};
odir='./slumpsmergedBanks/';

indir={'./slumps/','../../banks/pgcnew/slumpsbanks_band1vl_mp_fixr/'};
odir='./slumpsbanks_band1vl_mp_fixr_2class/';

indir={'./slumps/','../../banks/pgcnew/slumpsbanks_band1vl_mp_wofixr_curv/'};

% 11 regions
indir={'./slumps_region09/','./slumps_region07/','./slumps_region08/','./slumps_region10/','./slumps_region11/','./slumps_region15/','./slumps_region18/','./slumps_region21/','./slumps_region25/','./slumps_region27/','./slumps_region30/'};
odir='./slumps_com/';
% augmentation
%indir={'./slumpsregion09aug/','./slumpsregion07aug/','./slumpsregion08aug/','./slumpsregion10aug/','./slumpsregion11aug/','./slumpsregion15aug/','./slumpsregion18aug/','./slumpsregion21aug/','./slumpsregion25aug/','./slumpsregion27aug/','./slumpsregion30aug/'};
%odir='./slumps_comaug/';
[co]=mergejson(indir,odir);
end

% Do this first for each group
% remove bad images from json files
if 0
odir='./slumps/';
odir='./slumpsbanks_band1vl_mp_wofixr_curv_selectionb2/';
odir='./slumpsbanks_band1vl_mp_wofixr_curv_rot_selectionb/';
jsonfile_t=[odir,'/train/via_region_data.json'];
jsonfile_val=[odir,'/val/via_region_data.json'];

%badfiles={'slump2g1.jpg','slump4g1.jpg'};
%Tuk Hugues
badfiles={'slump976.jpg','slump806.jpg','slump926.jpg','slump927.jpg','slump924.jpg','slump111.jpg','slump108.jpg','slump666.jpg','slump668.jpg','slump656.jpg','slump654.jpg','slump55.jpg','slump551.jpg','slump253.jpg','slump326.jpg','slump465.jpg','slump470.jpg','slump252.jpg','slump243.jpg','slump1.jpg'};

%Ingmar
badfiles={'slump94.jpg','slump95.jpg','slump84.jpg','slump178.jpg','slump188.jpg','slump191.jpg','slump187.jpg','slump182.jpg','slump201.jpg','slump200.jpg','slump213.jpg','slump214.jpg','slump208.jpg','slump250.jpg','slump262.jpg','slump272.jpg','slump256.jpg','slump235.jpg','slump197.jpg','slump34.jpg','slump25.jpg','slump13.jpg','slump11.jpg','slump314.jpg','slump317.jpg','slump318.jpg','slump319.jpg','slump288.jpg','slump287.jpg','slump487.jpg','slump491.jpg','slump490.jpg','slump492.jpg','slump493.jpg','slump494.jpg','slump495.jpg','slump496.jpg','slump497.jpg'};

%banks
filename='badfilelistb.txt';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
%badfiles=cell(n,1);
count=0;
for i=1:n
   ifile=[fgetl(fid)];
   %badfiles{i}=ifile;
   count=count+1;
   badfiles{count}=ifile;
	dirfilename=ifile
   for iaug=1:100
	   count=count+1;

        strg{5}=['rotation',num2str(iaug)];
        strcom=[strg{5},'.jpg'];
        dirfilename_aug=strrep(dirfilename,'.jpg',strcom);
	badfiles{count}=dirfilename_aug;
   end
end

[co]=modifyjson(jsonfile_t,badfiles);
[co]=modifyjson(jsonfile_val,badfiles);
end
