
%Step 1: Prepare training images from polygon shapefiles and elevation change files
% Output: jpg images and json files that are compatible with Detectron2.

%To do manually: ln -fs /home/chunli/chunliproject/work/landslide/arcticdem_19_russia_magadanskaya/mat0.mat .
%		 Update regiondir in constant.m
%		 check resr2m=2; in box2mosaic.m; note when slumps area too big, this may cause out of memory error. Use 10 m in that case.

constant
addpath(genpath([codedir]));

%Banks
imagefiles={'/Users/chunlidai/Downloads/banksnew/merge_jump.tif'};
%shapefiles={'Toni_merge_polygons.shp'}; %true slumps negatvie area
shapefiles={'manualperfect/perfect_merge_polygons.shp'}; %true slumps negatvie area
%shapefiles={'slumpsbanks_band1vl_1p_fixr_np_shp/train.shp'}; %matched with positive clusters (debris).
%shapefiles={'change2polygons_russia.shp'};
% 11 regions
shapefiles={'traindata/perfect_merge_polygons.shp'};
%shapefiles={'traindata/change2polyg_07_ai.shp'};
%shapefiles={'traindata/change2polyg_08_m_ai.shp'};
shapefiles={'traindata/change2polyg_10_ai_wopeel.shp'}; %shapefiles={'traindata/change2polyg_10_ai.shp'};
%shapefiles={'traindata/change2polyg_11_m_ai.shp'};
%shapefiles={'traindata/change2polyg_15_ai.shp'}; % check 47 (deleted), 12090041.out; 
%shapefiles={'traindata/change2polyg_18_ai.shp'};
%shapefiles={'traindata/change2polyg_21_m_ai_cri1to3.shp'};
%shapefiles={'traindata/change2polyg_25_m_ai.shp'}; %failed (due to missing shx files)
%shapefiles={'traindata/change2polyg_27_ai.shp'};
%shapefiles={'traindata/change2polyg_30_ai.shp'};

%Ingmar polygons
if 0
shapefiles={'merged_thaw_slumps.shp','Slumps_2017_18_SCAR_FINAL.shp'}; %raw ground truth 
shapefiles={'change2polygons_ingmar.shp'};
end

if 0
%[trainingImages1,trainingLabels1,slumpbox]=prepareDetectron2_getnp(shapefiles,imagefiles);
%for svm
[featlabels,featg,category_id_feat]=prepareSVM(shapefiles,imagefiles);

save test2_scenario2.mat -v7.3
end


[trainingImages1,trainingLabels1,slumpbox]=prepareDetectron2(shapefiles,imagefiles);

