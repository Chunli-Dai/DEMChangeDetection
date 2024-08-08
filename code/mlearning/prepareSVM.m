function [featlabels,featg,category_id_feat]=prepareSVM(shapefiles,imagefiles)
% Modified from prepareDetectron2.m 
% % Prepare training images for Detectron2 (upgraded MaskRCNN)
% Given shapefiles, and imagefiles to get training image matrices.
% shapefiles: shapefile of known slump outlines.
% imagefiles: elevation change files list
% output: slumpbox: the boundary of slumps

% trainingLabels1=categorical({'slump', 'non slump'}');

scenario=[];
featlabels={''};
constant
scenario

flagplot=0; %1 plot;0 not plot
flagaug=0; %if use data augment;
maxtrain=1725-173;%461; %917; %1725-173;%9e9; %maximum number of images used as training images, and rest will be used as validation.
flagorthotif=0; % save mosacked ortho image in Geotiff
flagrescale=1;  %1, fixed range [e.g. -20 10] to keep the information of positive/negative.
		%0, re-adjust the range for each image.

odir='slumps/';
%odir='slumpsbanks_band2dem/'; %test
%odir='slumpsbanks_band3ortho/'; %test
odir='slumpsbanks_3bands_uint8/';
odir='slumpsbanks_band1change_uint8/';
%odir='slumpsbanks_band2dem_uint8/';
%odir='slumpsbanks_band3ortho_uint8/';
odir='slumpsbanks_band1change_uint8jpg/';
odir='slumpsbanks_3bands_uint8jpg/'; %test hi
odir='slumpsbanks_band2dem_uint8jpg/';
odir='slumpsbanks_band3ortho_uint8jpg/';
odir='slumpsbanks_band1val/';
odir='slumpsbanks_band1vlaug/';
odir='slumps/'
odir1=[odir,'/train/'];odir2=[odir,'/val/'];odir3=[odir,'/test/'];
if ~exist(odir,'dir')
    mkdir(odir)
    mkdir(odir1)
    mkdir(odir2)
end
if ~exist(odir1,'dir')
    mkdir(odir1)
end
if ~exist(odir2,'dir')
    mkdir(odir2)
end
if ~exist(odir3,'dir')
    mkdir(odir3)
end

if flagorthotif==1
	%save orthoimage in geotiff 
 	odirt1=['slumps_ortho/'];	
	if ~exist(odirt1,'')
	    mkdir(odirt1)
	end
end
        

% S1=shaperead('/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Ashley/shapefiles/Slumps_2017_18_SCAR_FINAL.shp');
S1=shaperead(shapefiles{1});% figure;mapshow(S1)
%/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_06_21_jump.tif

flagband=3; %3; %1 writing 1 band data for training images; 3 writing 3 band data;
flagdata=2; %1 use the merged elevation change file; 2; use elevation change file in each sub tile.
% changefile='/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_merge.tif';
changefile=imagefiles{1};
if flagdata==1&&exist(changefile,'file') %use the merged file
data0=readGeotiff(changefile);
resr=mean(diff(data0.x));
else %use all subtiles
%construct a merged file; to be compatible with old version. Might to disable this.
rang0=[-2032100.000 -1627890.000 -638110.000 -299900.000];
resr=10; %
resr=100;
xout=rang0(1):resr:rang0(2);yout=rang0(4):-resr:rang0(3);
data0.x=xout;data0.y=yout;
ny=length(yout);nx=length(xout);
data0.z=zeros(ny,nx);;
end
resr2m=2; % 2m resolution;

%Generate mask for different classes: stable, slumps, noise, gullies based on given input shapefiles by manual classification.

f=[];fdir=[];
%get the f fdir list
load mat0.mat

%Get slumps
n=length(S1);
count=0; %clear trainingImages1 trainingLabels1 collecti
countaug=0; %count with data augmentation
xot={}; yot={}; fileot={}; 
xov={}; yov={}; fileov={}; 
xote={}; yote={}; fileote={}; %test
count_tr=0;count_v=0;count_te=0;
countf=0;
absfilenames={''};filenames=[];

buff=100; %buffer the polygon by 100m

%80 Different number of images yield incoherent mosaic.
%27 half of the polygon in old version(bp1). good in new version
%184 long shape
% 904 dark: anomaly in change map
for k=1:length(shapefiles)
S1=shaperead(shapefiles{k});
n=length(S1);

featgt_par=cell(n,1); category_id_feat_par=cell(n,1);
poolsize=20;
poolobj=parpool(20);
parfor i=1:n
%for i=1:n
	[featlabels,feati,classi]=prepareSVM_sub(S1,i,buff,shapefiles{k},data0,flagdata,flagband,f,fdir,scenario);
        featgt_par{i}=feati;
        category_id_feat_par{i}=classi;
	featlabels=featlabels;	
end %for i
delete(poolobj)

save test_svm2.mat -v7.3

%collect results
countf=0;
featgt=[];category_id_feat=[];
for i=1:n
%countf=countf+length(category_id_feat_par{i});
featgt=[featgt,featgt_par{i}];
category_id_feat=[category_id_feat,category_id_feat_par{i}];
end %i

end %for k
featg=featgt';


%write shapefiles of training images 
idd=find(cellfun(@isempty,xot)); %fix bug 12
xot(idd)=[];yot(idd)=[];fileot(idd)=[];
idd=find(cellfun(@isempty,xov)); %fix bug 12
xov(idd)=[];yov(idd)=[];fileov(idd)=[];
idd=find(cellfun(@isempty,xote)); %fix bug 12
xote(idd)=[];yote(idd)=[];fileote(idd)=[];

shp1 = struct('Geometry', 'PolyGon','X', xot, 'Y', yot,'filename',fileot);
ofile=[odir,'train.shp'];
if ~isempty(shp1)
shapewrite(shp1, ofile);
end

ofile=[odir,'val.shp'];
shp1 = struct('Geometry', 'PolyGon','X', xov, 'Y', yov,'filename',fileov);
if ~isempty(shp1)
shapewrite(shp1, ofile);
end

ofile=[odir,'test.shp'];
shp1 = struct('Geometry', 'PolyGon','X', xote, 'Y', yote,'filename',fileote);
if ~isempty(shp1)
shapewrite(shp1, ofile);
end

return
end
