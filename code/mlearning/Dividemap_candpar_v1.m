constant

%New versoin May 2023: Use the given candidates from Change2poly6classespar.m to get test images.
%new version Jan, 2022: use constraint 1 to get candidate clusters, and then get images for all clusters.

%codedir='/home/dai.56/arcticdemapp/landslide/code1/';  %Directory of codes.
flagorthotif=0; % save mosacked ortho image in Geotiff
flagdemtif=0; % save mosacked dem in Geotiff
flagrescale=0;  %1, fixed range [e.g. -20 10] to keep the information of positive/negative.
                %0, re-adjust the range for each image.
flagband=3; %3; %1 writing 1 band data for training images; 3 writing 3 band data; 2: elevation change and DEM curvature, two bands.

[demdir,name,ext] =fileparts([strtrim(regiondir)]);
regid=name(11:12);

%regionstr='eureka';
%regionstr='russia';
%regionstr='peel';

regionstr=['region',regid]; %'region09';
shapefile=['change2polyg_',regid,'.shp';] %'../region09/v3con/change2polyg_09.shp'; %candidates
mergedfile=[regiondir,'/region',regid,'_jump.tif'];

if 1 %flagorthotif==1
        %save orthoimage in geotiff 
        odirt1=['slumps_ortho/'];
        if ~exist(odirt1,'')
            mkdir(odirt1)
        end
end
odirt2='slumps_dem/';
%mkdir slumpsub_3bandsuint8jpg
%mkdir slumpsub_band1changeuint8jpg
%mkdir slumpsub_band2demuint8jpg
%mkdir slumpsub_band3orthouint8jpg

addpath(genpath([codedir]));

%get fdir f
matfile=[regiondir,'/mat0.mat'];
load(matfile)
%load mat0.mat

%% Read the candidate shapefile (output from Change2poly6classespar.m)
% % Old: Divide target image to ~1 km size;
S1=shaperead(shapefile);
shpc=struct2cell(S1)';
areaj=[shpc{:,6}]';
M1=areaj>=2000;
S1s=S1(M1);S1=S1s;
ns=length(S1);

imagefiles_val={mergedfile};
%imagefiles_val={'32_33_2_2_jump.tif'}; %eureka
%imagefiles_val={'49_60_2_2_jump.tif'}; %nina
changefile=mergedfile; %imagefiles_val{1};

%to get the grids
data0=readGeotiff(changefile);
data0.z=[]; %clear data in data0.z
resrx=mean(diff(data0.x));resry=mean(diff(data0.y));
%make sure the data0.z has the resolution of 10 m; 
%100m to 10m;
if resrx~=10;
resrc=10 ; %40;
rang0=[min(data0.x) max(data0.x) min(data0.y) max(data0.y) ]; 

ranget=round(rang0/resrc)*resrc;rang0=ranget;
tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);
xoutr=tx;youtr=ty;
data0.x=xoutr;data0.y=youtr;
%data0.z=false(length(data0.y),length(data0.x));
end

%data0.z=false(size(data0.z));
resrx=mean(diff(data0.x));resry=mean(diff(data0.y));

%step: 500 m ; buffer 200 m on each size. for Peel Plateau
% try 800 pixels (2m ): 1600 m ; For Eureka
%buff=100; dx=300; %-> Chandi uses 250 pixels by 250 pixels;
%dy=-dx;

buff=100; %buffer the polygon by 100m % consistent with training images prepareDetectron2.m

odir1=['slumpsub_',regionstr,'/']; 
%odir1='slumpsub_peel_timeband3/';
if ~exist(odir1,'')
   mkdir(odir1)
end


fprintf(['There are a total of ',num2str(ns),' candidate clusters.\n'])

%[X,Y]=meshgrid(data0.x,data0.y);

%ids={};
ids=cell(ns,1);Mpcg=cell(ns,1);

sz = getenv('SLURM_NTASKS');
sz=str2num(sz);
fprintf(['\n ',num2str(sz),' worker(s) allocated in job.slurm.\n'])
poolobj=parpool(sz);
%poolsize=20;
%poolobj=parpool(39); %hi
%poolobj=parpool(190);
%poolobj=parpool(119);  %Error: Failed to bind to endpoint; This is probably the following known problem in R2023a;
%poolobj=parpool(100); %success
parfor ix=1:ns
%for ix=1:ns

        fprintf(['\n Working on ',num2str(ix),'/',num2str(ns),' candidate cluster.\n'])
	tic

        [Mpc,idsi]=dividemapsub(ix,S1,buff,data0,flagrescale,flagband,f,fdir,odir1);
	Mpcg{ix}=Mpc;

        ids{ix}=idsi;
	fprintf(['\n ',num2str(ix),' takes ',num2str(toc),' seconds.\n']);
end
delete(poolobj)

save test1.mat -v7.3

%Move this to later to save memory issue.
for ix=1:0 %ns
        Mpc=Mpcg{ix};
        Mx=Mpc.Mx;My=Mpc.My;
        
        Mpcpre.z=data0.z(My,Mx); 

        %data0.z(My,Mx)=Mpc.z|Mpcpre.z;
	tmp=Mpc.z|Mpcpre.z;
        data0.z(My,Mx)=tmp;

end

%save ixy.mat data0 ids buff -v7.3
%save ixy_eureka.mat data0 ids buff dx rangx rangy ns -v7.3
ofile=['ixy_',regionstr,'_cand.mat'];
save(ofile,'data0', 'ids', 'ns','Mpcg', '-v7.3')
