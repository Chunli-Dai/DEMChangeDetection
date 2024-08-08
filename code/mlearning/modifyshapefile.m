codedir='/home/dai.56/arcticdemapp/landslide/code1/';  %Directory of codes.

addpath(genpath([codedir]));


% Modify shapefiles from Jurjen
% 
%output
o_shapefile={'/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Jurjen_Segal_series/slumpfromchange_JS/maskgroundtruth_Jurjen_v2.shp'}; 
o_shapefile={'./Jurjen/maskgroundtruth_Jurjen_v2.shp'}; 
%version 2; class noise is the same; map individual gullies in gully area; 
%select scar with elevation change <-1m, and match the corresponding positive area;

%input
%Sj=shaperead('/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Jurjen_Segal_series/slumpfromchange_JS/maskgroundtruth.shp');
Sj=shaperead('Jurjen/maskgroundtruth.shp');
for i=1:length(Sj)
    if strcmp(Sj(i).zone,'scar');Mscar(i)=1;
    elseif strcmp(Sj(i).zone,'nois');Mscar(i)=2;
    elseif strcmp(Sj(i).zone,'gull');Mscar(i)=3;
    elseif strcmp(Sj(i).zone,'quar');Mscar(i)=4;%1 area
    elseif strcmp(Sj(i).zone,'debr');Mscar(i)=5;%2 areas
    elseif strcmp(Sj(i).zone,'nogo');Mscar(i)=6;%1 area
    end
end
figure;histogram(Mscar)

if 0
% 484 scars (true positives) by Jurjen;but lots of them have no elevation change signal
Sj_tp=Sj(Mscar==1);
ofile1='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Jurjen_Segal_series/slumpfromchange_JS/maskgroundtruth_scar.shp';
shapewrite(Sj_tp, ofile1);
Sj_tp=Sj(Mscar==2);
ofile1='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Jurjen_Segal_series/slumpfromchange_JS/maskgroundtruth_nois.shp';
shapewrite(Sj_tp, ofile1);
Sj_tp=Sj(Mscar==3);
ofile1='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Jurjen_Segal_series/slumpfromchange_JS/maskgroundtruth_gull.shp';
shapewrite(Sj_tp, ofile1);
end

%% noise; keep it the same.
Sj_noise=Sj(Mscar==2);


%semi-manually get the ground truth based on field data and elevation change data
%data0=readGeotiff('/Users/chunlidai/Downloads/41_16_1_1_jump.tif');
data0=readGeotiff('41_16_1_1_jump.tif');
smlarea=1e3*2; %1000 m^2
resx=mean(data0.x(2:end)-data0.x(1:end-1));resy=mean(data0.y(2:end)-data0.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);
%Get the examed slump area (as ground truth).
Mf4=data0.z<=-1; %abs(jump)>=2; %jump significant;Only negative
Mf4clean = bwareaopen(Mf4, smlarea/resr/resr); %remove small clusters
data0x=data0.x;data0y=data0.y;

Mp1=data0.z>0.5; %>1m;
Mpclean = bwareaopen(Mp1, smlarea/resr/resr/2); %remove small clusters
clear Mp
Mp.x=data0.x;Mp.y=data0.y;Mp.z=Mpclean;

%% scars
if 0
Sj_scar=Sj(Mscar==1);

% Get the mask drawn manually that contains slumps.
% S1=shaperead('/Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/peelmanual.shp');% figure;mapshow(S1)
S1=Sj_scar;% figure;mapshow(S1)
n=length(S1);
[ny,nx]=size(data0.z);
mask=false(size(data0.z));
for i=1:n
        [s1x,s1y]=polarstereo_fwd(S1(i).Y,S1(i).X,[],[],70,-45);
        s1x(isnan(s1x))=[];s1y(isnan(s1y))=[];
        Xb=round((s1x-data0.x(1))/resr)+1; 
        Yb=round((s1y-data0.y(1))/(-resr))+1; %descending    
        maski=poly2mask(Xb,Yb, ny,nx); 
        mask=mask|maski;
end
%add 50 m buffer, since some of the shapefiles cut off the slumps edges.
mask=(imdilate(mask,ones(round(50*2/resr))));

Mf5=Mf4clean&mask;
Mf5clean = bwareaopen(Mf5, smlarea/resr/resr); %remove small clusters

%find the matched positive clusters.
Mn=data0;Mn.z=Mf5clean;

fprintf('\n Start matching.\n')
tic
[Mouts, Mout]=matchclusters(Mn,Mp); %400 clusters, 15 minutes
toc 

%Generate polygon shape files.
Mf5clean=Mout.z;
Modfil=~(Mf5clean);  %Modfil (logical; 1 No change, 0 Change) 
Medgs1=false(size(Mf5clean));
[Sj_scar_np]=mask2boundary(data0x,data0y,(Modfil),Medgs1,'t1.shp');
Sj_scar_np.id=1:length(Sj_scar_np);
shapewrite(Sj_scar_np, 't1.shp'); %lots of bad ones.
end %if 0

%% gullies
Sj_gully=Sj(Mscar==3);

% Get the mask drawn manually that contains slumps.
S1=Sj_gully;% figure;mapshow(S1)
n=length(S1);
[ny,nx]=size(data0.z);
mask=false(size(data0.z));
for i=1:n
        [s1x,s1y]=polarstereo_fwd(S1(i).Y,S1(i).X,[],[],70,-45);
        s1x(isnan(s1x))=[];s1y(isnan(s1y))=[];
        Xb=round((s1x-data0.x(1))/resr)+1; 
        Yb=round((s1y-data0.y(1))/(-resr))+1; %descending    
        maski=poly2mask(Xb,Yb, ny,nx); 
        mask=mask|maski;
end

Mf5=Mf4clean&mask;
Mf5clean = bwareaopen(Mf5, smlarea/resr/resr); %remove small clusters

%find the matched positive clusters.
Mn2=data0;Mn2.z=Mf5clean;

 CCp = bwconncomp(Mn2.z);
 fprintf(['\n Number of negative clusters within gullies:',num2str(CCp.NumObjects)])

fprintf('\n Start matching.\n')
tic
[Mouts, Mout2]=matchclusters(Mn2,Mp);
toc
save test1.mat -v7.3

%Generate polygon shape files.
Mf5clean=Mout2.z;
Modfil=~(Mf5clean);  %Modfil (logical; 1 No change, 0 Change) 
Medgs1=false(size(Mf5clean));
[Sj_gully_np]=mask2boundary(data0x,data0y,(Modfil),Medgs1,'t2.shp');

save test1.mat -v7.3
 
