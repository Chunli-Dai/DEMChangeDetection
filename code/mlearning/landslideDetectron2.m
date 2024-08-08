% https://github.com/facebookresearch/detectron2
% see https://colab.research.google.com/drive/1n_0KNE0XCu6AwgEfhwbr6tHuY_bp152d#scrollTo=PIbAM2pv-urF

addpath(genpath('/Users/chunlidai/Documents/MATLAB/Examples/R2020b/'));
% addpath(genpath('/Users/chunlidai/surge/home/dai.56/arcticdemapp/landslide/code1/'));
addpath(genpath('/Users/chunlidai/ArchivedBoxSync/ESI2019/landslide/code1'));

%Step 0: prepare training images from point shapefiles and elevation change files
% Output: polygon shapefiles containing the outlines of slumps.

%step0.m

pointfiles='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Toni.gmt';
imagefiles={'/Users/chunlidai/Downloads/banksnew/36_21_2_1_jump.tif'};
% [shapefiles]=changemap2polygon(pointfiles,imagefiles);

imagefiles={'/Users/chunlidai/Downloads/banksnew/merge_jump.tif'};
imagefilesall={'34_23_2_2_jump.tif' '34_24_1_2_jump.tif' '35_21_1_2_jump.tif' '35_21_2_2_jump.tif' '35_22_1_2_jump.tif' '35_22_2_2_jump.tif' '35_23_1_1_jump.tif' '35_23_1_2_jump.tif' '35_23_2_1_jump.tif' '35_23_2_2_jump.tif' '35_24_1_1_jump.tif' '36_20_2_2_jump.tif' '36_21_1_1_jump.tif' '36_21_1_2_jump.tif' '36_21_2_1_jump.tif' '36_21_2_2_jump.tif' '36_22_1_1_jump.tif' '36_22_2_1_jump.tif' '36_23_1_1_jump.tif' '36_23_2_1_jump.tif' '36_23_2_2_jump.tif' '36_24_1_1_jump.tif' '36_24_1_2_jump.tif' '36_24_2_2_jump.tif' '37_21_1_1_jump.tif' '37_21_2_1_jump.tif' '37_21_2_2_jump.tif' '37_22_1_1_jump.tif' '37_22_2_2_jump.tif' '37_23_2_2_jump.tif' '37_24_1_1_jump.tif' '37_24_1_2_jump.tif' '37_24_2_1_jump.tif'};
shp1g=[];
for i=30:length(imagefilesall)
    imagefiles={['/Users/chunlidai/Downloads/banksnew/',imagefilesall{i}]};
    [shapefiles]=changemap2polygon(pointfiles,imagefiles);
    if exist(shapefiles,'file')
    shp1{i}=shaperead(shapefiles);
    shp1g=[shp1g;shp1{i}];
    end
end
%merge shapefiles
ofile1='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Toni_merge_polygons.shp';
shapewrite(shp1g, ofile1);


%Step 1: Prepare training images from polygon shapefiles and elevation change files
% Output: jpg images and json files that are compatible with Detectron2.
% see step1.m step1b.m

% shapefiles='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Ashley/shapefiles/Slumps_2017_18_SCAR_FINAL.shp';
% imagefiles={'/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_merge.tif'};

% shapefiles='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Toni_polygons.shp';

imagefiles={'/Users/chunlidai/Downloads/banksnew/merge_jump.tif'};
shapefiles=ofile1;
[trainingImages1,trainingLabels1,slumpbox]=prepareDetectron2(shapefiles,imagefiles);
 
% Step 2: Run trainning on google colab: https://colab.research.google.com/drive/1n_0KNE0XCu6AwgEfhwbr6tHuY_bp152d#scrollTo=RF5qonX33Q05

% Step 3: run classification on validation images in Peel Plateau;

% Divide target image to ~1 km size;
% See dividemap_nopar.m
imagefiles_val={'/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_merge.tif'};
changefile=imagefiles_val{1};
data0=readGeotiff(changefile);

%step: 500 m ; buffer 200 m on each size.
buff=200; dx=500;dy=-dx;

rangx=data0.x(1):dx:data0.x(end);nx=length(rangx)-1;
rangy=data0.y(1):dy:data0.y(end);ny=length(rangy)-1;
ns=nx*ny;
for ix=1:nx
    for iy=1:ny
        ixy=iy+(ix-1)*ny;
        % rangi=[min(x)-buff max(x)+buff min(y)-buff max(y)+buff];
        rangi=[rangx(ix)-buff rangx(ix+1)+buff rangy(iy+1)-buff rangy(iy)+buff];
        idx=find(data0.x>=rangi(1)&data0.x<=rangi(2));
        idy=find(data0.y>=rangi(3)&data0.y<=rangi(4));
        data.x=data0.x(idx);data.y=data0.y(idy);
        data.z=data0.z(idy,idx);
        A1=uint8(rescale(data.z,0,255));
        ofile=['slumppeelsub/slumppeelsubi',num2str(ixy),'.jpg'];
        imwrite(A1,ofile);
        ids{ixy}.idy=idy;
        ids{ixy}.idx=idx;
    end
end
save ixy.mat data0 ids buff -v7.3

load ixy.mat

%merge masks together 
resr2m=2;
xout=min(data0.x(:))-buff:resr2m:max(data0.x(:))+buff;
yout=max(data0.y(:))+buff:-resr2m:min(data0.y(:))-buff;
mask=false(length(yout),length(xout)); %false
maskx=xout;masky=yout;
for ixy=1:ns %length(ids) %13%
% idx=ids{ixy}.idx;idy=ids{ixy}.idy;
% rang0=[min(data0.x(idx)) max(data0.x(idx)) min(data0.y(idy)) max(data0.y(idy))];
% ixy=iy+(ix-1)*ny;
ix=floor(ixy/ny)+1;iy=ixy-(ix-1)*ny;if(iy==0);ix=ix-1;iy=ny;end
rangi=[rangx(ix)-buff rangx(ix+1)+buff rangy(iy+1)-buff rangy(iy)+buff];
rang0=rangi;
x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
% hold on;plot(x0*1e-3,y0*1e-3,'r-','linewidth',3)
% ofile=['slumppeelsub/slumppeelsubi',num2str(ixy),'.jpg'];
% A1=imread(ofile);
% mask(idy,idx)=A1;
% figure;imagesc(data0.x(idx)*1e-3,data0.y(idy)*1e-3,A1);colorbar; colormap jet
idx=find(xout>=rangi(1)&xout<=rangi(2));
idy=find(yout>=rangi(3)&yout<=rangi(4));

maskiold=mask(idy,idx);
ofile=['slumppeelsubmask/slumppeelsubmaski',num2str(ixy),'.mat'];
if ~exist(ofile,'file')
    continue
end
A1=load(ofile);
A1=A1.arr;
% A1=imresize(A1.arr,size(maskiold));
mask(idy,idx)=maskiold|A1;
end %ixy
figure;imagesc(maskx*1e-3,masky*1e-3,mask);colorbar;

figure;imagesc(data0.x*1e-3,data0.y*1e-3,data0.z,'alphadata',mask);colorbar;
caxis([-3 3]);colormap jet

%Get polygon shape files.
Modfil=~(mask);  %Modfil (logical; 1 No change, 0 Change) 
Medgs1=false(size(mask));
% [Co]=mask2boundary(maskx,masky,(Modfil),Medgs1,'maskthres90_1732train2m_3bands_thres95.shp');%thres 90 vs 95. No big changes.
[Co]=mask2boundary(maskx,masky,(Modfil),Medgs1,'maskthres90_1552train_band1jpg.shp');%thres 90 vs 95. No big changes.


% Step 4: Assessement (estimating precision and recall).
data0=readGeotiff('/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_jump.tif');

smlarea=1e3*2; %1000 m^2
resx=mean(data0.x(2:end)-data0.x(1:end-1));resy=mean(data0.y(2:end)-data0.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);

%Get the examed slump area (as ground truth).
Mf4=data0.z<=-1; %abs(jump)>=2; %jump significant;Only negative
Mf4clean = bwareaopen(Mf4, smlarea/resr/resr); %remove small clusters
data0x=data0.x;data0y=data0.y;
projstr='polar stereo north';
save t1.mat data0x data0y Mf4clean projstr -v7.3
writeGeotiff('maskslumpraw1.tif',data0x,data0y,uint8(Mf4clean),1,0,projstr)

% Get the mask drawn manually that contains slumps.
S1=shaperead('/Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/peelmanual.shp');% figure;mapshow(S1)
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
%Get polygon shape files.
Modfil=~(Mf5clean);  %Modfil (logical; 1 No change, 0 Change) 
Medgs1=false(size(Mf5clean));
[Co]=mask2boundary(data0x,data0y,(Modfil),Medgs1,'maskgroundtruth.shp');
save maskgroundtruth.mat data0x data0y Mf5clean -v7.3

% Get precision and recall
S2=shaperead('/Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/maskgroundtruth.shp');% figure;mapshow(S1)
Sml=shaperead('maskthres90_1552train_band1jpg.shp');
Sml=shaperead('maskthres90_48ktrain2m_band1jpg.shp');
% Sml=shaperead('maskthres90_1732train2m_band1_uint8jpg.shp');
% Sml=shaperead('maskthres90_1732train2m.shp');
[precision,recall]=poly2precision(Sml,S2);

%% Manually crop sub images.
changefile=imagefiles{1};
data0=readGeotiff(changefile);
M=data0.z<-50|data0.z>50;
M=data0.z<-10|data0.z>10;
M=data0.z<-5|data0.z>5;
figure;imagesc(data0.x*1e-3,data0.y*1e-3,data0.z,'alphadata',M);colorbar;colormap jet;
data0.z(M)=0;
A1=uint8(rescale(data0.z,0,255));
imwrite(A1,'slumppeel50km_tr50.jpg'); %truncate all values > 50 or <-50
imwrite(A1,'slumppeel50km_tr10.jpg'); %truncate all values > 10 or <-10
imwrite(A1,'slumppeel50km_tr5.jpg');

mp3=getline;
x=mp3(:,1)*1e3;y=mp3(:,2)*1e3;
hold on;plot3(x*1e-3,y*1e-3,9e3*ones(size(mp3(:,1))),'k-')
buff=100; %buffer the polygon by 100m
rangi=[min(x)-buff max(x)+buff min(y)-buff max(y)+buff];
idx=find(data0.x>=rangi(1)&data0.x<=rangi(2));
idy=find(data0.y>=rangi(3)&data0.y<=rangi(4));
data.x=data0.x(idx);data.y=data0.y(idy);
data.z=data0.z(idy,idx);
A1=uint8(rescale(data.z,0,255));
imwrite(A1,'slumppeelslump13big.jpg');
imwrite(A1,'slumppeelslump13big2.jpg');
imwrite(A1,'slumppeelslump13big3.jpg');
imwrite(A1,'slumppeelslump13new.jpg');
imwrite(A1,'slumppeelslump13new1km.jpg');


data=readGeotiff('/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_04_04_jump.tif');
data=readGeotiff('/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_03_05_jump.tif');
A1=uint8(rescale(data.z,0,255));
imwrite(A1,'41_16_1_1_03_05_jump.jpg');
imwrite(A1,'slumppeelslump13tile.jpg');
imwrite(A1,'slumppeel50kmzone1s1.jpg');
imwrite(A1,'slumppeel50kmzone1.jpg'); 


