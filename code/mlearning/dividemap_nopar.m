constant

%new version Jan, 2022: use constraint 1 to get candidate clusters, and then get images for all clusters.

%codedir='/home/dai.56/arcticdemapp/landslide/code1/';  %Directory of codes.
flagorthotif=0; % save mosacked ortho image in Geotiff
flagdemtif=0; % save mosacked dem in Geotiff
flagrescale=0;  %1, fixed range [e.g. -20 10] to keep the information of positive/negative.
                %0, re-adjust the range for each image.
flagband=2; %3; %1 writing 1 band data for training images; 3 writing 3 band data; 2: elevation change and DEM curvature, two bands.

regionstr='eureka';
%regionstr='russia';
regionstr='peel';

if flagorthotif==1
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

load mat0.mat

% Divide target image to ~1 km size;
imagefiles_val={'/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_merge.tif'};
imagefiles_val={'41_16_1_1_merge.tif'}; %peel
imagefiles_val={'41_16_1_1_jump.tif'}; %peel
%imagefiles_val={'32_33_2_2_jump.tif'}; %eureka
%imagefiles_val={'49_60_2_2_jump.tif'}; %nina
changefile=imagefiles_val{1};
data0=readGeotiff(changefile);

%step: 500 m ; buffer 200 m on each size. for Peel Plateau
buff=200; dx=500;dy=-dx;
% try 800 pixels (2m ): 1600 m ; For Eureka
buff=100; dx=300;dy=-dx; %-> Chandi uses 250 pixels by 250 pixels;
buff=50; dx=200;
dy=-dx;

buff=100; %buffer the polygon by 100m % consistent with training images prepareDetectron2.m

%odir1=['slumpsub_band1changeuint8jpg_',regionstr,'_dx',num2str(dx),'/']; 
odir1=['slumpsub_band3_',regionstr,'_dx',num2str(dx),'/']; 
odir1='slumpsub_band1_peel_dx200_fixr_curv/';
odir1='slumpsub_band1_peel_dx200_wofixr_curv/';
odir1='slumpsub_band1_peel_fixr_curv_cand/';
odir1='slumpsub_band1_peel_wofixr_curv_cand/';
odir1='slumpsub_peel_timeband3/';
if ~exist(odir1,'')
   mkdir(odir1)
end

%get candidates of clusters based on constraint 1: 
smlarea=1e3*2; %1000 m^2
resx=mean(data0.x(2:end)-data0.x(1:end-1));resy=mean(data0.y(2:end)-data0.y(1:end-1));
resr=mean([abs(resx),abs(resy)]) 
%Get the examed slump area (as ground truth).
Mf4=data0.z<=-1; %Only negative
Mf4clean = bwareaopen(Mf4, round(smlarea/resr/resr)); %remove small clusters
clear Mf4

CC = bwconncomp(Mf4clean);
ns=CC.NumObjects;
fprintf(['There are a total of ',num2str(CC.NumObjects),' candidate clusters.\n'])

[X,Y]=meshgrid(data0.x,data0.y);

%ids={};
ids=cell(ns,1);
%poolsize=20;
%poolobj=parpool(10);
%parfor ix=1:ns
for ix=1:ns
	fprintf(['\n Working on ',num2str(ix),'/',num2str(ns),' candidate cluster.\n'])

	ki=ix; ixy=ix;
    %get the individual cluster.
	BW3=[];Mx=[];My=[];
    BW3=false(size(Mf4clean));
    BW3(CC.PixelIdxList{ki})=Mf4clean(CC.PixelIdxList{ki});

    xmin=min(X(BW3));xmax=max(X(BW3));ymin=min(Y(BW3));ymax=max(Y(BW3));
    rang2=[xmin-buff xmax+buff ymin-buff ymax+buff];
    Mx=data0.x>=rang2(1)&data0.x<=rang2(2);My=data0.y>=rang2(3)&data0.y<=rang2(4);
    
    %cropped clusters
	Mpc=[];
    Mpc.x=data0.x(Mx);Mpc.y=data0.y(My);
    Mpc.z=data0.z(My,Mx);

    rang2=[min(Mpc.x) max(Mpc.x) min(Mpc.y) max(Mpc.y)]; %make sure the boundary is consistent with ids.

%idxs=[]; idxe=[]; idys=[]; idye=[];
	idxs=find(Mx,1,'first');idxe=find(Mx,1,'last');
	idys=find(My,1,'first');idye=find(My,1,'last');

        %ids{ix}.My=My;
        %ids{ix}.Mx=Mx;
        ids{ix}=[idxs idxe idys idye];

        rangi=rang2;

%fprintf('hi1')
 	[datao]=box2mosaic(rangi,flagband,f,fdir); %test hi
%[datao]=box2mosaic(rangi,1,f,fdir); %test hi
%	fprintf('hi2')
%	save testhi1.mat -v7.3

        %save orthoimage in geotiff 
        if flagorthotif==1
        filenamet1i=[odirt1,'/','slump',num2str(ixy),'.tif'];
        projstr='polar stereo north';
        writeGeotiff(filenamet1i,datao.x,datao.y,uint16(datao.z(:,:,3)),2,0,projstr)
        end
        if 1 %flagdemtif==1
        filenamet1i=[odirt2,'/','slump',num2str(ixy),'.tif'];
        projstr='polar stereo north';
        %writeGeotiff(filenamet1i,datao.x,datao.y,int16(datao.z(:,:,2)),2,0,projstr)
        writeGeotiff(filenamet1i,datao.x,datao.y,int16(datao.z),2,0,projstr)
        end

	%A1=uint16(rescale(datao.z(:,:,1),0,65535)); %Detectron2 only works with uint8!
	A1=[];T2=[];
	if flagrescale==0
	    A1=uint8(rescale(datao.z(:,:,1),0,255));
	elseif flagrescale==1
	    inmin=-20;inmax=10; %truncate data to this range [inmin inmax], then rescale.
	    A1=uint8(rescale(datao.z(:,:,1),0,255,'InputMin',inmin,'InputMax',inmax));
	end

	if flagband==2 
          %get curvature
          totalc=datao.z(:,:,2);
          A1(:,:,2)=uint8((totalc(:,:))); %DEM;
          A1(:,:,3)= A1(:,:,1); %  %Data with 2 components not supported for JPEG files. So add 3rd band.
	elseif flagband ==4
	  T2=datao.z(:,:,2);
          A1(:,:,2)=uint8((T2(:,:))); %
          A1(:,:,3)= A1(:,:,1); %
	elseif flagband ==3
	  A1(:,:,2)=uint8((datao.z(:,:,2))); % %curvature
	  A1(:,:,3)=uint8((datao.z(:,:,3))); % %time
	end

        %ofile=['slumpsub_band1changeuint8jpg/slumppeelsubi',num2str(ixy),'.jpg']; 
        ofile=[odir1,'/slumppeelsubi',num2str(ixy),'.jpg']; 
        imwrite(A1,ofile);

%	if 0
%	A1(:,:,2)=uint8(rescale(datao.z(:,:,2),0,255)); %change DEM ORTHO
%	A1(:,:,3)=uint8(rescale(datao.z(:,:,3),0,255));
%
%        ofile=['slumpsub_3bandsuint8jpg/slumppeelsubi',num2str(ixy),'.jpg'];
%
%        imwrite(A1,ofile); %TIFF format
%        %imwrite(A1,ofile,'Mode','lossless'); %jpg; avoid compression ; cause error in Detectron2: 'NoneType' object has no attribute 'shape' 
%
%        %ids{ixy}.idx=idx;
%	
%	A0=A1; 
%	A1=A0(:,:,1);
%        ofile=['slumpsub_band1changeuint8jpg/slumppeelsubi',num2str(ixy),'.jpg']; 
%        imwrite(A1,ofile);
%
%	A1=A0(:,:,2);
%        ofile=['slumpsub_band2demuint8jpg/slumppeelsubi',num2str(ixy),'.jpg'];
%        imwrite(A1,ofile);
%
%	A1=A0(:,:,3);
%        ofile=['slumpsub_band3orthouint8jpg/slumppeelsubi',num2str(ixy),'.jpg'];
%        imwrite(A1,ofile);
%	end  %if 0
end
delete(poolobj)

%save ixy.mat data0 ids buff -v7.3
%save ixy_eureka.mat data0 ids buff dx rangx rangy ns -v7.3
ofile=['ixy_',regionstr,'_cand.mat'];
save(ofile,'data0', 'ids', 'ns', '-v7.3')
