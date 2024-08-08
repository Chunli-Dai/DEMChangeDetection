constant

%New versoin May 2023: Use the given candidates from Change2poly6classespar.m to get test images.
%new version Jan, 2022: use constraint 1 to get candidate clusters, and then get images for all clusters.

%codedir='/home/dai.56/arcticdemapp/landslide/code1/';  %Directory of codes.
flagorthotif=0; % save mosacked ortho image in Geotiff
flagdemtif=0; % save mosacked dem in Geotiff
flagrescale=0;  %1, fixed range [e.g. -20 10] to keep the information of positive/negative.
                %0, re-adjust the range for each image.
flagband=3; %3; %1 writing 1 band data for training images; 3 writing 3 band data; 2: elevation change and DEM curvature, two bands.

%regionstr='eureka';
%regionstr='russia';
%regionstr='peel';
regionstr='region09';
shapefile='../region09/v3con/change2polyg_09.shp'; %candidates
mergedfile=[regiondir,'/region09_jump.tif'];

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

%get fdir f
load mat0.mat

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
changefile=imagefiles_val{1};

%to get the grids
data0=readGeotiff(changefile);

data0.z=false(size(data0.z));
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

[X,Y]=meshgrid(data0.x,data0.y);

%ids={};
ids=cell(ns,1);
%poolsize=20;
%poolobj=parpool(10);
%parfor ix=1:ns
for ix=1:ns

    fprintf(['\n Working on ',num2str(ix),'/',num2str(ns),' candidate cluster.\n'])

    ixy=ix; i=ix;

    % see prepareDetectron2.m changemap2polygon.m
    %get the box in polar coordinates
    ll=S1(ix).BoundingBox(1,:);ru=S1(ix).BoundingBox(2,:);
    [llx,lly]=polarstereo_fwd(ll(2),ll(1),[],[],70,-45); %lat lon
    [rux,ruy]=polarstereo_fwd(ru(2),ru(1),[],[],70,-45); %lat lon
    [rlx,rly]=polarstereo_fwd(ll(2),ru(1),[],[],70,-45); %right low
    [lux,luy]=polarstereo_fwd(ru(2),ll(1),[],[],70,-45); %left up
    x=[llx,rlx,rux,lux];y=[lly,rly,ruy,luy];
    rang2=[min(x)-buff max(x)+buff min(y)-buff max(y)+buff];

    Mx=data0.x>=rang2(1)&data0.x<=rang2(2);My=data0.y>=rang2(3)&data0.y<=rang2(4);
    
    %cropped clusters
    Mpc=[];
    Mpc.x=data0.x(Mx);Mpc.y=data0.y(My);
    Mpcpre.z=data0.z(My,Mx);

    %store the candidate mask to data0.z;
    [s1x,s1y]=polarstereo_fwd(S1(ix).Y,S1(ix).X,[],[],70,-45);
    Xb=round((s1x-Mpc.x(1))/resrx)+1; 
    Yb=round((s1y-Mpc.y(1))/(resry))+1; %descending        
    XYbi=[Xb(:),Yb(:)]; %n by 2
    id1=find(any(isnan(XYbi),2));
    XYbi(id1,:)=[];
    [nyi,nxi,nzi]=size(Mpcpre.z);
    Mpc.z=poly2mask(XYbi(:,1),XYbi(:,2), nyi,nxi); 
    data0.z(My,Mx)=Mpc.z|Mpcpre.z;

    %get the test images

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
        if flagdemtif==1
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

end
%delete(poolobj)

%save ixy.mat data0 ids buff -v7.3
%save ixy_eureka.mat data0 ids buff dx rangx rangy ns -v7.3
ofile=['ixy_',regionstr,'_cand.mat'];
save(ofile,'data0', 'ids', 'ns', '-v7.3')
