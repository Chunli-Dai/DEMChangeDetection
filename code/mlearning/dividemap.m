constant
codedir='/home/dai.56/arcticdemapp/landslide/code1/';  %Directory of codes.
flagorthotif=1; % save mosacked ortho image in Geotiff

regionstr='eureka';

if flagorthotif==1
        %save orthoimage in geotiff 
        odirt1=['slumps_ortho/'];
        if ~exist(odirt1,'')
            mkdir(odirt1)
        end
end
%mkdir slumpsub_3bandsuint8jpg
%mkdir slumpsub_band1changeuint8jpg
%mkdir slumpsub_band2demuint8jpg
%mkdir slumpsub_band3orthouint8jpg


addpath(genpath([codedir]));

load mat0.mat

% Divide target image to ~1 km size;
imagefiles_val={'/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_merge.tif'};
imagefiles_val={'41_16_1_1_merge.tif'};
imagefiles_val={'32_33_2_2_jump.tif'};
changefile=imagefiles_val{1};
data0=readGeotiff(changefile);

%step: 500 m ; buffer 200 m on each size. for Peel Plateau
buff=200; dx=500;dy=-dx;
% try 800 pixels (2m ): 1600 m ; For Eureka
%buff=200; dx=1200;dy=-dx;
%buff=100; dx=300;dy=-dx; %-> Chandi uses 250 pixels by 250 pixels;

rangx=data0.x(1):dx:data0.x(end);nx=length(rangx)-1;
rangy=data0.y(1):dy:data0.y(end);ny=length(rangy)-1;
ns=nx*ny;
ids=[];
poolsize=20;
poolobj=parpool(20);
parfor ix=1:nx
    for iy=1:ny
        ixy=iy+(ix-1)*ny;
	fprintf(['\n Working on ',num2str(ixy),'/',num2str(ns),' tiles.\n'])
        % rangi=[min(x)-buff max(x)+buff min(y)-buff max(y)+buff];
        rangi=[rangx(ix)-buff rangx(ix+1)+buff rangy(iy+1)-buff rangy(iy)+buff];

	[datao]=box2mosaic(rangi,3,f,fdir); %test hi

        %save orthoimage in geotiff 
        if flagorthotif==1
        filenamet1i=[odirt1,'/','slump',num2str(ixy),'.tif'];
        projstr='polar stereo north';
        writeGeotiff(filenamet1i,datao.x,datao.y,uint16(datao.z(:,:,3)),2,0,projstr)
        end

	%A1=uint16(rescale(datao.z(:,:,1),0,65535)); %Detectron2 only works with uint8!
	A1=uint8(rescale(datao.z(:,:,1),0,255));
	A1(:,:,2)=uint8(rescale(datao.z(:,:,2),0,255)); %change DEM ORTHO
	A1(:,:,3)=uint8(rescale(datao.z(:,:,3),0,255));

        %ofile=['slumppeelsub/slumppeelsubi',num2str(ixy),'.tif'];
        %ofile=['slumppeelsub_3bandsuint8/slumppeelsubi',num2str(ixy),'.tif']; %TIFF format
        %ofile=['slumppeelsub_3bandsuint8jpg/slumppeelsubi',num2str(ixy),'.jpg'];
        ofile=['slumpsub_3bandsuint8jpg/slumppeelsubi',num2str(ixy),'.jpg'];

        imwrite(A1,ofile); %TIFF format
        %imwrite(A1,ofile,'Mode','lossless'); %jpg; avoid compression ; cause error in Detectron2: 'NoneType' object has no attribute 'shape' 

        %ids{ixy}.idy=idy;
        %ids{ixy}.idx=idx;
	
	A0=A1; 
	A1=A0(:,:,1);
        %ofile=['slumppeelsub_band1change/slumppeelsubi',num2str(ixy),'.tif'];
        %ofile=['slumppeelsub_band1changeuint8/slumppeelsubi',num2str(ixy),'.tif'];
        %ofile=['slumppeelsub_band1changeuint8jpg/slumppeelsubi',num2str(ixy),'.jpg']; 
        ofile=['slumpsub_band1changeuint8jpg/slumppeelsubi',num2str(ixy),'.jpg']; 
        imwrite(A1,ofile);

	A1=A0(:,:,2);
%       ofile=['slumppeelsub_band2dem/slumppeelsubi',num2str(ixy),'.tif'];
        %ofile=['slumppeelsub_band2demuint8/slumppeelsubi',num2str(ixy),'.tif'];
        %ofile=['slumppeelsub_band2demuint8jpg/slumppeelsubi',num2str(ixy),'.jpg'];
        ofile=['slumpsub_band2demuint8jpg/slumppeelsubi',num2str(ixy),'.jpg'];
        imwrite(A1,ofile);

	A1=A0(:,:,3);
        %ofile=['slumppeelsub_band3ortho/slumppeelsubi',num2str(ixy),'.tif'];
        %ofile=['slumppeelsub_band3orthouint8/slumppeelsubi',num2str(ixy),'.tif'];
        %ofile=['slumppeelsub_band3orthouint8jpg/slumppeelsubi',num2str(ixy),'.jpg'];
        ofile=['slumpsub_band3orthouint8jpg/slumppeelsubi',num2str(ixy),'.jpg'];
        imwrite(A1,ofile);
    end
end
delete(poolobj)

%save ixy.mat data0 ids buff -v7.3
save ixy_eureka.mat data0 ids buff dx rangx rangy ns -v7.3
