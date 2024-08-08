
function [filenames,filesize,XYb,absfilenames,countaug]=applyaug(filenames,filesize,XYb,absfilenames, countaug,dirfilename,A1,XYbi)

%Data augmentation: blur crop scaling(resize) rotation flip
%Increase sample by ~30 (10+10+10) times 
% see https://www.mathworks.com/help/deeplearning/ug/image-augmentation-using-image-processing-toolbox.html

    filenameorg=dirfilename; %e.g. ./train/slump1.jpg
    imOriginal=A1;
    %get the mask of original polygon XYbi
    [ny,nx,nz]=size(A1);
    maskorg=false(size(ny,nx));
    for k=1:length(XYbi) %cell
    maskorg_k=poly2mask(XYbi{k}(:,1),XYbi{k}(:,2), ny,nx); 
    maskorg=maskorg|maskorg_k;
    end
    
    %initializde 
    strg={[],[],[],[],[]};
    
    multi=30; % 30 times of the original
    multi=100; % 30 times of the original
    %blur
    for iaug=[] %1:multi %10 %100
        countaug=countaug+1;
                
        sigma = 1+5*rand; 
        imAugmented = imgaussfilt(imOriginal,sigma); 
        maskaug=maskorg;XYbi_aug=XYbi;
%         B = bwboundaries(maskaug); %assume one single polygon
%         k=1;xid=B{k}(:,2); yid=B{k}(:,1); XYbi_aug=[xid(:),yid(:)];
        
        strg{5}=['blur',num2str(iaug)];
        strcom=[strg{5},'.jpg'];
        dirfilename_aug=strrep(dirfilename,'.jpg',strcom);
        imwrite(imAugmented,dirfilename_aug);  
        
        %figure; hold on;imagesc(imAugmented);hold on; hold on;plot(XYbi_aug(:,1),XYbi_aug(:,2),'k.-'); title(dirfilename_aug);axis equal;pause;
        
        [filenames,filesize,XYb,absfilenames]=saveattributes(filenames,filesize,XYb,absfilenames,countaug,dirfilename_aug,XYbi_aug);
    end
    
    %crop
    for iaug=[] %1:multi  %[] %1:multi %10 %100
        countaug=countaug+1;
        
        a=0.7;b=0.99;p=a + (b-a).*rand(1);
        targetSize=round(p*[ny,nx]);
        win = randomCropWindow2d(size(imOriginal),targetSize);
        imAugmented = imcrop(imOriginal,win); 
        maskaug = imcrop(maskorg,win); 
%       B = bwboundaries(maskaug); %assume one single polygon
%       k=1;xid=B{k}(:,2); yid=B{k}(:,1); XYbi_aug=[xid(:),yid(:)];
        
	Mf4clean = bwareaopen(maskaug, 25); %remove small clusters<25*4 m^2
        B = bwboundaries(Mf4clean); %assume one single polygon
        XYbi_aug=[];idd=[];
        for k=1:length(B)
        xid=B{k}(:,2); yid=B{k}(:,1);
	%avoid bug ValueError: Failed to use mask_format=='polygon' from the given annotations!
	if (max(xid)-min(xid))==0 || (max(yid)-min(yid))==0;idd=[idd(:);k]; end  %point or line
	XYbi_aug{k}=[xid(:),yid(:)];
        end
	XYbi_aug(idd)=[];
        
        strg{5}=['crop',num2str(iaug)];
        strcom=[strg{5},'.jpg'];
        dirfilename_aug=strrep(dirfilename,'.jpg',strcom);
        imwrite(imAugmented,dirfilename_aug);  
        
        %figure; hold on;imagesc(imAugmented);hold on; hold on;plot(XYbi_aug(:,1),XYbi_aug(:,2),'k.-'); title(dirfilename_aug);axis equal;pause;
        
        [filenames,filesize,XYb,absfilenames]=saveattributes(filenames,filesize,XYb,absfilenames,countaug,dirfilename_aug,XYbi_aug);
    end
    
    %scaling
    for iaug=1:2; %[] %1:multi % []% 1:multi %[] %1:multi %10 %100
        countaug=countaug+1;
                
        tform = randomAffine2d('Scale',[0.75,1.25]);
        outputView = affineOutputView(size(imOriginal),tform);
        imAugmented = imwarp(imOriginal,tform,'OutputView',outputView);
        maskaug = imwarp(maskorg,tform,'OutputView',outputView);

%       B = bwboundaries(maskaug); %assume one single polygon
%       k=1;xid=B{k}(:,2); yid=B{k}(:,1); XYbi_aug=[xid(:),yid(:)];
        
	Mf4clean = bwareaopen(maskaug, 25); %remove small clusters<25*4 m^2
        B = bwboundaries(Mf4clean); %assume one single polygon
        XYbi_aug=[];idd=[];
        for k=1:length(B)
        xid=B{k}(:,2); yid=B{k}(:,1);
	%avoid bug ValueError: Failed to use mask_format=='polygon' from the given annotations!
	if (max(xid)-min(xid))==0 || (max(yid)-min(yid))==0;idd=[idd(:);k]; end  %point or line
	XYbi_aug{k}=[xid(:),yid(:)];
        end
	XYbi_aug(idd)=[];
        
        strg{5}=['scale',num2str(iaug)];
        strcom=[strg{5},'.jpg'];
        dirfilename_aug=strrep(dirfilename,'.jpg',strcom);
        imwrite(imAugmented,dirfilename_aug);  
        
        %figure; hold on;imagesc(imAugmented);hold on; hold on;plot(XYbi_aug(:,1),XYbi_aug(:,2),'k.-'); title(dirfilename_aug);axis equal;pause;
        
        [filenames,filesize,XYb,absfilenames]=saveattributes(filenames,filesize,XYb,absfilenames,countaug,dirfilename_aug,XYbi_aug);
    end
    
    %rotation
    for iaug=1:18; %1:multi %[] %1:10 %100
        countaug=countaug+1;
        
        tform = randomAffine2d('Rotation',[0 359]); 
        outputView = affineOutputView(size(imOriginal),tform);
        imAugmented = imwarp(imOriginal,tform,'OutputView',outputView); 
        maskaug = imwarp(maskorg,tform,'OutputView',outputView);
%       B = bwboundaries(maskaug); %assume one single polygon
%       k=1;xid=B{k}(:,2); yid=B{k}(:,1); XYbi_aug=[xid(:),yid(:)];
	Mf4clean = bwareaopen(maskaug, 25); %remove small clusters<25*4 m^2
        B = bwboundaries(Mf4clean); %assume one single polygon
        XYbi_aug=[];idd=[];
        for k=1:length(B)
        xid=B{k}(:,2); yid=B{k}(:,1);
	%avoid bug ValueError: Failed to use mask_format=='polygon' from the given annotations!
	%if (max(xid)-min(xid))==0 || (max(yid)-min(yid))==0;idd=[idd(:);k]; continue;end  %point or line
	if (max(xid)-min(xid))==0 || (max(yid)-min(yid))==0;idd=[idd(:);k]; end  %point or line
	XYbi_aug{k}=[xid(:),yid(:)];
        end
	XYbi_aug(idd)=[];
        
        strg{5}=['rotation',num2str(iaug)];
        strcom=[strg{5},'.jpg'];
        dirfilename_aug=strrep(dirfilename,'.jpg',strcom);
        imwrite(imAugmented,dirfilename_aug);  
        
        %figure; hold on;imagesc(imAugmented);hold on; hold on;plot(XYbi_aug(:,1),XYbi_aug(:,2),'k.-'); title(dirfilename_aug);axis equal;pause;
        
        [filenames,filesize,XYb,absfilenames]=saveattributes(filenames,filesize,XYb,absfilenames,countaug,dirfilename_aug,XYbi_aug);
    end
        
    %flip random
    for iaug=1:2; %[] %1:2 %100
        countaug=countaug+1;
                
        tform = randomAffine2d('XReflection',true,'YReflection',true);
        outputView = affineOutputView(size(imOriginal),tform);
        imAugmented = imwarp(imOriginal,tform,'OutputView',outputView);
        maskaug = imwarp(maskorg,tform,'OutputView',outputView);

%        B = bwboundaries(maskaug); %assume one single polygon
%        k=1;xid=B{k}(:,2); yid=B{k}(:,1); XYbi_aug=[xid(:),yid(:)];
	Mf4clean = bwareaopen(maskaug, 25); %remove small clusters<25*4 m^2
        B = bwboundaries(Mf4clean); %assume one single polygon
        XYbi_aug=[];idd=[];
        for k=1:length(B)
        %xid=B{k}(:,2); yid=B{k}(:,1); XYbi_aug{k}=[xid(:),yid(:)];
        xid=B{k}(:,2); yid=B{k}(:,1);
	%avoid bug ValueError: Failed to use mask_format=='polygon' from the given annotations!
	if (max(xid)-min(xid))==0 || (max(yid)-min(yid))==0;idd=[idd(:);k]; end  %point or line
	XYbi_aug{k}=[xid(:),yid(:)];
        end
	XYbi_aug(idd)=[];
        
        strg{5}=['flip',num2str(iaug)];
        strcom=[strg{5},'.jpg'];
        dirfilename_aug=strrep(dirfilename,'.jpg',strcom);
        imwrite(imAugmented,dirfilename_aug);  
        
        %figure; hold on;imagesc(imAugmented);hold on; hold on;plot(XYbi_aug(:,1),XYbi_aug(:,2),'k.-'); title(dirfilename_aug);axis equal;pause;
        
        [filenames,filesize,XYb,absfilenames]=saveattributes(filenames,filesize,XYb,absfilenames,countaug,dirfilename_aug,XYbi_aug);
    end
    

    %translation
    for iaug=1:2; %1:multi %[] %1:10 %100

        countaug=countaug+1;

        tform = randomAffine2d('XTranslation',[-100 100],'YTranslation',[-100 100]);
        outputView = affineOutputView(size(imOriginal),tform);
        imAugmented = imwarp(imOriginal,tform,'OutputView',outputView);
        maskaug = imwarp(maskorg,tform,'OutputView',outputView);
	Mf4clean = bwareaopen(maskaug, 25); %remove small clusters<25*4 m^2
        B = bwboundaries(Mf4clean); %assume one single polygon
%         k=1;xid=B{k}(:,2); yid=B{k}(:,1); XYbi_aug=[xid(:),yid(:)];
        XYbi_aug=[];idd=[];
        for k=1:length(B)
        xid=B{k}(:,2); yid=B{k}(:,1);
	%avoid bug ValueError: Failed to use mask_format=='polygon' from the given annotations!
	if (max(xid)-min(xid))==0 || (max(yid)-min(yid))==0;idd=[idd(:);k]; end  %point or line
	XYbi_aug{k}=[xid(:),yid(:)];
        end
	XYbi_aug(idd)=[];

        strg{5}=['transl',num2str(iaug)];
        strcom=[strg{5},'.jpg'];
        dirfilename_aug=strrep(dirfilename,'.jpg',strcom);
        imwrite(imAugmented,dirfilename_aug);

        [filenames,filesize,XYb,absfilenames]=saveattributes(filenames,filesize,XYb,absfilenames,countaug,dirfilename_aug,XYbi_aug);


%    tform = randomAffine2d('XTranslation',[-5 5],'YTranslation',[-5 5]);
%    outputView = affineOutputView(size(imOriginal),tform);
%    imAugmented = imwarp(imOriginal,tform,'OutputView',outputView);
%    imshow(imAugmented)
    end
    
%         strcom=[strg{1},strg{2},strg{3},strg{4},strg{5},'.jpg']; %combined string

%Moved to saveattributes.m
%         [~,fname,ext]=fileparts(dirfilename_aug);
%         filenamei=[fname,ext];
%         
%         % get file attributes
%         str1=imfinfo(dirfilename_aug);
%         filesizei=str1.FileSize;
% 
%         filenames{countaug}=filenamei;
%         filesize(countaug)=filesizei;
%         XYb{countaug}=XYbi_aug;
%         absfilenames{countaug}=dirfilename_aug;


end
