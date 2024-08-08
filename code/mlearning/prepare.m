function [trainingImages1,trainingLabels1,slumpbox]=prepare(shapefiles,imagefiles)
% % Prepare training images for RCNN
% Given shapefiles, and imagefiles to get training image matrices.
% shapefiles: shapefile of known slump outlines.
% imagefiles: elevation change files list
% output: slumpbox: the boundary of slumps

% trainingLabels1=categorical({'slump', 'non slump'}');

flagplot=1; %1 plot;0 not plot
flagaug=1; %if use data augment;

% S1=shaperead('/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Ashley/shapefiles/Slumps_2017_18_SCAR_FINAL.shp');
S1=shaperead(shapefiles);% figure;mapshow(S1)
%/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_06_21_jump.tif

% changefile='/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_merge.tif';
changefile=imagefiles{1};
data0=readGeotiff(changefile);
resr=mean(diff(data0.x));
Maskslump=false(size(data0.z));

%Get slumps
n=length(S1);
count=0; %clear trainingImages1 trainingLabels1 collecti
buff=100; %buffer the polygon by 100m
for i=1:n
    if contains(S1(i).zone,'debris');continue;end %only detect scars
    ll=S1(i).BoundingBox(1,:);ru=S1(i).BoundingBox(2,:);
    [llx,lly]=polarstereo_fwd(ll(2),ll(1),[],[],70,-45); %lat lon
    [rux,ruy]=polarstereo_fwd(ru(2),ru(1),[],[],70,-45); %lat lon
    x=[llx,rux];y=[lly,ruy];
    rangi=[min(x)-buff max(x)+buff min(y)-buff max(y)+buff];
%     data=readGeotiff(changefile,'map_subset',rangi);
    idx=find(data0.x>=rangi(1)&data0.x<=rangi(2));
    idy=find(data0.y>=rangi(3)&data0.y<=rangi(4));
    data.x=data0.x(idx);data.y=data0.y(idy);
    data.z=data0.z(idy,idx);
    [mi,ni]=size(data.z); minm=min([mi,ni]);
    if (minm<3);continue
    else
        Maskslump(idy,idx)=1; %store all the slump masks;
        
        count=count+1;
        imagei=single(imresize(data.z,[32,32]));
        trainingImages1a(1:32,1:32,1,count)=imagei;
        trainingLabels1a(count)=categorical({'slump'});
        collecti(count)=i;
        
        pixx=round((rangi(1)-data0.x(1))/resr)+1; 
        pixy=round((rangi(4)-data0.y(1))/(-resr))+1; %descending
        pixdx=round((rangi(2)-rangi(1))/resr); pixdy=round((rangi(4)-rangi(3))/resr);
        slump.slump(count)={[pixx pixy pixdx pixdy]}; %pixel cooridnate: x y dx dy

        if 0 %for plotting
            [s1x,s1y]=polarstereo_fwd(S1(i).Y,S1(i).X,[],[],70,-45);
            figure;imagesc(data.x*1e-3,data.y*1e-3,data.z);colormap jet; colorbar
            caxis([-1 1])
            hold on;plot3(rux*1e-3,ruy*1e-3,9e9,'k>','linewidth',3)
            axis(rangi*1e-3)
            hold on;plot(llx*1e-3,lly*1e-3,'k>','linewidth',3)
            hold on;plot(s1x*1e-3,s1y*1e-3,'k.-','linewidth',3)
            
%             pixx=round((s1x-data0.x(1))/resr)+1; 
%             pixy=round((s1y-data0.y(1))/(-resr))+1; %descending
%             hold on;plot(pixx,pixy,'k.-','linewidth',3)

            figure;imagesc(imagei)
            figure;mapshow(S1(i))
            
        end
    end
end
nslump=length(slump.slump);
% slumpbox=slump;

for i=1:nslump
slump.imageFilename(i,1)={changefile};
end

T = table(slump.imageFilename,slump.slump','VariableNames',{'imageFilename','slump'});
slumpbox=T;

%% flagaug=1; %if use data augment;
% https://www.mathworks.com/help/deeplearning/ug/image-augmentation-using-image-processing-toolbox.html
% https://benanne.github.io/2015/03/17/plankton.html
%given trainingImages1 trainingLabels1
if flagaug==1
%     figure
    count=0;
    for j=1:nslump
        imOriginal=trainingImages1a(1:32,1:32,1,j);
        
    for i=1:20 %100
    tform = randomAffine2d('Rotation',[0 360]); 
    outputView = affineOutputView(size(imOriginal),tform);
    imAugmented = imwarp(imOriginal,tform,'OutputView',outputView);  
%     hold on;imagesc(imAugmented);title(num2str(i));pause;

        count=count+1;
        trainingImages1(1:32,1:32,1,count)=imAugmented;
        trainingLabels1(count)=trainingLabels1a(j);
        
    end
    
    end

    %translation
    if 0
    tform = randomAffine2d('XTranslation',[-5 5],'YTranslation',[-5 5]);
    outputView = affineOutputView(size(imOriginal),tform);
    imAugmented = imwarp(imOriginal,tform,'OutputView',outputView);
    imshow(imAugmented)
    end
else
        trainingImages1=trainingImages1a;
        trainingLabels1=trainingLabels1a;
end

[~,~,~,countslump]=size(trainingImages1);

%% get non slumps
[ny,nx]=size(data0.z);

if flagplot==1
 figure;imagesc(data0.z);colorbar;caxis([-3 3])
end
% count=0;
count2=0;
maxc=4000;%4000; %40;
flag=0;%
for i=1:32:nx-32
for j=1:32:ny-32
%         idy=1:32;idx=(1:32)+(j-1)*32; %j=1:40;
        idy=j+(1:32)-1;idx=i+(1:32)-1;
        imagei=single((data0.z(idy,idx)));
        
        %check it's not a slump
        M1=false(size(data0.z));
        M1(idy,idx)=1;
        overlap=sum(sum(M1&Maskslump));
        
        if overlap==0 %no slump
            count=count+1;
            count2=count2+1;
            
            trainingImages1(1:32,1:32,1,count)=imagei;
            trainingLabels1(count)=categorical({'non slump'});
%           collecti(count)=j;

            if flagplot==1 %plot
                xid=[idx(1) idx(end) idx(end) idx(1) idx(1) ];
                yid=[idy(1) idy(1) idy(end) idy(end) idy(1) ];
                hold on;plot(xid,yid,'m.-')
            end
            if count2 >=maxc
                flag=1;
                break
            end
        end
end %j
     if flag==1
         break
     end
end %i

%check non slumps
if 0
idy=1:32;idx=1:(32+(j-1)*32);
figure;imagesc(data0.x(idx)*1e-3,data0.y(idy)*1e-3,data0.z(idy,idx))
colormap jet; colorbar
caxis([-10 10])
% plot slump polygons
for i=1:n
    [s1x,s1y]=polarstereo_fwd(S1(i).Y,S1(i).X,[],[],70,-45);
    hold on;plot(s1x*1e-3,s1y*1e-3,'k.-','linewidth',3)
end
end
            
%numImageCategories = 2;

categories(trainingLabels1)
% categories(testLabels)

if flagplot==1 % plot 
    
figure
thumbnails = trainingImages1(:,:,:,1:28);
montage(thumbnails);title('slumps')
colorbar;colormap jet;caxis([-10 10])

figure;
for j=1:28
        subplot(4,7,j)
        imagesc(trainingImages1(:,:,:,j));
        axis equal;colormap jet;caxis([-10 10]);axis off
end
title('slumps')

k=countslump+1;
figure
thumbnails = trainingImages1(:,:,:,k:k+40-1);
montage(thumbnails);title('non slumps')
colorbar;colormap jet;caxis([-10 10])

figure;
for j=k:k+40-1
        subplot(5,8,j-countslump)
        imagesc(trainingImages1(:,:,:,j));
        axis equal;colormap jet;caxis([-10 10]);axis off
end
title('non slumps')
end % if plot



return
end