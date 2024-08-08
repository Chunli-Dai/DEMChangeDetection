function [Precision,Recall]=poly2precision(Sml,S2)
%Given two shapefiles, get the precision and recall
%S2 is ground truth
%Sml is the shapefiles from machine learning

% S2=shaperead('/Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/maskgroundtruth.shp');% figure;mapshow(S1)
% Sml=shaperead('maskthres90_1552train_band1jpg.shp');

constant  %get smlarea
smlarea=1e3*2; %1000 m^2

if ischar(Sml)
    flagfmt=1; %input is a shapefile
    file=Sml;
    Sml=shaperead(file);
else isstruct(Sml) %input is a struct of mask matrix
    flagfmt=2;
end

IoUthre=0.1; %0.5; %0.1 ^;  %0.4;  0.5 0.6  
buff=1e3; %assume the ground truth polygon that overlaps with the target is no wider than 1 km
resr=10; %2; %2m grid size to creat mask from given polygon
resrmask=2; % 66sec (10m) vs 64 sec (2m)

if flagfmt==1 %input is a shapefile

TPFP=length(Sml);
s1xg=[];s1yg=[];
for i=1:length(Sml) %target shapefiles
    [s1x,s1y]=polarstereo_fwd(Sml(i).Y,Sml(i).X,[],[],70,-45);
    s1x(isnan(s1x))=[];s1y(isnan(s1y))=[];
    s1xg=[s1xg(:);s1x(:)];s1yg=[s1yg(:);s1y(:)];
end

    rang0=round([min(s1xg) max(s1xg) min(s1yg) max(s1yg)]/resr)*resr;
    mask1.x=rang0(1):resr:rang0(2);mask1.y=rang0(4):-resr:rang0(3);
    maskx=mask1.x;masky=mask1.y;
    ny=length(mask1.y);nx=length(mask1.x);
    mask1.z=false(ny,nx);

for i=1:length(Sml) %target shapefiles
    [s1x,s1y]=polarstereo_fwd(Sml(i).Y,Sml(i).X,[],[],70,-45);
    s1x(isnan(s1x))=[];s1y(isnan(s1y))=[];

    Xb1=round((s1x-maskx(1))/resr)+1; 
    Yb1=round((s1y-masky(1))/(-resr))+1; %descending    
    maski=poly2mask(Xb1,Yb1, ny,nx);     
    mask1.z=mask1.z|maski; 
    
%     %tmp
%     figure;imagesc(mask1.x,mask1.y,mask1.z);colorbar
%     hold all; plot(s1x,s1y,'r.-')
%             rang0=round([min(s1x)-100 max(s1x)+100 min(s1y)-100 max(s1y)+100]);
%     axis(rang0)
%     pause
end

%%Poly2mask in dividemapsub.m may not fully recover the original mask, and unintentionally separate clusters.
Mf4clean = bwareaopen(mask1.z, round(smlarea/resr/resr)); %remove small clusters
mask1.z=Mf4clean;

CC = bwconncomp(mask1.z);
TPFP=CC.NumObjects;

elseif flagfmt==2 %input is a mask matrix
    
mask1=Sml;clear Sml
resx=mean(mask1.x(2:end)-mask1.x(1:end-1));resy=mean(mask1.y(2:end)-mask1.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);
if resr < resrmask
    %reduce mask to 10m resolution for faster computation
    mask1org=mask1;
    mask1.x=min(mask1org.x):resrmask:max(mask1org.x);
    mask1.y=max(mask1org.y):-resrmask:min(mask1org.y);
    mask1.z = interp2(mask1org.x,mask1org.y,mask1org.z,mask1.x,mask1.y','*nearest',0);    
    resr=resrmask;
end

Mf4=mask1.z;mask1.z=[];
Mf4clean = bwareaopen(Mf4, round(smlarea/resr/resr)); %remove small clusters
mask1.z=Mf4clean;

BW=Mf4clean;clear Mf4clean
[n1,m1]=size(BW);
if min([n1,m1])> 100e3
    %reduce matrix size
    BW=BW(1:2:end,1:2:end);
end
CC = bwconncomp(BW);
TPFP=CC.NumObjects;
end %flagfmt

lenS2=length(S2);
for j=1:lenS2 %ground truth
        [s2x,s2y]=polarstereo_fwd(S2(j).Y,S2(j).X,[],[],70,-45);
        s2x(isnan(s2x))=[];s2y(isnan(s2y))=[];
        s2xg{j}=s2x;s2yg{j}=s2y;
end

FP=0;TP=0;FN=0;
idp=[]; %ids of all ground truth slumps that paired with target slumps.
s2xp=[];s2yp=[]; %for plot
ngto=0; %number of ground truth that outside of study area.
    for j=1:lenS2 %ground truth
%         [s2x,s2y]=polarstereo_fwd(S2(j).Y,S2(j).X,[],[],70,-45);
%         s2x(isnan(s2x))=[];s2y(isnan(s2y))=[];
        s2x=s2xg{j};s2y=s2yg{j};
        resr1m=1; %for small temporary areas; better resolution yields better mask from polygon 
        rang0=round([min(s2x)-buff max(s2x)+buff min(s2y)-buff max(s2y)+buff]/resr1m)*resr1m;
        maskx=rang0(1):resr1m:rang0(2);masky=rang0(4):-resr1m:rang0(3);
        ny=length(masky);nx=length(maskx);
        x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
        
        Xb1=round((s2x-maskx(1))/resr1m)+1; 
        Yb1=round((s2y-masky(1))/(-resr1m))+1; %descending    
        maskj=poly2mask(Xb1,Yb1, ny,nx);     
               
        %crop mask1 for faster interpolation
        rang2=round([min(s2x)-2*buff max(s2x)+2*buff min(s2y)-2*buff max(s2y)+2*buff]/resr1m)*resr1m;
        Mx=mask1.x>=rang2(1)&mask1.x<=rang2(2);My=mask1.y>=rang2(3)&mask1.y<=rang2(4);
        mask1r.x=mask1.x(Mx);mask1r.y=mask1.y(My);
        mask1r.z=mask1.z(My,Mx);
        
        [n1,m1]=size(mask1r.z);
        
        if min([n1,m1])<2
            fprintf(['\n Groundtruth is not within the boundary of testing images, remove this, at j=:',num2str(j)]);
            %remove this slump groundtruth list.
            ngto=ngto+1;
            continue
        end
        
        try
        maskk=interp2(mask1r.x,mask1r.y,mask1r.z,maskx,masky','*nearest',0); 
        catch e
        fprintf(['Error at j=:',num2str(j)]);
        end
        maskkj=maskk&maskj;
        
        overlap=sum(sum(maskk&maskj)); %pixels

        if overlap>0
            idp=[idp(:);j]; %id of ground truth that intersects with target.

            % get another box rang1 larger than rang0;
            rang1=round([min(s2x)-1.1*buff max(s2x)+1.1*buff min(s2y)-1.1*buff max(s2y)+1.1*buff]/resr1m)*resr1m;
            maskx1=rang1(1):resr1m:rang1(2);masky1=rang1(4):-resr1m:rang1(3);
            ny1=length(masky1);nx1=length(maskx1);
            %convert the box of rang0 to mask in rang1.
            Xb1=round((x0-maskx1(1))/resr1m)+1; 
            Yb1=round((y0-masky1(1))/(-resr1m))+1; %descending    
            mask_box=poly2mask(Xb1,Yb1, ny1,nx1); 
            maskj_big=interp2(maskx,masky,maskj,maskx1,masky1','*nearest',0);
            maskk_big=interp2(mask1r.x,mask1r.y,mask1r.z,maskx1,masky1','*nearest',0);
            
%             CC = bwconncomp(maskkj);
%             nA=CC.NumObjects; %total number of polygons that intersects with the ground truth;
            %recover the polygons in the original matrix;
%             CC = bwconncomp(maskk);
%             maskrecover=false(size(mask1.z));
%             maskks=false(size(maskk)); % all polygons that intersect with the ground truth
            CC = bwconncomp(maskk_big);
            maskks=false(size(maskk_big)); % all polygons that intersect with the ground truth
            for ki=1:CC.NumObjects
%                 BW3=false(size(maskk));
%                 BW3(CC.PixelIdxList{ki})=maskk(CC.PixelIdxList{ki});
%                 overlap_ki=sum(sum(BW3&maskj));
                BW3=false(size(maskk_big));
                BW3(CC.PixelIdxList{ki})=maskk_big(CC.PixelIdxList{ki});
                overlap_ki=sum(sum(BW3&maskj_big));
                if overlap_ki>0
                    maskks=maskks|BW3;
%                     mask_kir=interp2(maskx,masky,double(BW3),mask1.x,mask1.y','*nearest',0);
%                     maskrecover=maskrecover|mask_kir;
                end
            end
%             mask_rest=mask1.z-maskrecover;
%             CC = bwconncomp(mask_rest);
%             nB=CC.NumObjects; %rest of the polygons that do not intersects with the ground truth or cropped out of the box.
%             ApB=nA+nB;
            
%             s2xp=[s2xp(:);nan;s2x(:)];s2yp=[s2yp(:);nan;s2y(:)];
            
            %check if the all polygons that intersect with ground truth is within the box of rang0.
%             if ApB > TPFP
%                 %or caused by pixel difference from interpolation to different resolution
%                 fprintf(['\n Attention: the ground truth polygon extends beyond the buffer zone! Increase the parameter buff! \n'])
%             end

            %check if  all polygons (maskks) that intersect with ground truth is within the box of rang0.
            
            mask_rest=maskks&(~mask_box); %the part of maskks that's outsize of the box rang0.
            if sum(mask_rest(:))>1
                fprintf(['\n Attention: the ground truth polygon extends beyond the buffer zone! Increase the parameter buff! \n'])
            end
            
            overlapi=overlap; 
%             union=sum(sum(maskks|maskj));
            union=sum(sum(maskks|maskj_big));
            IoU=overlapi/union;

            if IoU>=IoUthre
                TP=TP+1;
            else
%                 FP=FP+1;
            end
            
        else %overlap <= 0
            FN=FN+1;
 
        end
    end
    
idp=unique(idp);
ng=length(S2)-ngto;
FN=ng-length(idp);
FP=TPFP-TP;

% FP=TPFP-TP;
if TPFP == (TP+FP)
    str1='Yes.';
else
    str1='No.';
end
fprintf(['\n Double check TPFP= TP+FP?  ',num2str(TPFP) ,' ?= ',num2str(TP),' + ',num2str(FP), ' ',str1]);

Precision=(TP)/(TP+FP)*100;
Recall=(TP)/(TP+FN)*100;
fprintf(['\n IoUthre=',num2str(IoUthre),'.\n'])
fprintf(['\n True Positive: ',num2str(TP),', False Positive: ',num2str(FP),', False Negative: ',num2str(FN), ', Total ground truth:',num2str(ng),'.\n'])
fprintf(['\n Precision: ',num2str(Precision),'%%; Recall: ',num2str(Recall),'%%\n'])
fprintf(['\n Total ground truth, original: ',num2str(length(S2)),'; number that is out of study area ',num2str(ngto),'\n'])


end
