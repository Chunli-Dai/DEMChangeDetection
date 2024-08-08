function [trainingImages1,trainingLabels1,slumpbox]=prepareDetectron2(shapefiles,imagefiles)
% % Prepare training images for Detectron2 (upgraded MaskRCNN)
% Given shapefiles, and imagefiles to get training image matrices.
% shapefiles: shapefile of known slump outlines.
% imagefiles: elevation change files list
% output: slumpbox: the boundary of slumps

constant

trainingLabels1=categorical({'slump', 'non slump'}');
slumpbox=[];trainingImages1=[];

flagplot=0; %1 plot;0 not plot
flagaug=1; %1, if use data augment; 0, do not use augment
maxtrain=1004-189; %1000*0.9; %50000*0.9; %1725-173;%461; %917; %1725-173;%9e9; %maximum number of images used as training images, and rest will be used as validation.
flagorthotif=0; % save mosacked ortho image in Geotiff
flagrescale=0;  %1, fixed range [e.g. -20 10] to keep the information of positive/negative.
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

[demdir,name,ext] =fileparts([strtrim(regiondir)]);
regid=name(11:12);
regionstr=['region',regid]; %'region09';
if flagaug==0
odir=['slumps_',regionstr,'/'];
else
odir=['slumps',regionstr,'aug/'];
end

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

%if flagorthotif==1
	%save orthoimage in geotiff 
 	odirt1=['slumps_ortho/'];	
	% 3band geotiff files , dh, dem curvature *1e6, (T2 -2000/1/1)/365.25*10
	if flagaug==0
          odirt1=['./slumptif_',regionstr]; %e.g., slumptif_region09
	else
          odirt1=['./slumptif_',regionstr,'aug/']; %e.g., slumptif_region09aug
	end
	if ~exist(odirt1,'')
	    mkdir(odirt1)
	end
%end
        

% S1=shaperead('/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Ashley/shapefiles/Slumps_2017_18_SCAR_FINAL.shp');
S1=shaperead(shapefiles{1});% figure;mapshow(S1)
%/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_06_21_jump.tif

flagband=3; %3; %1 writing 1 band data for training images; 3 writing 3 band data; 2: elevation change and DEM curvature, two bands. 4: 2 band: elevation change and time of change. 
flagdata=3; %1 use the merged elevation change file; 2; use elevation change file in each sub tile.
	    %3 used the saved files, e.g., e.g., slumptif_region09_vl/slump447.tif
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
    data0.z=zeros(ny,nx);
end
resr2m=2; % 2m resolution;

%get the f fdir list
%load mat0.mat
matfile=[regiondir,'/mat0.mat'];
load(matfile)

%Get slumps
n=length(S1);
count=0; %clear trainingImages1 trainingLabels1 collecti
countaug=0; %count with data augmentation
xot={}; yot={}; fileot={}; 
xov={}; yov={}; fileov={}; 
xote={}; yote={}; fileote={}; %test
count_tr=0;count_v=0;count_te=0;

buff=100; %buffer the polygon by 100m

%80 Different number of images yield incoherent mosaic.
%27 half of the polygon in old version(bp1). good in new version
%184 long shape
% 904 dark: anomaly in change map
for k=1:length(shapefiles)
S1=shaperead(shapefiles{k});
n=length(S1);
for i=1:n
    if isfield(S1, 'zone')
    %if contains(S1(i).zone,'debris');continue;end %only detect scars
    %if ~(contains(S1(i).zone,'scar')|| contains(S1(i).zone,'nois'));continue;end   % two classes
    if ~(contains(S1(i).zone,'scar'));continue;end 
    end
    ll=S1(i).BoundingBox(1,:);ru=S1(i).BoundingBox(2,:);
    [llx,lly]=polarstereo_fwd(ll(2),ll(1),[],[],70,-45); %lat lon
    [rux,ruy]=polarstereo_fwd(ru(2),ru(1),[],[],70,-45); %lat lon
    [rlx,rly]=polarstereo_fwd(ll(2),ru(1),[],[],70,-45); %right low
    [lux,luy]=polarstereo_fwd(ru(2),ll(1),[],[],70,-45); %left up
    x=[llx,rlx,rux,lux];y=[lly,rly,ruy,luy];
    rangi=[min(x)-buff max(x)+buff min(y)-buff max(y)+buff];
    
    %Get multiple the polygons within the buffer zone
    x0=[rangi(1) rangi(2) rangi(2) rangi(1) rangi(1) ];
    y0=[rangi(4) rangi(4) rangi(3) rangi(3) rangi(4) ];
    [lat,lon]=polarstereo_inv(x0,y0,[],[],70,-45);
    bbox=[min(lon) min(lat); max(lon) max(lat)]; %left bottom ; right top
    tic;S1sub=shaperead(shapefiles{k},'BoundingBox',bbox);toc %10 minutes

    %gather catagory ids.
    idd=[];category_idsub=zeros(length(S1sub),1);
    for i1=1:length(S1sub)
	if contains(S1sub(i1).zone,'scar')
	    category_idsub(i1)=0;
%	elseif contains(S1sub(i1).zone,'nois') %two classes
%	    category_idsub(i1)=1;
	else
	    idd=[idd(:);i1];
	end
    end
    S1sub(idd)=[];category_idsub(idd)=[];
	if length(S1sub) ==0; continue;end

%     data=readGeotiff(changefile,'map_subset',rangi);
    idx=find(data0.x>=rangi(1)&data0.x<=rangi(2));
    idy=find(data0.y>=rangi(3)&data0.y<=rangi(4));
    if flagdata==1
      data.x=data0.x(idx);data.y=data0.y(idy);
      data.z=data0.z(idy,idx);
    elseif flagdata==2
     fprintf(['\n Working on i polygon:',num2str(i),'\n'])
	close all
     %Given a box (rangi), to get the mosaicked files
     [datao]=box2mosaic(rangi,flagband,f,fdir,i);
%    [datao]=box2mosaic(rangi,3,f,fdir); %test hi
     data=datao;data.z(:,:,2:end)=[]; %to be compatible with older version
     M1=isnan(data.z);
     ratio=sum(M1(:))/length(data.z(:))*100; %ratio of void pixels
     if ratio>10 % more than 10% of the box has no data.
	%fprintf(['\n There is no elevation change data for this slump! \n '])
	fprintf(['\n There are more than 10%% of the box has no elevation change data. Skip this slump! \n '])
	continue
     end
     	
    elseif flagdata==3
      datao=getslumptif(S1(i));
      data=datao;data.z(:,:,2:end)=[]; %to be compatible with older version
      M1=isnan(data.z);
      ratio=sum(M1(:))/length(data.z(:))*100; %ratio of void pixels
      if ratio>10 % more than 10% of the box has no data.
        %fprintf(['\n There is no elevation change data for this slump! \n '])
        fprintf(['\n There are more than 10%% of the box has no elevation change data. Check this slump: id=',num2str(i),' ',shapefiles,' \n '])
        continue
      end
    end %flagdata==1

    [mi,ni]=size(data.z); minm=min([mi,ni]);
    if (minm<3);continue
    else
        
        count=count+1;
	    countaug=countaug+1;
                
        %For Detectron2
        %Generate filenames,filesize,XYb
        %Slump polygon
        %change data to 2 m resolution;-> 2m yields better accuracy for Detecton2
	%move this part to box2mosaic.m
%       data2m.x=data.x(1):resr2m:data.x(end);
%       data2m.y=data.y(1):(-resr2m):data.y(end);
%       data2m.z=interp2(data.x,data.y,double(data.z),data2m.x,data2m.y','*linear');
%       data=data2m;
        resrx=mean(diff(data.x));resry=mean(diff(data.y));
        if 0 %single polygon
            [s1x,s1y]=polarstereo_fwd(S1(i).Y,S1(i).X,[],[],70,-45);
            Xb=round((s1x-data.x(1))/resrx)+1; 
            Yb=round((s1y-data.y(1))/(resry))+1; %descending        
            XYbi=[Xb(:),Yb(:)]; %n by 2
            id1=find(any(isnan(XYbi),2));
            XYbi(id1,:)=[];
        else %multiple polygons
            XYbi=[]; %cell
	    category_idg=[];% array
            nsub=length(S1sub);idd=[];
	    countmp=0; %count of multiple polygons
            for i1=1:nsub
            [s1x,s1y]=polarstereo_fwd(S1sub(i1).Y,S1sub(i1).X,[],[],70,-45);
            Xb=round((s1x-data.x(1))/resrx)+1; 
            Yb=round((s1y-data.y(1))/(resry))+1; %descending        
            XYbi1=[Xb(:),Yb(:)]; %n by 2
            id1=find(any(isnan(XYbi1),2));
            XYbi1(id1,:)=[]; 
            %XYbi{i1}=XYbi1;
         
            %get rid of the out of boundary points.
            [nyi,nxi,nzi]=size(data.z);
            maskorg=false(size(nyi,nxi));
        
            maskorg_k=poly2mask(XYbi1(:,1),XYbi1(:,2), nyi,nxi); 
            maskorg=maskorg|maskorg_k;
       
	    Mf4clean = bwareaopen(maskorg, 25); %remove small clusters<25*4 m^2
            B = bwboundaries(Mf4clean); 

            for j1=1:length(B)
            xid=B{j1}(:,2); yid=B{j1}(:,1);
	    %avoid bug ValueError: Failed to use mask_format=='polygon' from the given annotations!
	    if (max(xid)-min(xid))==0 || (max(yid)-min(yid))==0; continue;end  %point or line
	    countmp=countmp+1;
	    %XYbi{countmp}=[xid(:),yid(:)];
          
	    if 0 
	    %use negative/positive pairs for training;
            %find the matching positive clusters.
	    maskjm=data;maskjm.z=Mf4clean;
	    Mpc2=data;Mpc2.z=data.z>0.5; %> 1m;
            [~, Mout]=matchclusters(maskjm,Mpc2);
	    Bnp=bwboundaries(Mout.z);
	    %get the outlines.
            j11=1; xid=Bnp{j11}(:,2); yid=Bnp{j11}(:,1);
	    end % if 0

	    XYbi{countmp}=[xid(:),yid(:)];
            category_idg(countmp)=category_idsub(i1);
            end
	    end %i1
        
        end
        
        filenamei=['slump',num2str(count),'.jpg']; %jpg format is compressed!! try imwrite(A1,'test4.jpg','Mode','lossless');
							%cause error in Detectron2: 'NoneType' object has no attribute 'shape' 
	%filenamei=['slump',num2str(count),'.tif']; %tiff format is not compressed; ->yields 
	fprintf(['\n Output filename: ',filenamei,' \n']);

	n1=length(XYbi);
        lati=[];loni=[];
	for i1=1:n1
        xid=XYbi{i1}(:,1);yid=XYbi{i1}(:,2);
        s1x=(xid-1)*resrx+data.x(1);s1y=(yid-1)*resry+data.y(1);
        [lati1,loni1]=polarstereo_inv(s1x,s1y,[],[],70,-45);
	lati{i1}=lati1; loni{i1}=loni1;
	end


	%collect polygon shapefiles
        if count <= maxtrain %count<=27
%       if  rux >=-2485300 %for peel, right 70%
        dirfilename=[odir,'train/',filenamei];
	for i1=1:length(XYbi)
	count_tr=count_tr+1;
	%xot{count_tr}=S1(i).X;yot{count_tr}=S1(i).Y; %input shapefile
	%connected negatvie/positive clusters
	xot{count_tr}=loni{i1};yot{count_tr}=lati{i1};
	fileot{count_tr}=dirfilename;
	end % i1
%       elseif rux >=-2495300
	else
        dirfilename=[odir,'val/',filenamei];
	for i1=1:length(XYbi)
	count_v=count_v+1;
%	xov{count_v}=S1(i).X;yov{count_v}=S1(i).Y;
        xov{count_v}=loni{i1};yov{count_v}=lati{i1};
	fileov{count_v}=dirfilename;
        end % i1
%	else %test
%         dirfilename=[odir,'test/',filenamei];
%	  if contains(S1(i).zone,'scar') %only test scar, no noise
%       for i1=1:length(XYbi)
%	count_te=count_te+1;
%	xote{count_te}=S1(i).X;yote{count_te}=S1(i).Y;%fileote{count_te}=dirfilename;
%       xote{count_te}=loni{i1};yote{count_te}=lati{i1};
%	fileote{count_te}=dirfilename;
%	end % i1
%	% else; continue
%	  end
        end

        %convert the double matrix to uint16 matrix with values between 0 and 65535.
%         A1=uint16(rescale(data.z,0,65535));
%         imwrite(A1,dirfilename,'BitDepth',16);
        %Detectron2 seems to be only work with uint8 images. Need to be checked.
	if flagrescale==0
        A1=uint8(rescale(data.z(:,:),0,255));
	elseif flagrescale==1
	    inmin=-20;inmax=10; %truncate data to this range [inmin inmax], then rescale.
	    A1=uint8(rescale(data.z(:,:),0,255,'InputMin',inmin,'InputMax',inmax));
	end
%       A1=uint16(rescale(data.z(:,:),0,65535));%uint16 yields bad results!! Detectron2 works with uint8 [0 255]!
	if flagband==1
          %A1=uint8(rescale(data.z(:,:),0,255));
%         A1=uint16(rescale(datao.z(:,:,3),0,65535)); 
%         A1=uint8(rescale(datao.z(:,:,3),0,255)); %test hi
	elseif flagband==3
	  %get orthoimage and DEM	
	  dem.z=datao.z(:,:,2);ortho.z=datao.z(:,:,3);
          %A1(:,:,2)=uint8(rescale(dem.z(:,:),0,255)); %DEM;
          %A1(:,:,3)=uint8(rescale(ortho.z(:,:),0,255)); %image;
          A1(:,:,2)=uint8((dem.z(:,:))); % %curvature
          A1(:,:,3)=uint8((ortho.z(:,:))); % time map

	%save orthoimage in geotiff 
	if flagorthotif==1
	filenamet1i=[odirt1,'/','slump',num2str(count),'.tif'];
	projstr='polar stereo north';
	writeGeotiff(filenamet1i,datao.x,datao.y,uint16(ortho.z),2,0,projstr)
	end
	
	elseif flagband==2 
	  %get curvature
	  totalc=datao.z(:,:,2);
	%  A1(:,:,2)=uint8(rescale(totalc(:,:),0,255)); %DEM; one singular value may reduce the value of other pixels.
	  A1(:,:,2)=uint8((totalc(:,:))); %DEM;
	  A1(:,:,3)= A1(:,:,1); %  %Data with 2 components not supported for JPEG files. So add 3rd band.

	if 0 && count==6 %debug
	t1=uint16(A1(:,:,2));
	nanmax(totalc(:))
	max(max(A1(:,:,2)))
	max(t1(:))
	save t6.mat -v7.3
	exit
	end

	%save total curvature in geotiff 
	if 1
	filenamet1i=[odirt1,'/','totalc',num2str(count),'.tif'];
	projstr='polar stereo north';
	writeGeotiff(filenamet1i,datao.x,datao.y,uint16(A1(:,:,2)),2,0,projstr)
	end

	elseif flagband ==4
	   T2=datao.z(:,:,2);
	   A1(:,:,2)=uint8((T2(:,:))); %
          A1(:,:,3)= A1(:,:,1); %
	
	end

	filenamet1i=[odirt1,'/','slumpm',num2str(count),'.tif']; %multi-spectral band
	projstr='polar stereo north';
	writeGeotiff(filenamet1i,datao.x,datao.y,double(datao.z),4,nan,projstr)

        imwrite(A1,dirfilename);  %tif format
%       imwrite(A1,dirfilename,'Mode','lossless'); %jpg %cause error in Detectron2: 'NoneType' object has no attribute 'shape' 

        str1=imfinfo(dirfilename);
        filesizei=str1.FileSize;

        filenames{countaug}=filenamei;
        filesize(countaug)=filesizei;
        XYb{countaug}=XYbi;
	    category_idgc{countaug}=category_idg;
        absfilenames{countaug}=dirfilename;
        
        if flagplot==1
            figure;
            %plot image in coordinates of pixel ids.
            figure;imagesc(data.z);colorbar%imshow(A);
            caxis([-3 3]);colormap jet
            hold on;plot(Xb,Yb,'k.-')
            title(dirfilename)
            pause
            close all
        end

	%Apply Data augmentation to both image and polygons;
	if flagaug==1
	  id1=countaug;
 	  [filenames,filesize,XYb,absfilenames,countaug]=applyaug(filenames,filesize,XYb,absfilenames, countaug,dirfilename,A1,XYbi);
	  id2=countaug;
	  category_idgc(id1:id2)={category_idg};
	end %if flagaug
        
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

%             figure;imagesc(imagei)
            figure;mapshow(S1(i))
            
            %plot image in coordinates of pixel ids.
            A2=imread(dirfilename);
            figure;imagesc(A2);colorbar%imshow(A);
            hold on;plot(Xb,Yb,'r.-')
        end %if 0 %for plotting 
        
    end
	close all
end
end %for k

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

save test1.mat -v7.3

%% Creat json files
ofile=[odir,'train/via_region_data.json'];
nall=length(filenames);
%maxtrain=round(0.9*nall);
ntr=min(nall,maxtrain);
%id=1:ntr;
id=find(contains(absfilenames,'train')); 
filenames1=filenames(id);filesize1=filesize(id);XYb1=XYb(id); category_idgc1=category_idgc(id);
%remove empty polygons
M=false(length(XYb1),1);
for ti=1:length(XYb1); if isempty(XYb1{ti}); M(ti)=1;end;end
fprintf(['\n ', num2str(nansum(M==1)), ' polygons will be deleted due to voidness (by translation).'])
filenames1(M)=[];filesize1(M)=[];XYb1(M)=[];category_idgc1(M)=[];id(M)=[];

creatjson(filenames1,filesize1,XYb1,ofile,category_idgc1);

ofile=[odir,'val/via_region_data.json'];
%id=ntr+1:nall;
id=find(contains(absfilenames,'val')); 
filenames1=filenames(id);filesize1=filesize(id);XYb1=XYb(id); category_idgc1=category_idgc(id);

M=false(length(XYb1),1);
for ti=1:length(XYb1); if isempty(XYb1{ti}); M(ti)=1;end;end
fprintf(['\n ', num2str(nansum(M==1)), ' polygons will be deleted due to voidness (by translation).'])
filenames1(M)
filenames1(M)=[];filesize1(M)=[];XYb1(M)=[];category_idgc1(M)=[];id(M)=[];

if length(id)>=1
creatjson(filenames1,filesize1,XYb1,ofile,category_idgc1);
end

return
end
