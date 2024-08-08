
%infile1='WV02_20160827_103001005A735B00_103001005CAE7900_seg1_2m_meta.txt';
%%infile1='WV02_20160827214547_103001005CAE7900_16AUG27214547-M1BS-500890307070_01_P009_u16ns3413.xml';
%infile2='WV02_20100710213820_1030010005A70000_10JUL10213820-M1BS-052618580140_01_P004_u16ns3413.xml';
%infile1='WV02_20150714220327_10300100442BD900_15JUL14220327-M1BS-500515445020_01_P003_u16ns3413.xml';
%macdir='/Users/chunlidai/surge/';
%macdir='';
%addpath(genpath([macdir,'/home/dai.56/arcticdemapp/river/rivergithub2/']));
%load groad.mat

function [OutName1,OutName2,data1r,data2r]=prepareMJ(infile1,infile2,odir,groad,rang0,tag,ofmt)
% Modified from ~/arcticdemapp/river/rivergithub2/rivercore/coregimage/prepareMJ.m
% prepare the image to be tiff files with res=15m, same common grids, in same coordinates UTM.
% data.z be uint8, max=255 -> no need. MIMC works good with uint16.
% Seongsu suggest to rescale to  [0 1].
% accurate coordinate system
% check orthoimage and nadir-view image -> Aster images are also orthoimages.
% simplify the code. -> done by polygon
% bp2:
% bp3: ell2utm too slow, try another toolbox
% dataflag=1, infile1 is ASTER, infile2 is Worldview images
% dataflag=2, infile1 is ASTER, and infile2 is Aster images.
% dataflag=3, infile 1 and 2 are both digital globe images, input filename are metafiles.
% infile1: reference image metafile; infile2: target image metafile.
% tag: rock mask; 1 rock; 0 non rock;
% ofmt: output format: 1, Digital Numbers, uint16, (input of SETSM); 2, reflectance, double, (input of COSICorr)
macdir=[];
if 0
addpath([macdir,'/home/chunli/scripts/']);
addpath(genpath([macdir,'/home/chunli/scripts/Denali/geodetic299/']));
addpath('/Users/chunlidai/Google Drive/OpticalImages/ArcticDEM/Coregister');
addpath('/Users/chunlidai/Google Drive/OpticalImages/ArcticDEM/Coregister/MJ_coreg/');
addpath(genpath('/Applications/bda/Earthdata/ASTERDEM/'));
addpath(genpath('/Users/chunlidai/Google Drive/OpticalImages/MIMC2/MIMC2_demo_package_release_170412'));

zone=7;
infile1='AST14DMO_00302222002211959_20180119102946_29249_V3N.tif';
infile2='WV01_20160308_102001004D7D5100_102001004AAA7200_seg1_2m_ortho.tif';

infile2='AST14DMO_00302182003211315_20180201141944_8450_V3N.tif';
infile2='AST14DMO_00304112002211932_20180201141814_6623_V3N.tif';
end
constant

co=[];

% dataflag=2; %data type of infile2: 1 is Worldview images, and 2 is Aster images.
%flagplot=1;
znan=0;
wgsa=6378137.0;e=0.08181919;e2=e^2;%wgs84 in polarstereo_fwd.m
resr=15.; %for coregisteration
resr=2.; %for coregisteration

[~,filename1,ext1]=fileparts(infile1);name1=[filename1,ext1];
[~,filename2,ext2]=fileparts(infile2);name2=[filename2,ext2];
if strcmp(filename1(1:3),'AST') && strcmp(filename2(1:3),'AST') 
    dataflag=2; %data type of infile2: 1 is Worldview images, and 2 is Aster images.
elseif strcmp(filename1(1:3),'AST') && ~strcmp(filename2(1:3),'AST') 
    dataflag=1;
else %both images are not ASTER.
    dataflag=3;
end

if contains(name1,['_',demext]) %strcmp(name1,['_',demext]); %'_dem.tif')
    dataflag=4; %DEM files
end
  
%infile 1 is ASTER
if dataflag==1 || dataflag==2 
    
data=readGeotiff(infile1);
if flagplot==1
figure;imagesc(data.x,data.y,data.z);colorbar;axis equal
title('Image 1 in its own UTM coordinate')
end
[X1sv,Y1sv]=meshgrid(data.x,data.y);
%get zone
% metafile=[infile1(1:47),'.zip.met'];
metafile=strrep(infile1,'_V3N.tif','.zip.met');
c=textread(metafile,'%s','delimiter','\n');
% r=find(contains(c,'MAPZONECODE'));
r=find(~cellfun(@isempty,strfind(c,'MAPZONECODE')));
Xbs=c{r(1)+2};  %deblank(strrep(c{r(1)+2},'',''));
r = strfind(Xbs, '"');
zone = sscanf(Xbs((r(1)+1):(r(2)-2)), '%g', 1);

%Get the actual polygon boundary of image 1 that has valid data.
M=data.z~=znan;
Mb=false(size(M));
Mb(1,:)=M(1,:);Mb(end,:)=M(end,:);Mb(:,1)=M(:,1);Mb(:,end)=M(:,end);
Md1=imdilate(M,ones(3));
M=logical(Md1-M)|Mb; %include the horizontal boundary 
idx=X1sv(M); idy=Y1sv(M);
k=boundary(idx,idy,0); % 0 point boundary, 1 conv hull;

Xb1=idx(k);
Yb1=idy(k);
if flagplot==1;hold on;plot(Xb1,Yb1,'g>-');end
[nyov,nxov]=size(data.z);
idx=abs((Xb1-data.x(1))/(data.x(2)-data.x(1)))+1;
idy=abs((Yb1-data.y(1))/(data.y(2)-data.y(1)))+1;
Mb1 = poly2mask(idx,idy, nyov,nxov);

% boundary
[LAT,LON]=utm2ell( Yb1(:),Xb1(:), zone,wgsa,e2);
LAT=LAT*180/pi;LON=LON*180/pi;

if dataflag==1 %data type of infile2: 1 is Worldview images, and 2 is Aster images.
    
%Transform boundary of image 1 in image 2's coordinates
[x1,y1]=polarstereo_fwd(LAT,LON,[], [],70,-45);
ranget=[[min(x1(:)) max(x1(:)) min(y1(:)) max(y1(:))]/resr];
ranget=[floor(ranget(1)) ceil(ranget(2)) floor(ranget(3)) ceil(ranget(4))]*resr;

data2=readGeotiff(infile2,'map_subset', ranget);

if flagplot==1
figure;imagesc(data2.x,data2.y,data2.z);colorbar;axis equal
title('Image 2 in its own polar sterographic coordinates')
end
[X,Y]=meshgrid(data2.x,data2.y);

%Get the actual polygon boundary of image 2 that has valid data.
M=data2.z~=znan;
Mb=false(size(M));
Mb(1,:)=M(1,:);Mb(end,:)=M(end,:);Mb(:,1)=M(:,1);Mb(:,end)=M(:,end);
Md1=imdilate(M,ones(3));
M=logical(Md1-M)|Mb; %keep the horizontal boundary 
idx=X(M); idy=Y(M);
k=boundary(idx,idy,0); % 0 point boundary, 1 conv hull;

%Transform the boundary of image 2 to UTM ZONE 1
X2=idx(k);
Y2=idy(k);
[LAT2,LON2]=polarstereo_inv(X2,Y2,[], [],70,-45);

lcm=deg2rad(zone*6-183);
[Yb2,Xb2,zoneo]=ell2utm(LAT2*pi/180,LON2*pi/180,wgsa,e2,lcm);

idx=round(abs((Xb2-data.x(1))/(data.x(2)-data.x(1))))+1;
idy=round(abs((Yb2-data.y(1))/(data.y(2)-data.y(1))))+1;
Mb2 = poly2mask(idx,idy, nyov,nxov);
M=Mb1&Mb2;
rangeov=[min(X1sv(M)) max(X1sv(M)) min(Y1sv(M)) max(Y1sv(M))];

%crop the overlapping data
ranget=rangeov/resr;
ranget=[floor(ranget(1)) ceil(ranget(2)) floor(ranget(3)) ceil(ranget(4))]*resr;
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);

if flagplot==1
hold on;plot([ranget(1) ranget(2) ranget(2) ranget(1) ranget(1)],[ranget(3) ranget(3) ranget(4) ranget(4) ranget(3)],'r-','linewidth',4)
end

%interpolate ASTERDEM to regular 15 m overlapping grid
znan=0;
%data.z(data.z == znan) = NaN; % %convert nodata values to nan for interpolation
%NaN is a double, can't assign it to a uint8
tz1=double(data.z);tz1(tz1==znan)=NaN;data.z=tz1;
tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
tz(isnan(tz))=znan; %return toznan 
data1r= struct();
%data1r.x=tx;data1r.y=ty;  data1r.z=uint8(tz./max(tz(:))*255);
%data1r.x=tx;data1r.y=ty;  data1r.z=uint16(tz);

% custmorized scale %[min1 max1] to [0 max2] 
max2=1;
min1=min(tz(tz~=znan)); max1=max(tz(tz~=znan));
a=max2./(max1-min1);b=-a*min1;
tzs=tz*a+b;
tzs(tz==znan)=0;
data1r.x=tx;data1r.y=ty;  data1r.z=tzs;

if flagplot==1
figure;imagesc(data1r.x,data1r.y,double(data1r.z));colorbar;axis equal
hold on;plot([ranget(1) ranget(2) ranget(2) ranget(1) ranget(1)],[ranget(3) ranget(3) ranget(4) ranget(4) ranget(3)],'r-','linewidth',4)
title('Overlapped Image 1 in its own coordinates')
end

%get the overlapping grid in polar sterographic coordinates
[X,Y]=meshgrid(tx,ty);
%utm 2 ell
[LAT,LON]=utm2ell( Y,X, zone,wgsa,e2);
[x1,y1]=polarstereo_fwd(LAT*180/pi,LON*180/pi,[], [],70,-45);

if flagplot==1
hold on;plot(x1(1:2000:end),y1(1:2000:end),'ro')
end

%% % second data type
elseif dataflag==2 %data type of infile2: 1 is Worldview images, and 2 is Aster images.
%get zone 2
metafile=strrep(infile2,'_V3N.tif','.zip.met')
c=textread(metafile,'%s','delimiter','\n');
% r=find(contains(c,'MAPZONECODE'));
r=find(~cellfun(@isempty,strfind(c,'MAPZONECODE')));
Xbs=c{r(1)+2};  %deblank(strrep(c{r(1)+2},'',''));
r = strfind(Xbs, '"');
zone2 = sscanf(Xbs((r(1)+1):(r(2)-2)), '%g', 1);

%transform the boundary of image 1 to the UTM coordinates of image 2.
lcm=deg2rad(zone2*6-183);
[y1,x1,zoneo]=ell2utm(LAT*pi/180,LON*pi/180,wgsa,e2,lcm);

ranget=[[min(x1(:)) max(x1(:)) min(y1(:)) max(y1(:))]];

data2=readGeotiff(infile2,'map_subset', ranget);
[X,Y]=meshgrid(data2.x,data2.y);

% filter out clouds
maskfile=strrep(infile2,'_V3N.tif','_mask.mat');
if 0
mp=roipoly;
save( maskfile, 'mp');
elseif 0
    load maskfile
    data2.z(mp)=0;
end

%Get the actual polygon boundary of image 2 that has valid data.
M=data2.z~=znan;
Mb=false(size(M));
Mb(1,:)=M(1,:);Mb(end,:)=M(end,:);Mb(:,1)=M(:,1);Mb(:,end)=M(:,end);
Md1=imdilate(M,ones(3));
M=logical(Md1-M)|Mb; %get the horizontal boundary 
idx=X(M); idy=Y(M);
k=boundary(idx,idy,0); % 0 point boundary, 1 conv hull;

%Transform the boundary of image 2 in UTM zone2 to UTM ZONE 1
X2=idx(k);
Y2=idy(k);
if flagplot==1; figure;imagesc(data2.x,data2.y,data2.z);colorbar;axis equal
hold on;plot(X2,Y2,'g>-'); 
title('Image 2 in its own UTM coordinates')
end

[lat,lon]=utm2ell(Y2(:),X2(:), zone2,wgsa,e2);
lcm=deg2rad(zone*6-183);
[y1,x1,zoneo]=ell2utm(lat,lon,wgsa,e2,lcm);
XYUTM2=[x1(:),y1(:)];

if 0 % inpolygon slow
X1svr=X1sv(1:100:end,1:100:end);Y1svr=Y1sv(1:100:end,1:100:end);
in=inpolygon(X1svr,Y1svr,XYUTM2(:,1),XYUTM2(:,2));
in=logical(imresize(int16(in),size(X1sv)));
M2=(data.z~=znan);
rangeov=[min(X1sv(in&M2)) max(X1sv(in&M2)) min(Y1sv(in&M2)) max(Y1sv(in&M2))];
else %faster
Xb2=XYUTM2(:,1); Yb2=XYUTM2(:,2);
idx=round(abs((Xb2-data.x(1))/(data.x(2)-data.x(1))))+1;
idy=round(abs((Yb2-data.y(1))/(data.y(2)-data.y(1))))+1;
Mb2 = poly2mask(idx,idy, nyov,nxov);
M=Mb1&Mb2;
rangeov=[min(X1sv(M)) max(X1sv(M)) min(Y1sv(M)) max(Y1sv(M))];
end

%crop the overlapping data
ranget=rangeov/resr;
ranget=[floor(ranget(1)) ceil(ranget(2)) floor(ranget(3)) ceil(ranget(4))]*resr;
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);

if flagplot==1
figure (1); 
hold on;plot(X1sv(M),Y1sv(M),'yo-')
hold on;plot([ranget(1) ranget(2) ranget(2) ranget(1) ranget(1)],[ranget(3) ranget(3) ranget(4) ranget(4) ranget(3)],'r-','linewidth',4)
end

%interpolate ASTERDEM to regular 15 m overlapping grid
znan=0;
tz1=double(data.z);tz1(tz1==znan)=NaN;data.z=tz1;
tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
tz(isnan(tz))=znan; %return toznan 
data1r= struct();
%data1r.x=tx;data1r.y=ty;  data1r.z=uint8(tz./max(tz(:))*255);
%data1r.x=tx;data1r.y=ty;  data1r.z=uint16(tz);

% custmorized scale %[min1 max1] to [0 max2] 
max2=1;
min1=min(tz(tz~=znan)); max1=max(tz(tz~=znan));
a=max2./(max1-min1);b=-a*min1;
tzs=tz*a+b;
tzs(tz==znan)=0;
data1r.x=tx;data1r.y=ty;  data1r.z=tzs;

%get the overlapping grid in UTM zone2
[X,Y]=meshgrid(tx,ty);
if 0 %slow
%utm 2 ell
clear PRO
for k=1:length(X(:))
PRO(k,1)={[num2str(zone),'V',num2str([X(k) Y(k)])]}; %6 longitude zone; N latitude zone, as long as north and south zones are distinguished.
end
% PRO={'32U460000.123 5498765.123'};
ELL=utm2ell(PRO,'wgs84',[]);
XYUTM1c=ell2utm(ELL,'wgs84',[],zone2); %faster using toolbox by Mike Craymer; 289 sec vs 0.8 sec
XYUTM1=zeros(length(XYUTM1c),2);
for i=1:length(XYUTM1c)
    letters=find(and(double(XYUTM1c{i})>=65,double(XYUTM1c{i})<=91));
    str1=strtrim(XYUTM1c{i}(letters(end)+1:end));
    cstr1=textscan(str1,'%f %f','Delimiter',' ');
    XYUTM1(i,:)=cell2mat(cstr1);
end
x1=XYUTM1(:,1);y1=XYUTM1(:,2);
else % use another toolbox
%[lat,lon]=utm2ell( 7084711.083,497795.565, zone2,wgsa,e2);
%etrue=[-147.044892194175 63.8895840826692]
% df=[ lon lat]*180/pi-etrue=1.0e-07 *[ 0.0003   -0.1771];
%df=[y1,x1]-[7084711.083,497795.565]=[ -0.0016   -0.0000];
tic; [lat,lon]=utm2ell( Y,X, zone,wgsa,e2);toc 
lcm=deg2rad(zone2*6-183);
tic; [y1,x1,zoneo]=ell2utm(lat,lon,wgsa,e2,lcm);toc;
end

if flagplot==1
figure (2);hold on;plot(x1(1:2000:end),y1(1:2000:end),'ro')
end
end  %if dataflag =1 or 2

%write 15 m resolution image
% interpolate image 2 on the overlapping grid.
tz1=double(data2.z);tz1(tz1==znan)=NaN;data2.z=tz1;
tz =interp2(data2.x,data2.y,double(data2.z),x1,y1,'*linear');%toc;%4s
tz(isnan(tz))=znan; %return toznan 
data2r= struct();
%no need to scale to uint8, uint16 works just fine.-> Use double.
%data2r.x=tx;data2r.y=ty;  data2r.z=uint8(tz./max(tz(:))*255);
%data2r.x=tx;data2r.y=ty;  data2r.z=uint16(tz);

% custmorized scale %[min1 max1] to [0 max2] 
min1=min(tz(tz~=znan)); max1=max(tz(tz~=znan));
a=max2./(max1-min1);b=-a*min1;
tzs=tz*a+b;
tzs(tz==znan)=0;
data2r.x=tx;data2r.y=ty;  data2r.z=tzs;

if flagplot==1
figure;imagesc(data1r.x,data1r.y,double(data1r.z));colorbar;axis equal
hold on;plot([ranget(1) ranget(2) ranget(2) ranget(1) ranget(1)],[ranget(3) ranget(3) ranget(4) ranget(4) ranget(3)],'r-','linewidth',4)
title('Overlapped Image 1 in Image 1 UTM coordinates')
figure;imagesc(data2r.x,data2r.y,double(data2r.z));colorbar;axis equal
hold on;plot([ranget(1) ranget(2) ranget(2) ranget(1) ranget(1)],[ranget(3) ranget(3) ranget(4) ranget(4) ranget(3)],'r-','linewidth',4)
title('Overlapped Image 2 in Image 1 UTM coordinates')
end

save test1.mat

projstr = 'UTM';
% OutName=['AST14DMO_00302222002211959_20180119102946_29249_15m.tif'];
OutName=strrep(infile1,'_V3N.tif','_15m.tif');
%writeGeotiff(OutName,data1r.x,data1r.y,data1r.z,1,0,projstr,zone)
%writeGeotiff(OutName,data1r.x,data1r.y,data1r.z,12,0,projstr,zone) %uint16
writeGeotiff(OutName,data1r.x,data1r.y,data1r.z,5,0,projstr,zone)

% OutName=['WV01_20160308_102001004D7D5100_102001004AAA7200_seg1_2m_ortho_15m.tif'];
OutName=['WV01_20160308_102001004D7D5100_102001004AAA7200_seg1_2m_ortho_15m.tif'];
if dataflag ==1
    OutName=strrep(infile2,'_2m_ortho.tif','_15m.tif');
elseif dataflag==2
OutName=strrep(infile2,'_V3N.tif','_15m.tif');
end
%writeGeotiff(OutName,data2r.x,data2r.y,data2r.z,1,0,projstr,zone)
writeGeotiff(OutName,data2r.x,data2r.y,data2r.z,5,0,projstr,zone)


elseif dataflag==3 % Both image are digital globe images
% ofmt: output format: 1, Digital Numbers, uint16; 2, reflectance, double
    ofmtflag=1; %default
    if exist('ofmt','var')
	if ofmt==2
	   ofmtflag=2;
	end
    end

    %get the polygon boundary
    metafile1=infile1;
    metafile2=infile2;
%    odir='./outcoreg/'
        if ~exist(odir,'dir')
          mkdir(odir)
        end
    
    [XYb1,range1]=imagebd(metafile1); %ref
    [XYb2,range2]=imagebd(metafile2);
    %get the image filenames
    if contains(name1,'meta.txt')
    flagfmt1=3; %stereo files
    % file names to be used
    infile1 = strrep(metafile1,'meta.txt','ortho.tif');
    OutName1= strrep(filename1,'meta','ortho_prep.tif'); %filename1: no directory, no extension
    elseif strcmp(ext1,'.xml')
    flagfmt1=1; %xml files mono
    infile1 = strrep(metafile1,'.xml','.tif');
    OutName1= [filename1,'_prep.tif'];
    end
	OutName1= [odir,OutName1];
    %2
    if contains(name2,'meta.txt')
    flagfmt2=3; %stereo files
    % file names to be used
    infile2 = strrep(metafile2,'meta.txt','ortho.tif');
        OutName2= strrep(filename2,'meta','ortho_prep.tif');
    elseif strcmp(ext2,'.xml')
    flagfmt2=1; %xml files mono
    infile2 = strrep(metafile2,'.xml','.tif');
        OutName2= [filename2,'_prep.tif'];
    end
	OutName2= [odir,OutName2];
    
    rangtar=[min(XYb2(:,1)) max(XYb2(:,1)) min(XYb2(:,2)) max(XYb2(:,2))];
    rangref=[min(XYb1(:,1)) max(XYb1(:,1)) min(XYb1(:,2)) max(XYb1(:,2))];   
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
%   reduce resolution
    ranget=[rangeov/resr];
    ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
    rangeov=ranget;
    tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
    
    if exist('rang0','var') 
    if ~isempty(rang0)
        rangeov=rang0;
    end
    end
    data1i=readGeotiff(infile1,'map_subset',rangeov);
    data2i=readGeotiff(infile2,'map_subset',rangeov);
    
    %Find image edges; Use the pixels that have data in both images
    edge1=data1i.z(:,:,1) == 0;
    edge2=data2i.z(:,:,1) == 0;
        
    %Digital Numbers to Reflectance; Notice: 0 to offset.
	if ofmtflag==2
    data1=DN2reflpan(data1i,metafile1);
    data2=DN2reflpan(data2i,metafile2);
	else
	data1=data1i;data2=data2i;
	end
    
    %Multiband to one band
    [~,~,nb]=size(data1.z);
    if nb ==1 %do nothing
        iNIR1=1; %panband
	iRed=1;
    elseif nb==4; iNIR1=4;
	iRed=3;
    elseif nb==8; iNIR1=7;
	iRed=5;
    end
    ref=double(data1.z(:,:,iRed)); %chose Near Infrared band
    [~,~,nb]=size(data2.z);
    if nb ==1 %do nothing
        iNIR1=1; %panband
	iRed=1;
    elseif nb==4; iNIR1=4;
	iRed=3;
    elseif nb==8; iNIR1=7;
	iRed=5;
    end
    tar=double(data2.z(:,:,iRed)); %chose Near Infrared band

    %interpolate to 15 m resolution grids
    ref(edge1)=nan;    tar(edge2)=nan;
    ref =interp2(data1.x,data1.y,double(ref),tx,ty','*linear',nan);
    tar =interp2(data2.x,data2.y,double(tar),tx,ty','*linear',nan);
    edge1 =interp2(data1.x,data1.y,double(edge1),tx,ty','*nearest',1); %1 is edge
    edge2 =interp2(data2.x,data2.y,double(edge2),tx,ty','*nearest',1);
    edge=edge1|edge2; %Image edges in either image will not be used.
%     ref(edge)=nan; tar(edge)=nan;
	if 0
 figure;imagesc(tx*1e-3,ty*1e-3,tar);colorbar;view(0,-90);title('tar');axis equal
 figure;imagesc(tx*1e-3,ty*1e-3,ref);colorbar;view(0,-90);title('ref');axis equal
	end

    %Apply the control surfaces map
    if 0 %(nargin == 3) 
    mp =interp2(groad.x,groad.y,double(groad.z),tx,ty','*nearest',0);
    M=mp&~edge;%map with valid data; 
    else
    M=~edge;%map with valid data;
    end
	if ofmtflag==2
        ref(~M)=nan;tar(~M)=nan;
    else
    ref(~M)=0;tar(~M)=0;
    end
    data.x=tx;data.y=ty;data.z=ref;
    [data1r]=cropmatrix(data,M);
%   data1r=data;
    data.x=tx;data.y=ty;data.z=tar;
    [data2r]=cropmatrix(data,M); %image too small for MIMC2
%   data2r=data;
	if 0
    figure;imagesc(data1r.x*1e-3,data1r.y*1e-3,data1r.z);colorbar;view(0,-90);title('ref');axis equal
    figure;imagesc(data2r.x*1e-3,data2r.y*1e-3,data2r.z);colorbar;view(0,-90);title('tar');axis equal
	end


	%apply rock mask; map out data at non-rock surface;
    if exist('tag','var') %1 rock; 0 non-rock;
	[nt,mt]=size(tag.z);
	if min([nt,mt])>2	   
	    
            rocktag = interp2(tag.x,tag.y,tag.z,data1r.x,data1r.y','*nearest',0); %logical
            rocktag = logical(rocktag);
	    data1r.z(rocktag~=1)=0;

            rocktag = interp2(tag.x,tag.y,tag.z,data2r.x,data2r.y','*nearest',0); %logical
            rocktag = logical(rocktag);
	    data2r.z(rocktag~=1)=0;
	end
    end
    M=data1r.z~=0&data2r.z~=0;
    [data1r]=cropmatrix(data1r,M);
    [data2r]=cropmatrix(data2r,M);

%    save testprepmj.mat -v7.3

    %Write data in double format
    projstr='polar stereo north';
    if isfield(data1i, 'Tinfo')
        data1r.Tinfo=data1i.Tinfo;
    else
        Tinfo       = imfinfo(infile1);%see readGeotiff.m
        data1r.Tinfo=Tinfo;%MIMC need time info!
    end
    if isfield(data2i, 'Tinfo')
        data2r.Tinfo=data2i.Tinfo;
    else
        Tinfo       = imfinfo(infile2);%see readGeotiff.m
        data2r.Tinfo=Tinfo;%MIMC need time info!
    end

    if ~(isempty(data1r.z) || isempty(data2r.z))
	if ofmtflag==2
    writeGeotiff(OutName1,data1r.x,data1r.y,data1r.z,5,nan,projstr)
    writeGeotiff(OutName2,data2r.x,data2r.y,data2r.z,5,nan,projstr)
	else
 	%uint16;setsm Input image should be unsigned integer 16 bits as original DG tif format.
        writeGeotiff(OutName1,data1r.x,data1r.y,uint16(data1r.z),12,0,projstr)
        writeGeotiff(OutName2,data2r.x,data2r.y,uint16(data2r.z),12,0,projstr)
	end
    else
        fprintf('\n prepareMJ.m error; one of the images is empty; Set pzxy to be zeros.')
    end
    
    %Plot figures; and create gif animations.
%   co=testgif(data1r,data2r);

elseif dataflag==4 % Both image are DEM files
    %get the polygon boundary
    metafile1=infile1;
    metafile2=infile2;
        if ~exist(odir,'dir')
          mkdir(odir)
        end
    
    %get the image filenames
    OutName1= strrep(name1,demext,'dem_prep.tif'); %filename1: no directory, no extension
    OutName1= [odir,OutName1];
    OutName2= strrep(name2,demext,'dem_prep.tif'); %filename1: no directory, no extension
    OutName2= [odir,OutName2];

    if 1 %old
      %if dem.tif file does not exist, try dem_10m.tif.
      if strcmp(demext,'dem.tif') && ~exist(infile1,'file')
         infile1b= strrep(infile1,demext,'dem_10m.tif');
         fprintf(['\n dem.tif does not exist, use dem_10m.tif instead for ',infile1,'\n.'])
      else
	 infile1b=infile1;
      end
      if strcmp(demext,'dem.tif') && ~exist(infile2,'file')
         infile2b= strrep(infile2,demext,'dem_10m.tif');
         fprintf(['\n dem.tif does not exist, use dem_10m.tif instead for ',infile2,'\n.'])
      else
	 infile2b=infile2;
      end

    data1i=readGeotiff(infile1b,'map_subset',rang0);
    data2i=readGeotiff(infile2b,'map_subset',rang0);
    
    for ii=1:2
	if ii==1;datar=data1i;metafile=infile1;
	elseif ii==2; datar=data2i;metafile=infile2;
   	end
    %apply geoid, same code as readdem.m
    filegeoid='geoidegm96grid15mpolar.tif';
    flaggeoid =0; %1 apply geoid undulation; 0 do not apply geoid;
    %SPOT 5 DEM, ALOS, SRTM, ASTER; to add (AST)
    if contains(metafile,'SPI')||contains(metafile,'GES')||contains(metafile,'ALPSML')||contains(metafile,'SRTM')
       flaggeoid=1;
       fprintf('\n Apply Geoid correction. \n')
       geoid=readGeotiff(filegeoid);
    end

    if flaggeoid==1 %appy geoid undulation to convert geoid height to ellipsoid height.
    geoid=readGeotiff(filegeoid);
    geoidz =interp2(geoid.x,geoid.y,double(geoid.z),datar.x,datar.y','*linear');% 65.5m ~ 67m
    datar.z(datar.z == -9999) = NaN;
    datar.z=datar.z+geoidz;
    datar.z(isnan(datar.z))=-9999; %return to -9999
    end
	if ii==1;data1i=datar;elseif ii==2; data2i=datar;end
    end %for ii

    else %; apply geoid for some dems; 
	%Bad: 1 not compatable with the following code of applying bitmask.
	% 2\ prefer to keep the option of not apply bitmask.	
        metafile1=strrep(infile1,demext,'meta.txt');
        if ~exist(metafile1,'file')
          metafile1=strrep(infile1,demext,'mdf.txt');
	end
        [XYbi,~]=imagebd(metafile1);
        metafile2=strrep(infile2,demext,'meta.txt');
        if ~exist(metafile2,'file')
          metafile2=strrep(infile2,demext,'mdf.txt');
	end
        [XYbi2,~]=imagebd(metafile2);
	resr=2; %Hope SETSM works better on 2m rather than 40m. -> No
        %  From MJ: "For SETSM coregistration, Coarse resolution is better.I recommend more than 10 meter resolution."
	%  Because the DEM slope and aspect is more smooth with coarser resolution.
	[data1i,nptsubrt]=readdem(XYbi,metafile1,rang0,resr);
	[data2i,nptsubrt]=readdem(XYbi2,metafile2,rang0,resr);
    end

    %apply filter masks %refers to readdem.m
    if flagfiltercoreg==1
    
    str1=demext; %'dem.tif';
    %get the mask
    maskfile1= strrep(infile1,str1,'bitmask.tif');
    if exist(maskfile1,'file') 
        if ~isempty(rang0)
            mask=readGeotiff(maskfile1,'map_subset',rang0); %0 is good; 1 is edge; > 1 bad data.
        else
            mask=readGeotiff(maskfile1); %0 is good; 1 is edge; > 1 bad data.
        end
 	% make sure the grids are the same as data
        M1= interp2(mask.x,mask.y,double(mask.z),data1i.x,data1i.y','*nearest',1); %logical
        mask.z= (M1); %compatible with older version

        data1i.z(mask.z>0)=-9999;
    end
    maskfile2= strrep(infile2,str1,'bitmask.tif');
    if exist(maskfile2,'file') 
        if ~isempty(rang0)
            mask=readGeotiff(maskfile2,'map_subset',rang0); %0 is good; 1 is edge; > 1 bad data.
        else
            mask=readGeotiff(maskfile2); %0 is good; 1 is edge; > 1 bad data.
        end
 	% make sure the grids are the same as data
        M1= interp2(mask.x,mask.y,double(mask.z),data2i.x,data2i.y','*nearest',1); %logical
        mask.z= (M1); %compatible with older version

        data2i.z(mask.z>0)=-9999;
    end

    end %if flagfiltercoreg==1

	%apply rock mask
    if exist('tag','var') %1 rock; 0 non-rock;
	[nt,mt]=size(tag.z);
	if min([nt,mt])>2	   
            rocktag = interp2(tag.x,tag.y,tag.z,data1i.x,data1i.y','*nearest',0); %logical
            rocktag = logical(rocktag);
	    data1i.z(rocktag~=1)=-9999;

            rocktag = interp2(tag.x,tag.y,tag.z,data2i.x,data2i.y','*nearest',0); %logical
            rocktag = logical(rocktag);
	    data2i.z(rocktag~=1)=-9999;
	end
    end

    %data1r=data1i;data2r=data2i;

    % crop on the overlapping area, with the same resolution (VERY IMPORTANT);
    resr1=round(abs(mean((diff(data1i.x)))));
    resr2=round(abs(mean((diff(data2i.x)))));
    resr=min([resr1,resr2]);
    resr=40; %for faster computation in writing files; From MJ: "For SETSM coregistration, Coarse resolution is better.I recommend more than 10 meter resolution."
    %overlapping boundary: the common area of target DEM, reference DEM
    rangref=[min(data1i.x) max(data1i.x) min(data1i.y) max(data1i.y)];
    rangtar=[min(data2i.x) max(data2i.x) min(data2i.y) max(data2i.y)];
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
    rangeov=round(rangeov/resr)*resr;
	
	ranget=rangeov;
    tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
    data=data1i;
    data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
    tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
    tz(isnan(tz))=-9999; %return to -9999
    refdemr= struct();
    refdemr.x=tx;refdemr.y=ty;  refdemr.z=tz;
    data=data2i;
    data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
    tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
    tz(isnan(tz))=-9999; %return to -9999
    tardemr= struct();
    tardemr.x=tx; tardemr.y=ty;  tardemr.z=tz;

    data1r=refdemr; data2r=tardemr;
    
    %Write data in double format
    projstr='polar stereo north';
    if isfield(data1i, 'Tinfo')
        data1r.Tinfo=data1i.Tinfo;
    else
        Tinfo       = imfinfo(infile1);%see readGeotiff.m
        data1r.Tinfo=Tinfo;%MIMC need time info!
    end
    if isfield(data2i, 'Tinfo')
        data2r.Tinfo=data2i.Tinfo;
    else
        Tinfo       = imfinfo(infile2);%see readGeotiff.m
        data2r.Tinfo=Tinfo;%MIMC need time info!
    end

    if ~(isempty(data1r.z) || isempty(data2r.z))
    writeGeotiff(OutName1,data1r.x,data1r.y,data1r.z,4,0,projstr)
    writeGeotiff(OutName2,data2r.x,data2r.y,data2r.z,4,0,projstr)
    else
        fprintf('\n prepareMJ.m error; one of the images is empty; Set pzxy to be zeros.')
    end

end %if dataflag

close all

return
end
