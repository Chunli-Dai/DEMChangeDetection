
constant

%addpath(genpath('/Users/chunlidai/UFwork/ESI2019/landslide/'));
%addpath(genpath('/Users/chunlidai/UFwork/ESI2019/mlearning/'));
addpath(genpath('/home/chunlidai/blue/apps/landslide/code1/'));  %Directory of codes.

flagpre=0; %1, calculate precision; 0 don't calculate precision, just get the combined polygons.
flagcri1to3=0; %1 apply criteria 1to3; 0 do not apply.
flagcom=1; % 1, combine two ML results; 0, just use one ML results.

classtype_g = {'_vvl', '_vl', '_l', '_m', '_s', '_vs'};
[demdir,name,ext] =fileparts([strtrim(regiondir)]);
regid=name(11:12);
regionstr=['region',regid]; %'region09';

ncomall=0;
count_g=zeros(6,1);
for ic=1:6
%for ic=[3]
   clear data0  ids  ns  Mpcg

   classtype=classtype_g{ic};

% % load ./ixy_peel_cand.mat
%load ixy_region10_cand.mat
%maskdir='slump_case5band3_cand/';
%oshapefile='change2polyg_10_ai.shp'; %after Machine Learning

matfile=['ixy_',regionstr,'_cand',classtype,'.mat'] ;
load(matfile);
maskdir=['slumpsub_',regionstr,classtype,'_mask/']    %'slumpsub_region04_vvl_mask/';
%oshapefile=['change2polyg_',regid,classtype,'_ai.shp']; %after Machine Learning
oshapefile=['change2polyg_',regid,classtype,'_ai_com.shp']; %after Machine Learning

gtshape='~/blue/chunlidai/landslide/mlearning/region10onet/groundtruth_perfect/maskgroundtruth.shp';

%reconstruct data0, 10m; see hpg: /home/chunlidai/blue/chunlidai/landslide/mlearning/Dividemap_candpar.m
data0.z=false(length(data0.y),length(data0.x));
for ix=1:ns
        Mpc=Mpcg{ix};
        Mx=Mpc.Mx;My=Mpc.My;
        
        Mpcpre.z=data0.z(My,Mx); 

        %data0.z(My,Mx)=Mpc.z|Mpcpre.z;
        tmp=Mpc.z|Mpcpre.z;
        data0.z(My,Mx)=tmp;
end

% maskdir='slump_case2timeband3_thre70_cand/';
% maskdir='slump_case5/';'slump_case5band3_cand';
% maskdir='/Users/chunlidai/Library/CloudStorage/OneDrive-UniversityofFlorida/print/peel/slump_case5band3_cand/';
%maskdir='slump_case5band3_cand/';
%oshapefile='change2polyg_10peel_ai.shp'; %after Machine Learning

% reconstruct mask from machine learning 1.
flagold=0;[mask1,Mpcg_ml]=getmaskml(data0,Mpcg,ns,ids,maskdir,flagold);

if flagcom==1
% reconstruct mask from machine learning 2.

fprintf('\n Combine ML results.\n')

maskdir2=['case5b/',maskdir] 
flagold=0;[mask2i,Mpcg_ml2]=getmaskml(data0,Mpcg,ns,ids,maskdir2,flagold);
mask2.z=mask2i.z; %mask2.z&mask2i.z;
maskdir2=['case5aug/',maskdir] 
flagold=0;[mask2i,Mpcg_ml2]=getmaskml(data0,Mpcg,ns,ids,maskdir2,flagold);
mask2.z=mask2.z&mask2i.z;
maskdir2=['case7/',maskdir] 
flagold=0;[mask2i,Mpcg_ml2]=getmaskml(data0,Mpcg,ns,ids,maskdir2,flagold);
mask2.z=mask2.z&mask2i.z;
maskdir2=['case7aug/',maskdir] 
flagold=0;[mask2i,Mpcg_ml2]=getmaskml(data0,Mpcg,ns,ids,maskdir2,flagold);
mask2.z=mask2.z&mask2i.z;

%Combine mask 1 and mask2
mask1.z=mask1.z&mask2.z;
%mask1.z=mask2.z;
end


[mask1crop]=cropmatrix(mask1,mask1.z);
try
    ofilefig=['region',regid,classtype,'_ml'];

    if 0
    figure;imagesc(mask1crop.x*1e-3,mask1crop.y*1e-3,mask1crop.z);colorbar;title(ofilefig)
    saveas(gcf,ofilefig,'jpg')
    saveas(gcf,ofilefig,'fig')
    end
   
projstr='polar stereo north';
OutName=['region',regid,classtype,'_ml.tif'];
writeGeotiff(OutName,mask1crop.x,mask1crop.y,int8(mask1crop.z),11,0,projstr)
 
catch e
   fprintf('There was an error! The message was:\n%s',e.message);
end

if flagpre==1
%S2=shaperead('maskgroundtruth.shp');
S2=shaperead(gtshape);

fprintf('\n Precision for ML only:');
tic; [precision,recall]=poly2precision(mask1crop,S2);toc
%tic; [precision,recall]=poly2precision(Mpcg_ml,S2);toc  %Not going to work, because cases of multiple polygons for one images, causing the duplicates count of slumps.
end

%Get the candidates mask from Change2poly6classespar.m dividemap_cand.m
% Mf4clean = data0.z; 
% Get precision and recall
% S2=shaperead('/Users/chunlidai/UFwork/ESI2019/mlearning/data/maskgroundtruth.shp');% figure;mapshow(S1)
% S2=shaperead('/Users/chunlidai/UFwork/ESI2019old/machinelearning/groundtruth_perfect/maskgroundtruth.shp');
% S2=shaperead('/Users/chunlidai/Downloads/maskgroundtruth.shp');
% S2=shaperead('/Users/chunlidai/Library/CloudStorage/OneDrive-UniversityofFlorida/print/groundtruth_perfect/maskgroundtruth.shp');
%S2=shaperead('maskgroundtruth.shp');

[wm]=cropmatrix(data0,data0.z);

if flagpre==1
fprintf('\n Precision for data0.z:')
tic; [precision,recall]=poly2precision(wm,S2);toc
%tic; [precision,recall]=poly2precision(Mpcg,S2);toc
end

% tic; [precision,recall]=poly2precision('change2polyg_10.shp',S2);toc
try
    ofilefig=['region',regid,classtype,'_con'];
    figure;imagesc(wm.x*1e-3,wm.y*1e-3,wm.z);colorbar;title(ofilefig)
    saveas(gcf,ofilefig,'jpg')
    saveas(gcf,ofilefig,'fig')
projstr='polar stereo north';
OutName=['region',regid,classtype,'_con.tif'];
writeGeotiff(OutName,wm.x,wm.y,int8(wm.z),11,0,projstr)
catch e
   fprintf('There was an error! The message was:\n%s',e.message);
end

% Select candidates that have >70% areas classified as RTS by ML.
%maskconstr=interp2(data0.x,data0.y,double(Mf4clean),mask1.x,mask1.y','nearest',0);
%Merged candidates; -> this will remove the edge issues around 50 km tiles
%Poly2mask in dividemapsub.m may not fully recover the original mask, and unintentionally separate clusters.
% hence artificially increase the number of total test clusters (and Negative positives).
Mf4cleankeep=false(size(data0.z));
%BW=maskconstr;
%CC = bwconncomp(BW); %takes large memory, and slow.
fprintf(['\n The total number of candidates in data0.z: ',num2str(ns),'.']);
tic 
%for k=1:0 %CC.NumObjects %slow
%    BW3=false(size(BW));
%    BW3(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});
%    
%    ratio=sum(sum(BW3.*mask1.z))./sum(sum(BW3))*100;
%    
%    if ratio>70
%        Mf4cleankeep=Mf4cleankeep|BW3;
%    end
%end
 
ratiog=zeros(ns,1);
for k=1:ns
    fprintf(['\n Working on k: ',num2str(k),' / ',num2str(ns),'.\n'])

        Mpc=Mpcg{k};
        Mx=Mpc.Mx;My=Mpc.My;
%	Mpc.x=data0.x(Mx); % uncropped data0
 %       Mpc.y=data0.y(My);%
	% Mpc.z;
	
	BW3=Mpc.z;
	%mask1k.z=interp2(mask1crop.x,mask1crop.y,double(mask1crop.z),Mpc.x,Mpc.y','nearest',0); %slow
	rang2=[min(Mpc.x) max(Mpc.x) min(Mpc.y) max(Mpc.y)]; 
        Mx1=data0.x>=rang2(1)&data0.x<=rang2(2);My1=data0.y>=rang2(3)&data0.y<=rang2(4);
	mask1k.z=mask1.z(My1,Mx1);
	
	ratio=sum(sum(BW3&mask1k.z))./sum(sum(BW3))*100;
	ratiog(k)=ratio;

%   	if ratio>50
%	  BW3=false(size(data0.z));
%	  BW3(My,Mx)=Mpc.z;
%	  Mf4cleankeep=Mf4cleankeep|BW3;
%	end
end
toc

save ratiog.mat ratiog

fprintf('\n Time to collect the combined (ML+con) clusters:')
tic
thres=50;
idk=find(ratiog>thres);

if 1
% % quick way: load original shapefile refers to Dividemap_candpar_vvl.m
if strcmp(classtype,'_vvl')
   classtypein='_landslide';
else
   classtypein=classtype;
end
shapefile=['change2polyg_',regid,classtypein,'.shp';] 
S1=shaperead(shapefile);
shpc=struct2cell(S1)';
areaj=[shpc{:,6}]';
M1=areaj>=2000;
S1s=S1(M1);S1=S1s;
if length(S1)~=ns
   fprintf(['\n ',classtype, ': loaded shapefile does not match the ns in the mat file: ',matfile]);
   exit
end

Co=S1(ratiog>thres);

if ~isempty(Co)
 shapewrite(Co, oshapefile); % May 2023
end

fprintf(['Combined with ML, there are a total of ',num2str(length(Co)), '  ',classtype,' clusters.\n'])
count_g(ic)=length(Co);

if flagcri1to3==1
%see /home/chunlidai/blue/apps/landslide/code1/mlearning/Getcri1to3.m
oshapefile2=['change2polyg_',regid,classtype,'_ai_cri1to3.shp']; %after Machine Learning
crimatfile=['ixy_',regionstr,'_cand',classtype,'_cri1to3.mat'];
load(crimatfile); 
Co=S1(ratiog>thres & Mflag==1);
shapewrite(Co, oshapefile2);
fprintf(['Combined with ML and Criteria 1to3, there are a total of ',num2str(length(Co)), '  ',classtype,' clusters.\n'])
end


if flagpre==1
fprintf('\n Precision for combined:')
tic; [precision,recall]=poly2precision(Co,S2);toc
%tic; [precision,recall]=poly2precision(Mpcg_mlcom,S2);toc
end

% % %%%%%%%%%
else  % % Slower 
%M1good=ratiog>thres;
%Mpcg_mlcom=Mpcg;
%Mpcg_mlcom(~M1good)=[]; %to check
for i=1:length(idk)
    fprintf(['\n Working on i: ',num2str(i),' / ',num2str(length(idk)),'.\n'])
k=idk(i);
Mpc=Mpcg{k};
Mx=Mpc.Mx;My=Mpc.My;
%BW3=false(size(data0.z));
%BW3(My,Mx)=Mpc.z;
%Mf4cleankeep=Mf4cleankeep|BW3; % 5 hours; too slow
        Mpcpre.z=Mf4cleankeep(My,Mx);
        tmp=Mpc.z|Mpcpre.z;
        Mf4cleankeep(My,Mx)=tmp; % 5 hours reduced to 0.3 seconds.
end %for i
toc

mask1c=mask1; %get x y
mask1c.z=Mf4cleankeep;

[mask1cc]=cropmatrix(mask1c,mask1c.z);

if flagpre==1
fprintf('\n Precision for combined:')
tic; [precision,recall]=poly2precision(mask1cc,S2);toc
%tic; [precision,recall]=poly2precision(Mpcg_mlcom,S2);toc
end

Mf4cleankeep=mask1cc.z;

%save the selected RTS after ML.
Modfil=~(Mf4cleankeep);  %Modfil (logical; 1 No change, 0 Change) 
Medgs1=false(size(Mf4cleankeep));
[Co2]=mask2boundary(mask1cc.x,mask1cc.y,Modfil,Medgs1,oshapefile);
end  % if 1 ; fast or slow

ncomall=ncomall+length(Co);
%fprintf(['Combined with ML, there are a total of ',num2str(length(Co)), '  ',classtype,' clusters.\n'])

end %ic

fprintf(['Combined with ML, there are a total of ',num2str(sum(count_g(1:4))), ' for _vvl _vl _l _m types of clusters.\n'])
fprintf(['Combined with ML, there are a total of ',num2str(sum(count_g(5:6))), ' for _s _vs types of clusters.\n'])
fprintf(['Combined with ML, there are a total of ',num2str(ncomall), ' for all types of clusters.\n'])

return
