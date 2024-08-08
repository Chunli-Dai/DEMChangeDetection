function [datestro,flago]=strip2date(stripmetafile,flag)
%given strip file name, to get date YYYYMMDDHHMMSS
%input flag: 1 use the epoch of image 1; 2, use the epoch of image 2; 3, %use the mean epoch of two images
% output: flago: 1 use this data; 0 do not use this data
constant

flago=1; %default
[~,ifile,~] =fileparts(stripmetafile);

if ~exist(stripmetafile,'file')
  fprintf(['\n File does not exist:',stripmetafile])
  datestro='00000000000000';flago=0; return 
end

%c=textread(stripmetafile,'%s','delimiter','\n');
fid= fopen(stripmetafile);
c=textscan(fid,'%s','delimiter','\n');c=c{1,1};
fclose('all'); %fclose(fid);
% rs=find(~cellfun(@isempty,strfind(c,'scene')));
rs=find(contains(c,'scene')&contains(c,'dem'));

if ~isempty(rs) %strip DEM
%/fs/project/howat.4/EarthDEM/region_34_alaska_north/strips_unf/2m/W1W2_20110825_1020010014850E00_103001000C5BD600_2m_lsf/W1W2_20110825_1020010014850E00_103001000C5BD600_2m_lsf_seg1_meta.txt    
%or: /home/dai.56/chunliwork/landslide/site1Eureka/stripdata//SETSM_W1W2_20140801_1020010031B97700_10300100348BE800_seg3_2m_v3.0_mdf.txt
nsce=length(rs); %length(rs)-1; %number of scenes; 
for i=1 %nsce
% % c1=c{rs(1)+i};r1=strfind(c1,'dem.tif');c1(1:(r1+6))='';
% %new release scene name is '_dem_smooth.tif'
% c1=c{rs(1)+i};
% matchStr = regexp(c1,'.tif','split');c1=matchStr{end};
% tmp=sscanf(c1, '%g', 4);

c1=c{rs(0+i)+3};r1=[strfind(c1,'Image 1=');strfind(c1,'Image1')]; 
if(isempty(r1)) Warning('Image 1 not found') ; end
% r1=[strfind(c1,'/W');strfind(c1,'"W');]; %missing GE
str1='("|/)[G-W]';r1=regexp(c1,str1, 'start');
c1(1:(r1(end)))='';  

c11=deblank(c1);
satname=c11(1:4);%use this, since strip could be W1W2_20110825_1020010014850E00_103001000C5BD600_seg1_2m_meta.txt

%WV02_20160304214247_1030010052B75A00_16MAR04214247
datestro=c11(6:19);

%Image 2
c1=c{rs(0+i)+4};r1=[strfind(c1,'Image 2=');strfind(c1,'Image2')]; 
if(isempty(r1)) Warning('Image 2 not found') ; end
% r1=[strfind(c1,'/W');strfind(c1,'"W');];
str1='("|/)[G-W]';r1=regexp(c1,str1, 'start');
c1(1:(r1(end)))='';%missing GE

c11=deblank(c1);
satname2=c11(1:4);%use this, since strip could be W1W2_20110825_1020010014850E00_103001000C5BD600_seg1_2m_meta.txt

%WV02_20160304214247_1030010052B75A00_16MAR04214247
datestro2=c11(6:19);

end

else %rs is empty; scene DEM meta file
% both for  /fs/project/howat.4/EarthDEM/region_31_alaska_south/strips_v4/2m/W1W2_20120212_1020010018BC7700_1030010011A35E00_2m_lsf_v030208//W1W2_20120212_1020010018BC7700_1030010011A35E00_2m_lsf_seg1_meta.txt
% and for /fs/project/howat.4/dai.56/chunliwork/landslide/rtsyukon/newimages/setsmdir973b/setsmdir973b_meta.txt
%use 
%Image_1_satID=WV02
%Image_1_Acquisition_time=2012-02-12T22:21:21.357575Z

check={'Image_1_satID='};
Mcheck=contains(c,check{1});
c1all=c(Mcheck);
nsce=sum(Mcheck(:));

for i=1  %:nsce
c1=c1all{i};
r1=strfind(c1,'=');c1(1:(r1(end)))='';
c11=deblank(c1);
satname=c11(1:4);%
end

check={'Image_1_Acquisition_time='};
Mcheck=contains(c,check{1});
c1all=c(Mcheck);
c1=c1all{1};
r1=strfind(c1,'=');c1(1:(r1(end)))='';
c11=deblank(c1);
%2012-02-12T22:21:21.357575Z to 20120212222121
datestro=[c11(1:4),c11(6:7),c11(9:10),c11(12:13),c11(15:16),c11(18:19)];

%Image 2
check={'Image_2_satID='};
Mcheck=contains(c,check{1});
c1all=c(Mcheck);
nsce=sum(Mcheck(:));

for i=1  %:nsce
c1=c1all{i};
r1=strfind(c1,'=');c1(1:(r1(end)))='';
c11=deblank(c1);
satname2=c11(1:4);%
end

check={'Image_2_Acquisition_time='};
Mcheck=contains(c,check{1});
c1all=c(Mcheck);
c1=c1all{1};
r1=strfind(c1,'=');c1(1:(r1(end)))='';
c11=deblank(c1);
%2012-02-12T22:21:21.357575Z to 20120212222121
datestro2=[c11(1:4),c11(6:7),c11(9:10),c11(12:13),c11(15:16),c11(18:19)];

if 0 %rs is empty; scene DEM meta file %Old method 
% compatible for scene DEM meta file.
%/fs/byo/howat-data2/ArcticDEM/region_31_alaska_south/tif_results/2m/WV03_20161223_10400100278D9E00_10400100251E2C00_501077223020_01_P008_501075261100_01_P007_2_meta.txt
check={'Image 1=','/','.tif'};
Mcheck=contains(c,check{1})&contains(c,check{2})&contains(c,check{3});
c1all=c(Mcheck);
nsce=sum(Mcheck(:));

for i=1  %:nsce
c1=c1all{i};
r1=strfind(c1,'/');c1(1:(r1(end)))='';

c11=deblank(c1);
satname=c11(1:4);%use this, since strip could be W1W2_20110825_1020010014850E00_103001000C5BD600_seg1_2m_meta.txt
datestro=c11(6:19);
end

%Image 2
check={'Image 2=','/','.tif'};
Mcheck=contains(c,check{1})&contains(c,check{2})&contains(c,check{3});
c1all=c(Mcheck);
nsce=sum(Mcheck(:));

for i=1  %:nsce
c1=c1all{i};
r1=strfind(c1,'/');c1(1:(r1(end)))='';

c11=deblank(c1);
satname2=c11(1:4);%use this, since strip could be W1W2_20110825_1020010014850E00_103001000C5BD600_seg1_2m_meta.txt
datestro2=c11(6:19);
end
end % if 0 old method

end %if ~isempty(rs)

%Time apart between two images.
epoch1=datenum(datestro,'yyyymmddHHMMSS');
epoch2=datenum(datestro2,'yyyymmddHHMMSS');
dt=epoch2-epoch1; %days
%days to day hour minute seconds
if dt>=0; sign='+';else;sign='-';end; dt=abs(dt);
day=floor(dt);
hour=floor((dt-day)*24);
minute=floor( dt*24*60 - (day*24*60 + hour*60 ) );
second=( dt*24*60*60 - (day*24*60*60 + hour*60*60 + minute*60 ) );

if abs(dt)>maxdt 
flago=0;
fprintf(['\n Time span between Image 2 and Image 1 is: ',sign,' ', num2str(day),' day ', num2str(hour),' hour ', num2str(minute),' minute ', num2str(second),' seconds; ifile: ',ifile,'\n' ]);
  datestro='00000000000000'; return 
end

if exist('flag','var')
    if flag==1 %use the epoch of image 1
        %do nothing
    elseif flag==2 %use the epoch of image 2
        datestro=datestro2;
    elseif flag==3 %use the mean epoch of two images
        epoch1=datenum(datestro,'yyyymmddHHMMSS');
        epoch2=datenum(datestro2,'yyyymmddHHMMSS');
        c1=datestr(mean([epoch1, epoch2]),30 );
        r1=strfind(c1,'T');c1(r1)='';
        datestro=c1;
    end        
end

return
end
