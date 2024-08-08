% %%%% inputs needed
%see also /home/dai.56/chunliwork/ice/Quelccaya/prepare.m
% codetest/testotherdems.m 
macdir=[];

currentdir=pwd;
addpath(genpath([macdir,'/home/dai.56/arcticdemapp/landslide/code1/']));

%step 0: change file name use changenameASTER.sh;
%	 change coordinates use ortho2utm2m.sh

% %%% Preparation: create meta files;
filename='asterdemlist';
fprintf ('\n Step 0: geting the boundary for all files in the region.\n')
%READ INPUT PARAMETERS; getting the boundaries for all files
fid = fopen(filename);
n = linecount(fid);
fclose(fid);
fid11 = fopen(filename,'r');
range=zeros(n,4);XYbg=cell(n,1);
for i=1:n
   ifile=[macdir,fgetl(fid11)]; %name *dem.tif
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   if contains(demdir,'home') %absolute directory
        fdir{i}=[demdir,'/'];%working on two regions 
   else %relative directory
        fdir{i}=['./',demdir,'/'];%working on two regions 
   end
   satname=f{i}(1:4);
end

for i=1:n
    fprintf(['\n file ',num2str(i),': ','\n'])
   ifile=[fdir{i},f{i}] %name *dem.tif
   [demdir,name,ext] =fileparts([strtrim(ifile)]);

   % change the boundary of data to -9999
        %infile='20000211000000_SRTM_utm.tif';
        infile=ifile;
        projstr='polar stereo north';
        %OutName=strrep(infile,'_utm.tif','_utm_dem.tif');
        OutName=strrep(infile,'.tif','_dem.tif');
   if 1
        data=readGeotiff(infile);
        %SRTM NAN -32767; ASTER DEM edges: 0 and -9999
        %UAV DEM Edges:-10000
        data.z=double(data.z);
        clear M
        %M=data.z==0|data.z==-9999|data.z==32767|data.z==-32767|data.z<0;
        M=data.z==0|data.z==-9999|data.z==32767|data.z==-32767|data.z<-1000; %ellipsoidal height negative south of nepal
        data.z(M)=-9999;
        writeGeotiff(OutName,data.x,data.y,double(data.z),5,0,projstr)
   end

   % get the boundary from xml file
   ifile=OutName;
   %[XYbi,rangei]=imagebd(ifile);
   co=createmeta(ifile);
end

