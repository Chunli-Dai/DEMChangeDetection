% main program for getting coastline for each ArcticDEM tile
% Requirements: gdal software
% %%%% inputs needed

constant 

macdir='/Users/chunlidai/surge/';
macdir=[];

currentdir=pwd;
%addpath(genpath(currentdir));
%addpath(genpath([macdir,'/data/chunli/landslide/']));
%addpath(genpath([macdir,'/data/chunli/scripts/']));
%addpath(genpath([macdir,'/home/dai.56/arcticdemapp/landslide/code1/']));
addpath(genpath([codedir]));

shpname='./GSHHS/GSHHS_f_L1.shp';% a priori coastline shapefile

%General directory that contains the tile DEM files, such as /elev/dem/setsm/ArcticDEM/mosaic/v2.0/
%tiledir='/Users/chunlidai/surge/data/chunli/coastline/';%ArcticDEM mosaic tile directory. 
%tiledir=[macdir,'/fs/byo/howat-data3/ArcticDEMmosaics/']; %
%stripdir='/*/ArcticDEM/region*/strips/2m/';
%stripdir='/fs/byo/howat-data2/ArcticDEM/region*/strips/2m/';
%stripdir='/fs/project/howat.4/EarthDEM/region*/strips_unf/2m/';
%stripdir=currentdir;

% %%%% control parameters
width=2e3; %buffer width of the a priori coastline, e.g., 2km.

if ~exist('mat0.mat','file') %readding boundary; time 1 hour for 90751 strip files and 10130 mono xml files

% %%% Preparation: get the list of strip files and boundries
filename='boundaries_regall_strip.dat'; %'boundaries_reg31.dat';
filename='striplist.dat';
if ~exist(filename,'file')
   %str=sprintf('find  %s -name ''*mdf.txt'' > %s',deblank(stripdir),filename);
   str=sprintf('find  %s -name ''*meta.txt'' > %s',deblank(stripdir),filename);
  [status, cmdout]=system(str);
end
fprintf ('\n Step 0: geting the boundary for all files in the region.\n')
%READ INPUT PARAMETERS; getting the boundaries for all files
% filename='boundaries_reg31.dat';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
range=zeros(n,4);XYbg=cell(n,1);
for i=1:n
   ifile=[macdir,fgetl(fid)];
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
   satname=f{i}(1:4);

   % get the boundary from xml file
   [XYbi,rangei]=imagebd(ifile);
   range(i,1:4)=rangei;XYbg{i}=XYbi;
end

save mat0.mat -v7.3

else 
load mat0.mat
end

dx=100e3;x0=-4000e3;y0=-4000e3;
dxs=dx/2/10; %5km
dxs=dx/50; %2km
%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;

% ArcticDEM mosaic tile grids
%Eureka lon lat: -85.9418555555556	79.988633
%32_33_2_2_5m_v2.0 ; 
%32_34_1_2_10_10 ; %work on tile, size of 5 km by 5 km.
%name convention; yid_xid_xids_yids_xidss_yidss
% inputtype=3;
%get the input from input.txt
filename='input.txt';
fid = fopen(filename);
inputtype=fscanf(fid, '%d', [1, 1])';
if inputtype ==1
    str=fgetl(fid);tilefile=fgetl(fid);
elseif inputtype==2 || inputtype==3
   latlon=fscanf(fid, '%f', [2, 1])';
   lateq=latlon(1);loneq=latlon(2);
elseif inputtype ==4
   rang0=fscanf(fid, '%f', [4, 1])';
end

% for xid=1:80
%     for yid=1:80
%         for xids=1:2
%             for yids=1:2
            % 54_06_2_2_5m_v2 yid_xid_xids_yids %name convention
            %     xid=6;yid=54;xids=2;yids=2;
            
            %     tilefile=sprintf('%02d_%02d_%01d_%01d_5m_v2.0_reg_dem.tif',yid,xid,xids,yids) ; %'54_06_2_2_5m_v2.0_reg_dem.tif';
             
            
            switch inputtype
              case 1 %based on input xid etc.
%                   xid=33;yid=32;xids=2;yids=2;xidss=1;yidss=1;
%                   tilefile=sprintf('%02d_%02d_%01d_%01d_%02d_%02d.tif',yid,xid,xids,yids,xidss,yidss);
                  tilefile=deblank(tilefile);
              case 2 %find the block based on input coordinates;
%                   loneq=-85.9418555555556; lateq=79.988633; %Eureka
%                   lateq=61.143; loneq= -148.16;
                  [xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
                  %x=x0+(xid-1)*dx+(xids-1)*dx/2+(xidss-1)*dxs;y=y0+(yid-1)*dx+(yids-1)*dx/2+(yidss-1)*dxs;
                  xid=floor((xeq-x0)/dx)+1; yid=floor((yeq-y0)/dx)+1;
                  xids=floor((xeq-(x0+(xid-1)*dx))/(dx/2))+1;
                  yids=floor((yeq-(y0+(yid-1)*dx))/(dx/2))+1;
                  xidss=floor((xeq-(x0+(xid-1)*dx+(xids-1)*dx/2))/dxs)+1;
                  yidss=floor((yeq-(y0+(yid-1)*dx+(yids-1)*dx/2))/dxs)+1;
                  tilefile=sprintf('%02d_%02d_%01d_%01d_%02d_%02d.tif',yid,xid,xids,yids,xidss,yidss) ;
              case 3 %Use a rectangle box around the input coordinates;
                % get the data boundary, rang0, of this DEM tile 
%                   loneq=-151.1055; lateq=59.399; % 59.399, -151.1055; Kinnikinnick Landslide, which happened between Aug 23 and Sept 4 2017
                  [xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
                  exb=1e3; %2.5e3;
		  rang0=[xeq-exb xeq+exb yeq-exb yeq+exb ]; %for test only
                % rang0=[-711.4e3-500 -710.4e3+500 -821.9e3-500 -820.9e3+500]; %2km by 2km
		% rang0=[    -1087000     -995000    -2004000    -1849000]; % entire Barnes;
		% rang0=[267000	320000	-2597000	-2545000]; %Helheim

                  tilefile=rang0;
                % [Co]=Landslide(rang0,range,XYbg,f,fdir);
	      case 4 %use boundary
                  
                 % xeq=(rang0(1)+rang0(2))/2;yeq=(rang0(3)+rang0(4))/2;
		 % [lateq,loneq]=polarstereo_inv(xeq,yeq,[], [],70,-45);
                  tilefile=rang0;

              otherwise
                  disp('Unknown input zone! Try again!') 
            end    
            

            tic
            [Co]=Landslide(tilefile,range,XYbg,f,fdir);
            fprintf(['\n Tile ',tilefile])
            toc

	    if flagoutput==0 %remove all files except output
%           [status, cmdout]=system('rm -f [h-Z]*'); %save dX4Sg.mat %dangerous! possibly remove all files in other folder!
%          [status, cmdout]=system('rm -f [c-C]*');
%          [status, cmdout]=system('rm -f e*');
	   [status, cmdout]=system('rm -f striplist.dat mat0.mat Tilemain.m constant.m input.txt test3.mat outliertp[0-9].mat tp[0-9].fig outlierex1.mat sv1.mat jumpstd.pdf jumpstd.fig jump.pdf jump.fig eventtime.pdf eventtime.fig eventtimestd.pdf eventtimestd.fig eventtime_start.fig eventtime_end.fig OverlappingCount_Final.fig Chosenisel.fig jumpfiltered.pdf jumpfiltered.fig jumpfiltered.dat core.[0-9]*[0-9]'); 
           [status, cmdout]=system('rm -rf outcoregdem');
	    end % flagoutput
            
            exit
%             end %yids
%         end%xids
%     end % yid
% end %xid

