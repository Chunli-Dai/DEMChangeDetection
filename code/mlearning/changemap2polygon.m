function [shapefile]=changemap2polygon(pointfiles,imagefiles)
%given elevation change map, generate the slump outlines (saved as
%shapefiles)
constant
flagcon=1; % 1 use additonal contraint; 0 no constraint;
%In proposal, it sayes "disturbance areas commonly exceed 20 ha (200 000 m^2)"
% Ashley's polygons: minimum area 1e3 m^2, median area 11e3 m^2, mean area 24e3.
smlarea=1e3*2; %1000 m^2
flagpoint=1; %1, input point files; 0, input shapefiles; 2, no input point/shapefiles

if isempty(pointfiles)
   flagpoint=2; 
end

%[demdir,name,~] =fileparts([strtrim(pointfiles)]);
%shapefile=[demdir,'/',name,'_polygons.shp'];
shapefile=['change2polygon_tilei.shp'];

changefile=imagefiles{1};
data0=readGeotiff(changefile);
[ny,nx]=size(data0.z);
resx=mean(data0.x(2:end)-data0.x(1:end-1));resy=mean(data0.y(2:end)-data0.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);
Medgs1=false(size(data0.z));

%Get the mask from point shapefiles, 300 m buffer box around all points.
if flagpoint==1
points=load(pointfiles);ptlat=points(:,2);ptlon=points(:,1);%longitude, latitude
[ptx,pty]=polarstereo_fwd(ptlat,ptlon,[],[],70,-45);
ptxid=round((ptx-data0.x(1))/resx)+1;
ptyid=round((pty-data0.y(1))/resy)+1;
M=ptxid>=1&ptxid<=nx&ptyid>=1&ptyid<=ny;
ptxid(~M)=[];ptyid(~M)=[];
ptmask=false(size(data0.z));
for j=1:length(ptxid)
    ptmask(ptyid(j),ptxid(j))=1;
end

% 300m buffer of the points coordinates due to uncertainites of the field data!
widpix=round(50/resr); %50 m buffer
ptmask= imdilate(ptmask, ones(widpix*2)); 

elseif flagpoint==0 % polygon shapefiles of slumps from field data
ptmask=false(size(data0.z));

for k=1:length(pointfiles)
S1=shaperead(pointfiles{k});
n=length(S1);
for i=1:n
        if isfield(S1, 'zone')
          if contains(S1(i).zone,'debris');continue;end %only detect scars
        end

        [s1x,s1y]=polarstereo_fwd(S1(i).Y,S1(i).X,[],[],70,-45);
        s1x(isnan(s1x))=[];s1y(isnan(s1y))=[];
        Xb=round((s1x-data0.x(1))/resx)+1; 
        Yb=round((s1y-data0.y(1))/(resy))+1; %descending        

%	fprintf(['i, size of Xb Yb ny nx:',num2str([i length(Xb) length(Yb) ny nx])])
	%avoid out of memory in poly2mask
	if 0 %weird shape
	M=Xb>=1&Xb<=nx&Yb>=1&Yb<=ny;
	Xb(~M)=[];Yb(~M)=[];
	end
	maski=poly2mask(Xb,Yb, ny,nx); 

	ptmask=ptmask|maski;

end %for i

end % for k

%Do not apply buffer zone to ensure good quality.

elseif flagpoint==2 % no points or shapefiles

ptmask=true(size(data0.z));

end %if flagpoint

%End of ptmask

%Get the potential slump area
Mf4=data0.z<=-1; %abs(jump)>=2; %jump significant;Only negative
Mf4clean = bwareaopen(Mf4, smlarea/resr/resr); %remove small clusters

%Only keep the clusters that are marked as significant in bitmask.tif
%Or consider: keep the cluster that overlapps a point (300m buffer)from the pointfiles (from field measurements)
bitfile=strrep(changefile,'_jump.tif','_bitmask.tif');
bitmask=readGeotiff(bitfile);
T1name=strrep(changefile,'_jump.tif','_eventtimeT1.tif');
dataT1=readGeotiff(T1name);

Mf4cleankeep=false(size(data0.z));
BW=Mf4clean;
CC = bwconncomp(BW);
tic %slow: half an hour for 8000 out of 50,000 -> 3.5 hours for 56433.
yr =365.25; 
for k=1:CC.NumObjects
    BW3=false(size(BW));
    BW3(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});
    
    if flagcon==0 % no constraint
      overlappt=sum(sum(BW3.*ptmask));
    elseif flagcon==1  %to do : add water mask or glacier mask
      % This polygon must contain one pixel in the bitmask.
      overlappt=sum(sum(BW3.*logical(bitmask.z) ) ); 

      % get the std of T1 within the polygon.
      [nyT,nxT]=size(dataT1.z);
      if ny~=nyT || nx~=nxT
	fprintf(['\n Data size does not match, do interpolation: ',T1name])
        tz=interp2(dataT1.x,dataT1.y,double(dataT1.z),data0.x,data0.y','*nearest',0);
	dataT1.z=tz; dataT1.x=data0.x; dataT1.y=data0.y;
      end
      T1pt=double(dataT1.z(BW3&dataT1.z>2000e4)); %date, days
      if isempty(T1pt)
         T1std= NaN;
      else
      T1ptyr=datenum(num2str(T1pt),'yyyymmdd')/yr;
      T1std=std(T1ptyr);
      end
      % The std of T1 within this polygon must larger than 1 year.
      if T1std>1 %good 
      else %bad
	overlappt=0;
      end
	
      % remove water and glacier areas
	
    end
    
    if overlappt>0
        Mf4cleankeep=Mf4cleankeep|BW3;
    end

end
toc

if 0 %plot
figure;imagesc(data0.z,'alphadata',Mf4cleankeep);hold on;plot(ptxid,ptyid,'r.');title('M keep new')
colormap jet; colorbar;caxis([-3 3])
end %plot

%To do: Store Mf4cleankeep in changemap2polygon.m for each tile.
Modfil=~(Mf4cleankeep);  %Modfil (logical; 1 No change, 0 Change) 

%save t1.mat -v7.3

[Co]=mask2boundary(data0.x,data0.y,Modfil,Medgs1,shapefile);

% str1=['mv boundary.shp ',shapefile];

shapefile=Co; %added May 2023

return
end
