function [tag]=getrgi(rang0)
%Given rang0, get the off-ice mask from the Randolph Glacier Inventory (RGI 6.0).
% Input: rang0, 4 by 1 vector, e.g. [5738000     5891000     3857000     4020000]
% Output: tag.x tag.y tag.z (logical, 1 rock, 0 non rock)

%refer to /fs/project/howat.4/dai.56/chunliScripts/scripts/plotsrtmvs.m

resr=40;
resr=10;
exb=5e3;%expand the box by 5km for the given range.

%refers to /home/dai.56/arcticdemapp/river/rivergithub2/rivercore/getcl.m
rang0in=rang0;
rang0=[rang0in(1)-exb rang0in(2)+exb rang0in(3)-exb rang0in(4)+exb];

x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
[lat,lon]=polarstereo_inv(x0,y0,[],[],70,-45);
bbox=[min(lon) min(lat); max(lon) max(lat)];

range=[bbox(1,1) bbox(2,1) bbox(1,2) bbox(2,2)];
x0=[range(1) range(2) range(2) range(1) range(1) ];y0=[range(4) range(4) range(3) range(3) range(4) ];
% figure;hold all;plot(x0,y0,'-','linewidth',3)

%glacier polygon shapefile
%infile='glims_2019/glims_polygons_2D.shp';
infile='/blue/chunlidai/apps/horizontalv/code/glims_2023/RGI2000v7.0C0102_2D.shp'; 
% infile='glims_download_53175/glims_polygons_3.shp';
if 1 %old
S2=shaperead(infile,'BoundingBox',bbox);%10 minutes; tends to lose holes in shapefiles.
else %use ogr2ogr to subset, which keeps holes. ->mask still lose the holes.
%ogr2ogr -f "ESRI Shapefile" glims_polygons_himalayas.shp glims_2019/glims_polygons_2D.shp  -spat 79 30.6 79.6 31.1
%regionR=['-R',num2str(min(lon)),'/',num2str(max(lon)),'/',num2str(min(lat)),'/',num2str(max(lat))];
%-spat xmin ymin xmax ymax
regionR=num2str([min(lon) min(lat) max(lon) max(lat)]);
str=['ogr2ogr -f "ESRI Shapefile" glims_polygons_sub.shp ',infile,'  -spat ',regionR] 
[status, cmdout]=system(str);
S2=shaperead('glims_polygons_sub.shp');%
end
% figure;mapshow(S2);
% figure;plot(S2(:).X,S2(:).Y,'k-');
ns=length(S2);
%Bug: Error using openShapeFiles>readHeaderTypeCode(line 145)
%Unsupported shape type PolygonZ (type code = 15).
%Solution: in qgis, convert 3D shapefile to 2D.
%         step: saveas ESRI shapefile, select Geometric Type -> Polygon
%               uncheck "include z-dimension "


%get the off-ice mask
%refers to /home/dai.56/chunliwork/ice/Helheim1/prep_rangerock.m
ranget=round(rang0/resr)*resr;rang0=ranget;
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
xout=tx;yout=ty;

wm.x=xout; wm.y=yout;
nwy=length(wm.y);nwx=length(wm.x);
% Mcb=true(nwy,nwx);idblock=[];
Msv=false(nwy,nwx);%1 ice, 0, none ice
% Msv=zeros(nwy,nwx,'int8');%1 ice, 0, none ice
%refers to /home/dai.56/arcticdemapp/coastline/codec2/CoastTileMono.m
for j=1:ns
        XYbi=[S2(j).X(:),S2(j).Y(:)];
        lon1=XYbi(:,1);lat1=XYbi(:,2);
%         M1=isnan(lon1)|isnan(lat1); lon1(M1)=[];lat1(M1)=[];
        [Xb,Yb]=polarstereo_fwd(lat1,lon1,[],[],70,-45);
%          hold on;plot(lon1,lat1,'.-')
        
        % check whether this polygon intersect with the block
        idx=round((Xb-wm.x(1))/resr)+1;
        idy=round((Yb-wm.y(1))/(-resr))+1;
        
        %Oct 10, 2018:fix bug 14; separate polygons using the NaN;
        M=[idx(:),idy(:)];
        idx = any(isnan(M),2);
        idy = 1+cumsum(idx);
        idz = 1:size(M,1);
        C = accumarray(idy(~idx),idz(~idx),[],@(r){M(r,:)});
        if 1 %old way; holes are gone
        for k=1:length(C)
            idx=C{k}(:,1);idy=C{k}(:,2);   
            Mb=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
            %Msv=Msv|Mb;%lost the holes 
%             Msv=int8(Mb)-Msv; %keep the holes -> bad in Yahtse: -3336e3 -3184e3 334e3 400e3 
	    Msv=logical(abs(int8(Mb)-int8(Msv))); % March 2023
        end        
        else %try to save the holes in the mask. -> failed
%         XY={[idx(:),idy(:)]};
%         Mb = mpoly2mask(XY, [nwy,nwx]);
          Mb = mpoly2mask(C, [nwy,nwx]);
            Msv=Msv|Mb;
        end
        
%         figure;surf(lon,lat,double(Msv));shading interp;colorbar;colormap jet;view(0,90)
%          hold on;plot(lon1,lat1,'.-')
%         pause
%         
%         Mb = poly2mask(idx,idy, nwy,nwx); % build polygon mask       
%         overl=Mb&Mcb; %1, ice 
        %it's possible to work on one pixel block.
%         if(sum(sum(overl))>0);idblock=[idblock;j]; %idblock: the index of idregion, dX4Sg
%             fprintf(['\n Polygon within area of interest : ',num2str(j)])
%             hold all;plot(Xb*1e-3,Yb*1e-3,'.-');title([num2str(j)])
%             Msv=Msv|Mb;
%             pause
%         end
end
    
% Msv=abs(Msv)==1;
tag.x=wm.x;tag.y=wm.y;tag.z=~Msv; %1 rock; 0 non rock

%save('BarnesrockRGI.mat','tag','-v7.3')
%save('rockmask.mat','tag','-v7.3')

%% write output
projstr='polar stereo north';
OutName=['rockmask.tif']; % %1 is good (landslide); 0 is bad
writeGeotiff(OutName,tag.x,tag.y,int32(tag.z),3,0,projstr)

if 0 %plot
    
[X,Y]=meshgrid(tag.x,tag.y);
[lat,lon]=polarstereo_inv(X,Y,[], [],70,-45);
figure;surf(lon,lat,double(tag.z));shading interp;colorbar;colormap jet;view(0,90)
for j=1:ns
    hold all;plot(S2(j).X(:),S2(j).Y(:),'-')
end
figure;hold all;mapshow(S2);

hold on;plot(lon1,lat1,'.-')
figure;surf(lon,lat,double(Msv));shading interp;colorbar;colormap jet;view(0,90)
end

% figure; hold all;imagesc(tag.x*1e-3,tag.y*1e-3,tag.z);colorbar; colormap jet
%save testrgi.mat -v7.3

return
end

