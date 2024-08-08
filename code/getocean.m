function [tag]=getocean(rang0)
%Given rang0, get the ocean mask from GSHHS_f_L1.shp
% Input: rang0, 4 by 1 vector, e.g. [5738000     5891000     3857000     4020000]
% Output: tag.x tag.y tag.z (logical, 1 water, 0 land)

%refer to ~/arcticdemapp/coastline/codec2/CoastTileMonoMain.m CoastTileMono.m getrgibp2.m

resr=40;
resr=10;
exb=5e3;%expand the box by 5km for the given range.

%refers to /home/dai.56/arcticdemapp/river/rivergithub2/rivercore/getcl.m
rang0in=rang0;
rang0=[rang0in(1)-exb rang0in(2)+exb rang0in(3)-exb rang0in(4)+exb];

x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
[lat0,lon0]=polarstereo_inv(x0,y0,[],[],70,-45);

%[status , cmdout ]=system(['find ',codedir,' -name GSHHS_f_L1.shp']);
%shpname=deblank(cmdout);
shpname='/blue/chunlidai/apps/landslide/code1/GSHHS/GSHHS_f_L1.shp';

str1=[num2str(rang0(1)),'_',num2str(randi([0 1000]))]; %to avoid same filename in the case of parallel.
bb = geoshape(lat0,lon0,'Geometry','Polygon');
tileshape=['tile',str1,'.shp'];
%tilecoastname='tilegshhs.shp';
tilecoastname=['tilegshhs',str1,'.shp'];
shapewrite(bb,tileshape);
system(['rm ',tilecoastname])
system(['time ogr2ogr -overwrite -clipsrc ',tileshape,' ',tilecoastname,' ',shpname]);
%ogr2ogr -overwrite -clipsrc tile.shp tilegshhs.shp GSHHS/GSHHS_f_L1.shp
S2 = shaperead(tilecoastname);

% figure;mapshow(S2);
% figure;plot(S2(:).X,S2(:).Y,'k-');
ns=length(S2);

ranget=round(rang0/resr)*resr;rang0=ranget;
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
xout=tx;yout=ty;

wm.x=xout; wm.y=yout;
nwy=length(wm.y);nwx=length(wm.x);
Msv=false(nwy,nwx);%1 land, 0, ocean
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

        for k=1:length(C)
            idx=C{k}(:,1);idy=C{k}(:,2);   
            Mb=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
            Msv=Msv|Mb;%lost the holes 
        end        
        
end
    
tag.x=wm.x;tag.y=wm.y;tag.z=~Msv; %1 water; 0 land

if 0
save('watermask.mat','tag','-v7.3')

%% write output
projstr='polar stereo north';
OutName=['watermask.tif']; % %1 is good (landslide); 0 is bad
writeGeotiff(OutName,tag.x,tag.y,int32(tag.z),3,0,projstr)
end

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

return
end

