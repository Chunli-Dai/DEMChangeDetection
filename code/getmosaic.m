

function [data0]=getmosaic(rang0,flago,pubres)
%Modified from coreg2.m 
% Given a region, to merge one DEM from DEM mosaics from PGC website.
% flago: 0, do not merge DEMs; 1 merge DEMs (default);
% pubres: the published resolution of the products available at PGC. default is 2.

constant

%addpath(genpath(['/home/dai.56/arcticdemapp/landslide/code1/']));

if nargin==1; flago=1;end
if ~exist('pubres','var');pubres=2;end

data0=[];

%tiledir=['/home/dai.56/data3/ArcticDEMmosaics/'];
tiledirnew=tiledir  %to store new downloaded dems

flagplot=0;

% find the DEM mosaics within rang0.
xr0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];yr0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];

resrc=2 ; %40;
resr=2; %10; %2;
%resr=40;
if flago==1
ranget=round(rang0/resrc)*resrc;rang0=ranget;
tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);
xoutr=tx;youtr=ty;
data0r.x=xoutr;data0r.y=youtr;data0r.z=nan(length(youtr),length(xoutr));
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
xout=tx;yout=ty;
data0.x=xout;data0.y=yout;data0.z=nan(length(yout),length(xout));
end 

fprintf ('\n Getting the tile DEM. \n')
% get the data boundary, rang0, of this DEM tile 
dx=100e3;x0=-4000e3;y0=-4000e3;%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;
% x=x0+(xid-1)*dx+(xids-1)*dx/2;y=y0+(yid-1)*dx+(yids-1)*dx/2; %bottom left
% rang0b=[x x+dx/2 y y+dx/2]; %exact tile boundary;
icount=0;clear xc yc
for i=1:4 %[1,3]
    icount=icount+1;
    x=xr0(i);y=yr0(i);
    xid=floor((x-x0)/dx)+1;
    yid=floor((y-y0)/dx)+1;
    xids=floor((x-x0-(xid-1)*dx)/(dx/2))+1;
    yids=floor((y-y0-(yid-1)*dx)/(dx/2))+1;
    xc(icount)=x0+(xid-1)*dx+(xids-1)*dx/2; yc(icount)=y0+(yid-1)*dx+(yids-1)*dx/2; %bottom left
%bug chunli's used (before august 2023) yid_xid_xids_yids not consistent with PGC's convention, yid_xid_yids_xids
%     tilefile=sprintf('%02d_%02d_%01d_%01d_5m_v2.0_reg_dem.tif',yid,xid,xids,yids)
%     x=xc(icount);y=yc(icount);
%     rang0b=[x x+dx/2 y y+dx/2];
%     xb0=[rang0b(1) rang0b(2) rang0b(2) rang0b(1) rang0b(1) ];
%     yb0=[rang0b(4) rang0b(4) rang0b(3) rang0b(3) rang0b(4) ];
%     hold on;plot(xb0*1e-3,yb0*1e-3,'.-')
end
icount=0;
for x=min(xc):dx/2:max(xc)
    for y=min(yc):dx/2:max(yc)
        icount=icount+1;
    xid=floor((x-x0)/dx)+1;
    yid=floor((y-y0)/dx)+1;
    xids=floor((x-x0-(xid-1)*dx)/(dx/2))+1;
    yids=floor((y-y0-(yid-1)*dx)/(dx/2))+1;
    %tilefile=sprintf('%02d_%02d_%01d_%01d_2m_v3.0_reg_dem.tif',yid,xid,xids,yids);
    %tilefile=sprintf('%02d_%02d_%01d_%01d_2m_v4.1_dem.tif',yid,xid,xids,yids);
    if pubres==2
	    %17_27_2_1_2m_v4.1_dem.tif 
        tilefile=sprintf('%02d_%02d_%01d_%01d_2m_v4.1_dem.tif',yid,xid,yids,xids);
    elseif pubres==32
        tilefile=sprintf('%02d_%02d_32m_v4.1_dem.tif',yid,xid); %17_27_32m_v4.1_dem.tif 
    end
    fprintf(['\n Get tile ',tilefile,'. \n'])
    %check
%     xp=x0+(xid-1)*dx+(xids-1)*dx/2;yp=y0+(yid-1)*dx+(yids-1)*dx/2; %bottom left
%     rang0b=[xp xp+dx/2 yp yp+dx/2];
%     xb0=[rang0b(1) rang0b(2) rang0b(2) rang0b(1) rang0b(1) ];
%     yb0=[rang0b(4) rang0b(4) rang0b(3) rang0b(3) rang0b(4) ];
%     hold on;plot(xb0*1e-3,yb0*1e-3,'o-')
%     pause
    
        % find the dem tile file or download the data
    [status , cmdout ]=system(['find ',tiledir,' -name ',tilefile]); %status always 0, cmdout could be empty.
    if ~isempty(cmdout) && status ==0 % 
        tilefile=deblank(cmdout);
    else
        warning(['Tile file ',tilefile,' not found! Download it from website'])
        [dir,name,ext] =fileparts(tilefile);
    %http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/43_59/43_59_2_1_5m_v2.0.tar
    %2m: https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/latest/2m/17_27/17_27_2_1_2m_v4.1.tar.gz
    %32m: https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/latest/32m/17_27/17_27_32m_v4.1.tar.gz
    tarfile=[name(1:end-4),'.tar'];
%       webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/',name(1:5),'/',tarfile,'.gz']);
        %webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/2m/',name(1:5),'/',tarfile,'.gz   --no-check-certificate']);
     %  webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v4.1/2m/',name(1:5),'/',tarfile,'.gz   --no-check-certificate']);
        if pubres==2
        webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/latest/2m/',name(1:5),'/',tarfile,'.gz   --no-check-certificate']);
	elseif pubres==32
        webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/latest/32m/',name(1:5),'/',tarfile,'.gz   --no-check-certificate']);
	end
        system(['wget  ',webtile,' -o downloadlog'])
	system(['gunzip -f ',tarfile,'.gz']);
        system(['tar -xvf ',tarfile]);
        system(['rm  ',tarfile]);
        %collect all downloaded mosaic dems to tiledirnew/ % to do
        %system(['mv *dem_meta.txt  *reg.txt *.tar *_reg_dem.tif *_reg_matchtag.tif ',tiledirnew]);
	system(['mv ',tarfile(1:end-4),'* ', tiledirnew]);
        [status , cmdout ]=system(['find ',tiledirnew,' -name ',tilefile]);
        tilefile=deblank(cmdout);
        if ~exist(tilefile,'file')
          continue
        end
    end

    if flago==1
    try
    data=readGeotiff(tilefile,'map_subset',rang0);
    data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
    
    tz =interp2(data.x,data.y,double(data.z),xout,yout','*linear',nan);%toc;%4s
    M=~isnan(tz);
    data0.z(M)=tz(M);
    tz =interp2(data.x,data.y,double(data.z),xoutr,youtr','*linear',nan);%toc;%4s
    M=~isnan(tz);
    data0r.z(M)=tz(M);
%     dataref(icount)=data;
    catch e
        fprintf(['\n getmosaic.m There was an error! The message was:',e.message,'; rang0= ',num2str(rang0),'tilefile:',tilefile]);
    end
    end %if
    end 

end

if flago==1
data0.z(isnan(data0.z))=-9999; %return to -9999 % for 7 parameter coregistration.
data0r.z(isnan(data0r.z))=-9999;

if 0
ofile='demmosaic.tif';
projstr='polar stereo north';
writeGeotiff(ofile,data0r.x,data0r.y,data0r.z,4,0,projstr)

%convert cooridnates to longitude latitude. 
[status, cmdout]=system('gdalwarp demmosaic.tif demmosaic_lat.tif -t_srs epsg:4326');
end %if 0

end

if 0
ofile='demmosaic.dat';
data=data0;
[X,Y]=meshgrid(data.x,data.y);
[LAT,LON]=polarstereo_inv(X,Y,[],[],70,-45);
output=[LAT(:),LON(:),double(data.z(:))];
%save -ascii demTolbachik2m.dat output
save(ofile,'output','-ascii')
end

if flagplot==1
nsr=round(100/resr);
hills=hillshade(double(data0r.z(1:nsr:end,1:nsr:end)),data0r.x(1:nsr:end),data0r.y(1:nsr:end),'plotit');
set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
colorbar;%caxis([0 250])
hold on
title(['Reference DEM tile mosaic']);
axis equal
saveas(gcf,'RefDEM','fig')
save demtile.mat data0r
end
%save demtile.mat data0r -v7.3

return
end
