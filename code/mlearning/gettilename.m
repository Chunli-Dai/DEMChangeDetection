function [tilename]=gettilename(x,y,flag)
%Given coordinates, output tilename
%flag: 1, polar stereographic coordinates; 2; lon lat.
addpath(genpath('/Users/chunlidai/ArchivedBoxSync/ESI2019/landslide/code1'));

if ~exist('flag','var')
    flag=1; %default
else
    % flag =1 or 2;
    if length(flag)==1
        if ~(flag==1 ||flag==2)
           fprintf(['\n Input flag must be 1 (polar stereographic coordinates) or 2 (longitude/latitude coordinates).\n'])
           return
        end
    else
        fprintf(['\n Input flag must be a scalar, 1 (polar stereographic coordinates) or 2 (longitude/latitude coordinates).\n'])
        return
    end

end
if flag==1
    fprintf(['\n Input stereographic coordinates (x, y)=',num2str([x, y]),'\n'])
    xeq=x;yeq=y;
elseif flag==2
    lon=x;lat=y;
    [xeq,yeq]=polarstereo_fwd(lat, lon,[],[],70,-45);
    fprintf(['\n Input geographic coordinates (lon, lat)=',num2str([lon, lat]),'\n'])
end

dx=100e3;x0=-4000e3;y0=-4000e3;
dxs=dx/50; %2km

xid=floor((xeq-x0)/dx)+1; yid=floor((yeq-y0)/dx)+1;
xids=floor((xeq-(x0+(xid-1)*dx))/(dx/2))+1;
yids=floor((yeq-(y0+(yid-1)*dx))/(dx/2))+1;
xidss=floor((xeq-(x0+(xid-1)*dx+(xids-1)*dx/2))/dxs)+1;
yidss=floor((yeq-(y0+(yid-1)*dx+(yids-1)*dx/2))/dxs)+1;
tilefile=sprintf('%02d_%02d_%01d_%01d_%02d_%02d.tif',yid,xid,xids,yids,xidss,yidss) ;

[~,tilename,~]=fileparts(tilefile);

fprintf(['\n Tile name=',tilename,'\n'])

return
end
