% Get the shapefiles for all tiles in Greenland
addpath(genpath(['/home/dai.56/arcticdemapp/landslide/code1/']));
flag=2; %1 read box from file (e.g., 13_44/input.txt);
	%2 given a big region, get all the standard 50 km tiles.
rang0=[-2031000 -1628000 -637000 -300000]; % fill the minx maxx miny maxy if flag==2.

if flag==1
filename='tilenamelist.txt';
fid = fopen(filename);
n = linecount(fid);

xot={}; yot={}; fileot={};

fid = fopen(filename);
for i=1:n
   count=i;

   ifile=[fgetl(fid)]

   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   [demdir1,dirfilename,ext] =fileparts([demdir]);
   
   % get the boundary from xml file
   c=textread(ifile,'%s','delimiter','\n');
   Xbs=c{2};
   rangei=strread(Xbs,'%d');
   xi=[rangei(1) rangei(2) rangei(2) rangei(1) rangei(1) ];yi=[rangei(4) rangei(4) rangei(3) rangei(3) rangei(4) ];

   xot{count}=xi;yot{count}=yi; fileot{count}=dirfilename;

end

elseif flag==2;

dx=100e3;x0=-4000e3;y0=-4000e3;dxs=dx/2; %50km
xeq=rang0(1);yeq=rang0(3);
xid=floor((xeq-x0)/dx)+1; yid=floor((yeq-y0)/dx)+1;
xids=floor((xeq-(x0+(xid-1)*dx))/(dx/2))+1;
yids=floor((yeq-(y0+(yid-1)*dx))/(dx/2))+1;
xid_start=xid; yid_start=yid;

	xeq=rang0(2);yeq=rang0(4);
	xid=floor((xeq-x0)/dx)+1; yid=floor((yeq-y0)/dx)+1;
	xids=floor((xeq-(x0+(xid-1)*dx))/(dx/2))+1;
	yids=floor((yeq-(y0+(yid-1)*dx))/(dx/2))+1;
	xid_end=xid; yid_end=yid;

xot={}; yot={}; fileot={};
count=0;

	for xid=xid_start:xid_end
        for yid=yid_start:yid_end
         for xids=1:2
             for yids=1:2
count=count+1;
		tilefile=sprintf('%02d_%02d_%01d_%01d.tif',yid,xid,xids,yids);
                tilefile=deblank(tilefile);
	dirfilename=tilefile;

width=100;%2e3;%7e3; %buffer width of the a priori coastline, e.g., 2km.

    x=x0+(xid-1)*dx+(xids-1)*dx/2;y=y0+(yid-1)*dx+(yids-1)*dx/2;
    rang0=[x-width x+dxs+width y-width y+dxs+width]; %tile boundary with buffer width

xi=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];yi=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];


   xot{count}=xi;yot{count}=yi; fileot{count}=dirfilename;

	end
	end
	end
	end
end


shp1 = struct('Geometry', 'PolyGon','X', xot, 'Y', yot,'filename',fileot);
ofile=['tilenames.shp'];
if ~isempty(shp1)
shapewrite(shp1, ofile);
end

