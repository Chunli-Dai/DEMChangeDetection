function datao=getslumptif(Sti)
% find the corresponding slumptif file (e.g., saved in slumptif_region09_vvl), and read the data.
% Sti, the ith polygon in the input training shapefile.
% datao, output, same as output from box2mosaic.m, read from e.g., slumptif_region09_vl/slump447.tif

constant

[demdir,name,ext] =fileparts([strtrim(regiondir)]);
regid=name(11:12);

regionstr=['region',regid]; %'region09';
indir=['../',regionstr,'/'];

% Find the class type based on area. 

areaj=Sti.SHAPE_Area;

if areaj>=2000&areaj<5000;
   classtype='_vs';
elseif areaj >=5000 & areaj<10e3;
   classtype='_s';
elseif areaj >=10000 & areaj<50e3;
   classtype='_m';
elseif areaj >=50000 & areaj<100e3;
   classtype='_l';
elseif areaj >=100000 & areaj<500e3;
   classtype='_vl';
elseif areaj >=500e3 
   classtype='_vvl';
end

if strcmp(classtype,'_vvl')
   classtypein='_landslide';
else
   classtypein=classtype;
end

%load original shapefile refers to Dividemap_candpar_vvl.m

%shapefile=[indir,'change2polyg_',regid,classtypein,'.shp';]
shapefile=[indir,'change2polyg_',regid,classtypein,'.shp'];
S1=shaperead(shapefile);
shpc=struct2cell(S1)';
areaj=[shpc{:,6}]';
M1=areaj>=2000;
S1s=S1(M1);S1=S1s;

nx=length(Sti.X);
%find the matched slumptif file
id=find(areaj==Sti.SHAPE_Area);
matchid=[];
for j=1:length(id);
	k=id(j);
	if nx==length(S1(k).X)	
		%M2=Sti.X==S1(k).X; %NaN ~=NaN
		M2=abs(Sti.X-S1(k).X);
		if nansum(M2)< 0.1/3600 % all elements are the same. 3 m ~ 0.1 arcsecond
		   matchid=k;
   	fprintf(['\n getslumptif.m, find the matched id, [j matchid]: ',num2str([j matchid])]);
	        end
	end
end

if isempty(matchid)
   fprintf(['\n Error: getslumptif.m ',regionstr, ': cannot find the matched slumptif file. ']);
   datao=[];
   save test1.mat -v7.3
else
	% ../region09/slumptif_region09_vl/slump447.tif
   file=[indir,'slumptif_',regionstr,classtype,'/slump',num2str(matchid),'.tif'];
   datao=readGeotiff(file);
end

return
end
