%addpath(genpath(['/home/dai.56/arcticdemapp/landslide/code1/']));

%get mask from mask mask/GimpIceMask_90m.tif;
%rang0=[-268840  -228840 -2772800 -2732800];

%[tag]=getrgi(rang0);
function tag = getrockicegreenland(rang0);
infile='/blue/chunlidai/apps/landslide/code1/mask/GimpIceMask_90m.tif';
tag=readGeotiff(infile,'map_subset',rang0);
infile='/blue/chunlidai/apps/landslide/code1/mask/GimpOceanMask_90m.tif';
ocean=readGeotiff(infile,'map_subset',rang0);
tag_org=tag;

tag.z=~(tag.z)&~(ocean.z); %1 rock, non-ice and non-water; 0 non-rock, either ice or water.

%to be compatible with other codes.
%non-ice 1 ; ice 0;
tag.ice=~(tag_org.z);

%1 water; 0 land
tag.water=(ocean.z);

%save('BarnesrockRGI_sw.mat','tag','-v7.3')
%save('rockmask.mat','tag','-v7.3')

%% write output
if 1
try
projstr='polar stereo north';
OutName=['rockmask.tif']; % %1 is good (landslide); 0 is bad
writeGeotiff(OutName,tag.x,tag.y,int32(tag.z),3,0,projstr)
catch e
            fprintf('getrockicegreenland.m There was an error! The message was:\n%s',e.message);
	    save tag_wrong.mat -v7.3
end

end

return
end
