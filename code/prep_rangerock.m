function [Co]=prep_rangerock(rang0)
constant
flagicewater=1; %1 output would be both ice and water mask; 2, output would be ice map only.
Co=[];
%addpath(genpath(['/home/dai.56/arcticdemapp/landslide/code1/']));

%get rang0
%Only work on subset of all images around the scarp;
%rang0=[431 445 6772 6784]*1e3; %UTM '6N' coordinates;
rang0=rang0+[-4e3 +4e3 -4e3 +4e3]; 
%z1 ='6N';

%utm to polar stereographic coordinates;
		if 0
               %[ellipsoid,estr] = utmgeoid(z1);
                utmstruct = defaultm('utm'); 
                utmstruct.zone = z1; 
               %utmstruct.geoid = ellipsoid; 
                utmstruct.geoid = wgs84Ellipsoid;  %consistent with landsat image geoid
                utmstruct = defaultm(utmstruct);
		end

if strcmp(projgdal,'epsg:3413') %polar coordinates for input range
   rang0_polar=rang0;
else % UTM coordinates for input range
x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
[lat,lon]=minvtran(utmstruct,x0,y0);
[Xb,Yb]=polarstereo_fwd(lat,lon,[],[],70,-45);
rang0_polar=round([min(Xb) max(Xb) min(Yb) max(Yb)]/40)*40;
end

%get rock mask 
%refer to /fs/project/howat.4/dai.56/chunliScripts/scripts/plotsrtmvs.m

if flaggreenland==1 %1 use GIMP masks for greenland; 0 otherwise
  [tag]=getrockicegreenland(rang0_polar); % z %1 rock, non-ice and non-water; 0 non-rock, either ice or water.
  tag.z=tag.ice;
  water.x=tag.x;water.y=tag.y;water.z=tag.water; %1 ,water; 0 non-water
else

[tag]=getrgi(rang0_polar); %rock/non-ice 1 ; ice 0;

[water]=getocean(rang0_polar); %1 water; 0 land
end

%in case masks have only one dimension (failed), output all rock mask (no filter will be applied).

[nt,mt]=size(tag.z);
[nw,mw]=size(water.z);
if min([nt,mt,nw,mw])<=2
    data.x=rang0(1):10:rang0(2);data.y=rang0(4):-10:rang0(3);
    icewater.x=data.x;icewater.y=data.y;
    % 1 ice or water; 0 rock
    icewater.z=false(length(data.y),length(data.x));
    Co=icewater;
    fprintf('\n prep_rangerock.m ice/water mask retrieval failed.\n');
    save icewater_wrong.mat -v7.3
    return
end

if strcmp(projgdal,'epsg:3413') %polar coordinates for input range
	data=tag; % 0 ice; 1 non-ice;
	water2.z=interp2(water.x,water.y,double(water.z), data.x,data.y','*linear',1); %1 water;
%	tag=data; %1 rock, 0 ice;
else % UTM coordinates for input range
%Polar stereographic coordinates to UTM
%data.x=rang0(1):40:rang0(2);data.y=rang0(4):-40:rang0(3);
data.x=rang0(1):10:rang0(2);data.y=rang0(4):-10:rang0(3);
[X,Y]=meshgrid(data.x,data.y);
[lat,lon]=minvtran(utmstruct,X,Y);
[X,Y]=polarstereo_fwd(lat,lon,[],[],70,-45);
data.z=interp2(tag.x,tag.y,double(tag.z), X,Y,'*linear',0); %0 ice
water2.z=interp2(water.x,water.y,double(water.z), X,Y,'*linear',1); %1 water;
tag=data; %1 rock, 0 ice;
end

%get water mask;

icewater=data; %x y vectors
%flagicewater=2; %1 output would be both ice and water mask; 2, output would be ice map only.
if flagicewater==1
icewater.z=~tag.z|water2.z; % 1 ice or water; 0 rock
elseif flagicewater==2
	icewater.z=~tag.z; % 1 ice; 0 non-ice.
end
%save  icewater.mat   icewater tag water2   -v7.3  

Co=icewater;

%% write output
if 0
%projstr='polar stereo north';
projstr=projstrin; 
%projstr='UTM';%'Transverse_Mercator'; %'UTM_Zone_6N';
%zone='6N';
OutName=['icewater.tif']; % %1 is good (landslide); 0 is bad
writeGeotiff(OutName,icewater.x,icewater.y,int32(icewater.z),3,0,projstr)
end

return
end
