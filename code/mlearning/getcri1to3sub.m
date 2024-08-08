function [masko]=getcri1to3sub(ifile)

% refers to modifyshapefile.m

smlarea=1e3*2; %1000 m^2

data0=readGeotiff(ifile); %band 3: change (not scaled), curvature (scaled), time (scaled)
jump=data0;jump.z=[];jump.z=double(data0.z(:,:,1));

rang0=[min(data0.x) max(data0.x) min(data0.y) max(data0.y)]; 

try
[dem]=getmosaic(rang0);
catch e
   fprintf('\n getcri1to3sub.m There was an error! The message was:\n%s',e.message);
   masko=jump;masko.z=false(jump.z);
   return
end

demr.z=interp2(dem.x,dem.y,double(dem.z),jump.x,jump.y','nearest',nan);
demr.x=jump.x;demr.y=jump.y;

resx=mean(data0.x(2:end)-data0.x(1:end-1));resy=mean(data0.y(2:end)-data0.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);
% Get the examed slump area (as ground truth).
Mf4=jump.z<=-1; 
Mf4clean = bwareaopen(Mf4, smlarea/resr/resr); %remove small clusters
Mn=jump;Mn.z=Mf4clean;

Mp1=jump.z>0.5; %>1m; Positive clusters area >=1000 m^2
Mpclean = bwareaopen(Mp1, smlarea/resr/resr/2); %remove small clusters
clear Mp
Mp=jump;Mp.z=Mpclean;

[maskcri1to3]=dhelev(jump,demr,[],[],Mn,Mp);

masko=jump;masko.z=maskcri1to3;

return
end

