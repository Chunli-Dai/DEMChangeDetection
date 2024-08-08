function [profc,planc,Kt,totalc] = curvature2(DEM)

% 8-connected neighborhood curvature of a digital elevation model 
%
% Syntax
%
%     [profc,planc,Kt, totalc] = curvature2(DEM)
% Using the equations in Wilson & Gallant, 2000, pp. 57.
%  Reference: Wilson, J.P. and Gallant, J.C. eds., 2000. Terrain analysis: principles and applications. John Wiley & Sons.
% Note: results seem to be similar to curvature.m

if 0 % test example in DEM Surface Tools for ArcGIS_A4.pdf.
    x=linspace(0,2*pi,100); y=x;
    clear DEMt; DEMt.x=x;DEMt.y=y;
    [X,Y]=meshgrid(DEMt.x,DEMt.y);
    a=1;b=1;c=1;d=1;
    z=a.*sin(X+ b.*sin(Y))+ c.*sin(X)+ d;
    DEMt.z=z;
    figure;surf(DEMt.x,DEMt.y,DEMt.z); colorbar;title('DEM test')
    
    [profc,planc,Kt,totalc] = curvature2(DEMt);
    
    %compare with curvature.m
    resx=mean(DEMt.x(2:end)-DEMt.x(1:end-1));resy=mean(DEMt.y(2:end)-DEMt.y(1:end-1));
    resr=mean([abs(resx),abs(resy)]);
    [profc2,planc2] = curvature(DEMt.z,resr,'both');
    
%     [K,H,Pmax,Pmin] = surfature(X,Y,DEMt.z);%NOT compatible.
    
    figure;imagesc(DEMt.x,DEMt.y,profc);colorbar;title('profile curvature test'); view(0,-90)
    figure;imagesc(DEMt.x,DEMt.y,planc);colorbar;title('plan curvature test'); view(0,-90)
    figure;imagesc(DEMt.x,DEMt.y,(Kt));colorbar;title('tangential curvature test'); view(0,-90)
    figure;imagesc(DEMt.x,DEMt.y,totalc);colorbar;title('total curvature test'); view(0,-90)
    
end

% [X,Y]=meshgrid(DEM.x,DEM.y);
Z=DEM.z;
resx=mean(DEM.x(2:end)-DEM.x(1:end-1));resy=mean(DEM.y(2:end)-DEM.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);

% First Derivatives
% [Xu,Xv] = gradient(X);
% [Yu,Yv] = gradient(Y);
[Zu,Zv] = gradient(Z,resr);

% Second Derivatives
% [Xuu,Xuv] = gradient(Xu);
% [Yuu,Yuv] = gradient(Yu);
[Zuu,Zuv] = gradient(Zu,resr);

% [Xuv,Xvv] = gradient(Xv);
% [Yuv,Yvv] = gradient(Yv);
[Zuv,Zvv] = gradient(Zv,resr);

p=Zu.^2+Zv.^2;
q=p+1;

%Profile curvature 
Kp=(Zuu.*Zu.^2+2*Zuv.*Zu.*Zv+Zvv.*Zv.^2)./(p.*q.^1.5);
profc=Kp; % Profile curvature is positive for concave topography, negative for convex topography.

%Plan (contour) curvature
% Plan curvature is negative for diverging flow (on ridges) and positive for converging flow (valleys).
Kc=(Zuu.*Zv.^2-2*Zuv.*Zu.*Zv+Zvv.*Zu.^2)./p.^1.5;
planc=Kc;

%Tangential curvature
% Mitasova and Hofierka (1993) suggested Kt is better than plan curvature. Because it does not take on extremely large values when slope is small.
% 
Kt=(Zuu.*Zv.^2-2*Zuv.*Zu.*Zv+Zvv.*Zu.^2)./(p.*q.^0.5);

%Total curvature
totalc=Zuu.^2+2*Zuv.^2+Zvv.^2;

return
end