function [XYbi,rangei]=dembd(infile)
%input dem file name
% get XYbi, boundary of valid data;
%see also imagebd.m, mask2boundary.m

data=readGeotiff(infile);
resr=mean((diff(data.x)));

%reduce resolution of data
M=~(data.z==-9999|data.z==0|data.z==32767); %data valid
% M=~(data.z==-9999|data.z==32767); %data valid
datar.z=imresize(M, resr/200);%200m resolution
datar.x=imresize(data.x,resr/200);
datar.y=imresize(data.y,resr/200);
xout=datar.x;yout=datar.y;

Modfil=datar.z;
B = bwboundaries(Modfil,'noholes'); 

n=length(B);
if n>1; fprintf('dembd.m has more than one polygon (maybe multiple islands)! \n');end

XYbi=[];
for k=1:n 
    xid=B{k}(:,2); yid=B{k}(:,1);zid=ones(size(xid));
    x=xout(xid);y=yout(yid);  
    Xb=x;Yb=y;
%     XYbi=[XYbi;[x(:),y(:);[NaN NaN]]];
    XYbi=[XYbi;x(:),y(:);];
%     hold on;plot(XYbi(:,1)*1e-3,XYbi(:,2)*1e-3,'ro-')
end

    if 0 %get the convex hull of all points -> may have too many ocean for multiple islands;
        P = XYbi;
        [k,av] = convhull(P);
    %     figure;plot(P(:,1),P(:,2),'s');hold on;plot(P(k,1),P(k,2))
        XYbi=[P(k,1),P(k,2)];

    else %find the polygon with longest perimeter
        for k=1:n ; msize(k)=length(B{k}(:,1));end
        [~,k]=max(msize); 
        XYbi=[];
        xid=B{k}(:,2); yid=B{k}(:,1);zid=ones(size(xid));
        x=xout(xid);y=yout(yid);  
        Xb=x;Yb=y;
        XYbi=[XYbi;x(:),y(:);];
    end

        rangei=[min(Xb) max(Xb) min(Yb) max(Yb)];
return
end