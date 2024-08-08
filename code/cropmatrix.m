function [wm]=cropmatrix(data,M)
%crop matrix using the M matrix (1 valid data, 0 non valid data); 
%data.z and M should have the same size;
% jump=M;
xout=data.x;yout=data.y;

[n1,m1]=size(data.z);

if max([n1,m1]) < 1000 %small matrix
    [X,Y]=meshgrid(xout,yout);
    xmin=min(X(M==1));xmax=max(X(M==1));
    ymin=min(Y(M==1));ymax=max(Y(M==1));
else %big matrix
    [idy1,idx1]=find(M); %data0.z(idy1(1),idx1(1))
    xmin=xout(min(idx1));xmax=xout(max(idx1));
    ymax=yout(min(idy1));ymin=yout(max(idy1)); %assume y descending
end

if isempty(xmin)||isempty(xmax)||isempty(ymin)||isempty(ymax)
    wm.x=[];wm.y=[];wm.z=[];
    return;
end

idx=xout>=xmin&xout<=xmax;idy=yout>=ymin&yout<=ymax;
wm.x=xout(idx);wm.y=yout(idy);wm.z=data.z(idy,idx);

return
end