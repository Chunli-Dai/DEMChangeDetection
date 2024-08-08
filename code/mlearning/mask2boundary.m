function [shp1]=mask2boundary(xout,yout,Modfil,Medgs1,ofile1);
%Given water mask Modfil (logical; 1 water, 0 land) get its boundary.
%      or given change mask Modfil (logical; 1 No change, 0 Change) 
%Modified from /home/dai.56/arcticdemapp/coastline/codec2/CoastTileMono.m
%input: Medgs1, a matrix (logical), 1 is the void area to be removed.
% xout, yout: polar stereographic coordinates of the matrix Modfil

constant  %get cloudarea
% cloudarea=smlarea; %m
smlarea=smlarea 
resr=resr

%boundary of box
Medgstb=false(size(Modfil));
Medgstb(:,1)=1;Medgstb(:,end)=1;Medgstb(1,:)=1;Medgstb(end,:)=1;%boundary of the box
Medgs1=Medgstb|Medgs1;

Co=[];
% ofile1='boundary.shp';

if isempty(xout)||isempty(yout)
   fprintf('\n mask2boundary.m xout/yout is empty. \n')
   shp1=[]; return;
else
resx=mean(xout(2:end)-xout(1:end-1));resy=mean(yout(2:end)-yout(1:end-1));
resr=mean([abs(resx),abs(resy)]);
end

nsuby=length(yout);nsubx=length(xout);

nsml=round(smlarea/resr/resr)
%Modfil (logical; 1 water, 0 land) or (logical; 1 No change, 0 Change) 
Modfil = ~bwareaopen(~Modfil, round(smlarea/resr/resr)); %bwareaopen function will remove small clusters (small land areas).

% bwboundaries gets the boundary of 1s (which means land for coastline application).
B = bwboundaries(~Modfil,'noholes'); %^_^ This line will keep the unwanted lake islands (and this can be avoided by remove small land clusters (as shown in the previous line)).
% For change detection, we do not want the boundary of holes.
% B = bwboundaries(~Modfil); % Delta area separated from ocean by road is treated as a lake. We want to keep all large lakes. This line will have the boundary of holes (lakes). 
n=length(B);xo=cell(n,1);yo=cell(n,1);zo=cell(n,1);
BoundingBox=cell(n,1);zone=cell(n,1);SHAPE_Area=cell(n,1);Feature_Ty=cell(n,1);
%clear prob Medgsib Medgstb
tic
peri=sqrt(4*pi*smlarea/resr/resr); % shortest perimeter of a given area.
for k=1:n %can be slow if poor data quanlity, lots of scatterred points;workstereoprob2wocoregcrop/55_06_2_1_coast_v1.0.shp includes island of lakes!
    xid=B{k}(:,2); yid=B{k}(:,1);zid=ones(size(xid));
    if (length(xid)<peri);continue;end %to be immune to bug 11.
    %get the boundaries in terms of mask
    BWb=zeros(size(Modfil));BWb((xid-1)*nsuby+yid)=1;
    Mt=BWb&Medgs1;%filter out edges
    ne=sum(Mt(:));
    if ne~=0 %find the id and replace them with NaN;
       [idy,idx]=find(Mt==1);
       
       for j=1:ne
          Ml=xid==idx(j)&yid==idy(j);
          zid(Ml)=0;%fall into edges
       end
       
    end
    x=xout(xid);y=yout(yid);    
    [LAT,LON]=polarstereo_inv(x,y,[], [],70,-45);
    LAT(zid==0)=nan;LON(zid==0)=nan;%polygon but with only valid polylines displayed.%fix Bug 16:
    zidi=find(zid==0);zidid=zidi(2:end)-zidi(1:end-1);
    %id=find(zidid==1);%delete sequential nans
%    idx=zidi(id+1);
    idb=find(zidid~=1);M=zidi;M(idb+1)=nan;
    M(idb)=nan; %fix Bug 16: connecting lines undesirably
%   idxb=zidi(idb+1); LAT(idxb)=nan;LON(idxb)=nan;
    [idx,idn]=separatelines(M,10); %fix Bug 15
    %figure;plot(x,y,'go');hold on;plot(x(idx),y(idx),'ro-')
    LAT(zidi(idn))=nan;LON(zidi(idn))=nan;LON(idx)=[];LAT(idx)=[];
    
    %get area
    Mb=poly2mask(xid,yid,nsuby,nsubx); 
    area=sum(Mb(:))*resx*abs(resy); %poly2mask
    
    datamask.x=[]; datamask.y=[];datamask.z=[];
    xo{k}=LON;yo{k}=LAT; zo{k}=datamask;
    BoundingBox{k}=[min(LON) min(LAT);max(LON) max(LAT);];
    zone{k}='scar';
    OBJECTID{k}=[];%place holder
    Feature_Ty{k}='None';%place holder
    SHAPE_Leng{k}=[];%place holder
    SHAPE_Area{k}=area;%place holder
    Feature_ID{k}=[];%place holder
    Connectivi{k}=[];%place holder
    region{k}=[];%place holder
end
fprintf('Retrieve boundary')
toc
clear Mt BWb Medgs1
idd=find(cellfun(@isempty,xo)); %fix bug 12
xo(idd)=[];yo(idd)=[]; zo(idd)=[];
BoundingBox(idd)=[];zone(idd)=[];Feature_Ty(idd)=[];SHAPE_Area(idd)=[];

if 0
    out=[];
    for i=1:length(xo)
    out=[out;yo{i}(:), xo{i}(:)]; %lat lon
    end
    save -ascii outline.dat out

    %save xoyo.mat xo yo -v7.3
end

if 0
    figure;
    hold all
    for k=1:length(xo)
    plot(xo{k},yo{k},'.-')
    end
end

% shp1 = struct('Geometry', 'PolyGon', 'X', xo, 'Y', yo);
% shp1 = struct('Geometry', 'PolyLine', 'X', xo, 'Y', yo);
%shp1 = struct('Geometry', 'PolyGon','BoundingBox',BoundingBox, 'X', xo, 'Y', yo,'Z',zo ...
shp1 = struct('Geometry', 'PolyGon','BoundingBox',BoundingBox, 'X', xo, 'Y', yo ...
    , 'Feature_Ty',Feature_Ty,'SHAPE_Area',SHAPE_Area,'zone',zone);

if ~isempty(shp1)
 shapewrite(shp1, ofile1); % May 2023
end

return
end
