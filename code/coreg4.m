function [dX4Sg,idregion,data0r]=coreg4(rang0,idregion,XYbg,f,fdir,tag)
%Modified from /home/dai.56/arcticdemapp/landslide/code1/coreg3.m
%refers to /home/dai.56/arcticdemapp/volcano/code/volcano.m
%Input:
%       rang0 : the boundary of target zone;
%       idregion: index of strip DEMs within the target zone;
%       XYbg:     the boundary of the actual data of each strip DEM;
%       f,fdir:   the filename and directory of all available strip DEMs;
%       tag: a priori rock surfaces.
%
%Modification (January 2020): 
%          Put the coregistration part of code volcano.m to this subroutine;  
%	   This function 1) finds the best reference DEM; 2) coregister all rest DEMs to the reference DEM.

close all
constant

odircoreg='./outcoregdem/';
if ~exist(odircoreg,'dir')
    mkdir(odircoreg)
end

params.I = [1 0 0 0 0 0 0];
params.C = [50 0.001 0.001 0.05 0.0001];
%         params.G = [3000 20];
params.M = '3D';
params.V = [10 20 10];
params.G = [9000 20]; %Adjust the max height parameter for the 2012 Kamchatka Volcano
coregflag=8;%1 parameter (vertical), 3 parameter (Ian's); 7 parameter (MJ's); 8, MJ's setsm coregistration with 3 parameters.
flagplot=1;

% find the DEM mosaics within rang0.
%xr0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];yr0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];

resrc=200;
%resr=2;%10;%40; %2; % 2 m cause out of memory for 1.7GB DEM strip for 3 parameters (Nuth and Kabb) coregistration.
resr=40;

save test4.mat -v7.3
%exit

if 0
ranget=round(rang0/resrc)*resrc;rang0=ranget;
tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);
xoutr=tx;youtr=ty;
data0r.x=xoutr;data0r.y=youtr;data0r.z=nan(length(youtr),length(xoutr));
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
xout=tx;yout=ty;
data0.x=xout;data0.y=yout;data0.z=nan(length(yout),length(xout));
end

dx4=zeros(3,1); % do not use the tile reg.txt data.

if length(idregion)<2;fprintf(['/n There are less than 2 DEMs; hence there is no need for coregistration. /n']);return;end

fprintf ('\n Step 1.1b: Coregister all pairs of strip DEMs.')
tic
dzxyd=zeros(length(idregion),3);dzxy=cell(size(idregion));
% dxov=resr*10; %grid size for approximating overlapping area of polygons.
idreg=find(dzxyd(:,1)~=0);  %reg file
idregn=find(dzxyd(:,1)==0); %no reg file
pg=zeros(length(idregion),3);rmsreg2=zeros(length(idregion),1);
dX4Sg=pg;
if coregflag==7
%   clear dX4Sg
end
idd=[];idd2=[];

ndem=length(idregion);
count=0; %record the count of coregistered pairs;

id=idregion;
t=zeros(length(id),1);    
for j=1:length(id)%0%length(id)
        ymd=f{id(j)}(6:13);i=id(j); 
        t(j)=datenum(ymd,'yyyymmdd');
%display(['Overlapping Files : ',f{id(j)}])
end
[~,idsort]=sort(t);id=id(idsort); %sort the id based on time
idregion=idregion(idsort);

for j=1:length(id)%0%length(id)
display(['Overlapping Files ',num2str(j),': ',f{id(j)}])
end

for j=1:length(id)
XYbj=XYbg{id(j)};
Xb=XYbj(:,1);Yb=XYbj(:,2);
rangeid(j,1:4)=[min(Xb) max(Xb) min(Yb) max(Yb)];
end

%rang0sov=[max(range(id,1)) min(range(id,2)) max(range(id,3)) min(range(id,4))];
rang0sov=[max(rangeid(:,1)) min(rangeid(:,2)) max(rangeid(:,3)) min(rangeid(:,4))];
ranget=[rang0sov/resr];
rang0sov=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
rang0st=[max(rang0sov(1),rang0(1)) min(rang0sov(2),rang0(2)) max(rang0sov(3),rang0(3)) min(rang0sov(4),rang0(4))]; %overlapping grid
x0st=[rang0st(1) rang0st(2) rang0st(2) rang0st(1) rang0st(1) ];y0st=[rang0st(4) rang0st(4) rang0st(3) rang0st(3) rang0st(4) ];

fprintf (['\n Overlapping zone (rang0sov):',num2str(rang0sov),'.\n'])
mov=floor((rang0sov(2)-rang0sov(1))/resr);
nov=floor((rang0sov(4)-rang0sov(3))/resr);
if mov<=1||nov<=1
   fprintf (['\n One dimension of the overlapping zone is <=1; [mov nov]=:',num2str([mov, nov]),'.\n'])
   fprintf('\n Use rang0 as rang0sov.');
   rang0sov=rang0;
   rang0st=[max(rang0sov(1),rang0(1)) min(rang0sov(2),rang0(2)) max(rang0sov(3),rang0(3)) min(rang0sov(4),rang0(4))]; %overlapping grid
   x0st=[rang0st(1) rang0st(2) rang0st(2) rang0st(1) rang0st(1) ];y0st=[rang0st(4) rang0st(4) rang0st(3) rang0st(3) rang0st(4) ];
end

%display('Step 2.1: searching for the best reference DEM.')
%find the reference DEM: 1, in the overlapping zone maximum good points
% 2, Distance between center of overlapping and center of reference DEM
% shortest
ymd=[]; jref=[];
idd=[];
% nptsub=zeros(length(id),1);nptall=zeros(length(id),1);
nptsubraw=zeros(length(id),1);nptallraw=zeros(length(id),1);dist=nptsubraw;
score=zeros(length(id),3);rockn=zeros(length(id),1);
enlref=0.5;
tic 
p1=[(rang0sov(1)+rang0sov(2))/2 (rang0sov(3)+rang0sov(4))/2 ];
for j=1:length(id)
ymd=f{id(j)}(6:13);i=id(j); 
%t(j)=datenum(ymd,'yyyymmdd');
%data=datarsv((isv==i));

%read DEM 1
i=id(j);
XYbi=XYbg{i};
metafile=[fdir{i},'/',f{i}]; metafileref=metafile;
clear datar
[datar,nptsubrt]=readdem(XYbi,metafile,rang0,resr);
data=datar;

%checking rock coverage
if 0
rangeov=[min(data.x) max(data.x) min(data.y) max(data.y)];
idx=find(tag.x>=rangeov(1) & tag.x<=rangeov(2)); %assume rock data in 40m resolution
idy=find(tag.y>=rangeov(3) & tag.y<=rangeov(4));
rocktag=tag.z(idy,idx);
rockfilter=rocktag&data.z~=-9999;%&abs(dz)>1e-4;
rockn(j)=sum(rockfilter(:));
end
% Checking the overlapping coverage
idx=find(data.x>=rang0sov(1) & data.x<=rang0sov(2));
idy=find(data.y>=rang0sov(3) & data.y<=rang0sov(4));
demo=data.z(idy,idx);
%nptsubrt=sum(sum(demo~=-9999))*(str2double(res))^2/(dx*dx*4);
nptsubrt=sum(sum(demo~=-9999))/(length(demo(:)));
nptsubraw(j)=nptsubrt; %in cases, all got skipped
% if nptsubrt>=0.18;idd=[idd;j];end %2nd, the bad points ratio has to be less than 2%.
nptallraw(j)=sum(sum(data.z~=-9999));
p2=[(data.x(1)+data.x(end))/2 (data.y(1)+data.y(end))/2 ];
dist(j)=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);
score(j,1)=nptsubrt*10; %scale 0 to 10
%data=readGeotiff(infile);
if 0
hills=hillshade(double(data.z),data.x,data.y,'plotit');
set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
colorbar;%caxis([0 250])
[X,Y]=meshgrid(data.x,data.y);
hold on
plot(x0st, y0st,'r-','linewidth',6)
plot(xeq, yeq ,'r*','Markersize',12)
axis equal
% plot(X(M),Y(M),'r.')
% %saveas(gcf,['refDEMj',num2str(j)],'fig')
% saveas(gcf,'tt','fig')
end
%No edge within the subzone
% 
%if nptsubrt > 0.99 && ratio<1; jref=j;break;end  %hi, change 8% to 1%
end
% in cases, all got skipped
% nptsubraw(rockn<20)=0;%if no rock, not consider this dem
% do not use SRTM as reference
idd=[];
for j=1:length(id);
ymd=f{id(j)}(6:13);
if ymd=='20000211';idd=[idd;j]; end
if ymd=='20180611';nptsubraw(j)=2; end
% if ymd=='20120920'; idd=j;end
end
nptsubraw(idd)=0 

idmaxsub=find(abs(nptsubraw-max(nptsubraw)) <0.1); %0.01
% [~,jref]=min(dist(idmaxsub));
% [~,jref]=max(nptallraw(idmaxsub));
score(:,2)=(1-dist/max(dist))*10; %scale 0 to 10;
score(:,3)=(nptallraw/max(nptallraw))*10; %scale 0 to 10;
[~,idj]=max(score(idmaxsub,2)+score(idmaxsub,3));
jref=idmaxsub(idj)
jref=1; %hi

if isempty(jref);warning(['No reference DEM found; isel=',num2str(isel)]);jref=1;end
j=jref;  iref=id(j);  %store i before idd
id(idd)=[];%nptsub(idd)=[];nptall(idd)=[];
display('Searching reference DEM ...');
toc % 18 minutes

% 3rd considering the overlapping area with respect to others.

i=iref;
dzxyt=dzxy{idregion==i};
if isempty(dzxyt)
regflag=0; %0 bad 1 good
dxyzref=zeros(1,3);
else
dxyzref=[dzxyt(2:3) dzxyt(1)];
regflag=1;
end

    XYbi=XYbg{i};
    metafile=[fdir{i},'/',f{i}]; metafileref=metafile;
fprintf (['\n Selected reference DEM:',metafileref,'.\n'])
    clear datar
    [datar,nptsubrt]=readdem(XYbi,metafile,rang0,resr);
    data0r=datar;
    sat=filename2sat(metafile); ymd=filename2ymd(metafile);
    textref=[sat,'_',ymd];

    if nptsubrt<0.5 %data quality poor
fprintf (['\n Selected reference DEM has poor data quality:',metafileref,'.\n'])
    end

if 0 %plot the reference DEM
hills=hillshade(double(data0r.z),data0r.x,data0r.y,'plotit');
set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
colorbar;%caxis([0 250])
hold on
plot(x0st, y0st,'r-','linewidth',6)
plot(xeq, yeq ,'r*','Markersize',12)
axis equal
title('RefDEM')
saveas(gcf,'refDEM','fig')
end
%find the overlapping zone for this piece
%mx0=find(xout>=rang0st(1) & xout<=rang0st(2) );
%my0=find(yout>=rang0st(3) & yout<=rang0st(4) );
%reference dem at this piece
datadsx=data0r.x; datadsy=data0r.y ;
idx=find(datadsx>=rang0st(1) & datadsx<=rang0st(2));
idy=find(datadsy>=rang0st(3) & datadsy<=rang0st(4));
nyj=length(idy);nxj=length(idx);
%if (length(my0)~=nyj || length(mx0)~=nxj);warning(['Size mismatch, isel:',num2str(isel)]);end
%xout(1:nxj,isel)=datadsx(idx)';     yout(1:nyj,isel)=datadsy(idy)'; 
%yout(my0)=data0r.y(idy);xout(mx0)=data0r.x(idx);
%collecting reference dem 
[X,Y]=meshgrid(datadsx(idx),datadsy(idy));
[LAT,LON]=polarstereo_inv(X,Y,[], [],70,-45);
% wmp1 = roipoly;data0r.z(wmp1)=NaN; data0r.z(wmp1)=-9999; 
% temporary output of dem
demo=data0r.z(idy,idx);
% output=[output;LAT(:),LON(:),double(demo(:))];
% save -ascii demTyndalbig.dat output


idem1=jref;i=idregion(jref);
   
idall=1:ndem;
idselected=idall(~ismember(idall, jref));
    for idem2=idselected  %1:ndem;
        j=idregion(idem2);  
        XYbj=XYbg{j};
    
        fprintf(['\n Coregistering pair ',num2str([idem1  idem2]),': ', f{i}, ' vs ',f{j},'.\n'])
    
    %Add constraint: for each pair, the overlapping area has to be > 1% of each scene DEM.
    %all space of two DEMs;
    range1=[min([XYbi(:,1);XYbj(:,1);]) max([XYbi(:,1);XYbj(:,1);]) min([XYbi(:,2);XYbj(:,2);]) max([XYbi(:,2);XYbj(:,2);])];
    ranget=round(range1/resrc)*resrc;
    tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);

    idx=round((XYbi(:,1)-tx(1))/resrc)+1;
    idy=round(-(XYbi(:,2)-ty(1))/resrc)+1;
    inref= poly2mask(idx,idy, length(ty),length(tx)); % faster than inpolygon
    arearef=(sum(inref(:)));%pixels

    idx=round((XYbj(:,1)-tx(1))/resrc)+1;
    idy=round(-(XYbj(:,2)-ty(1))/resrc)+1;
    intar= poly2mask(idx,idy, length(ty),length(tx)); % faster than inpolygon
    areatar=(sum(intar(:)));%pixels
    overlaparea=sum(sum(inref&intar));
    ratioref=overlaparea/arearef*100; 
    ratiotar=overlaparea/areatar*100; 
    %e.g. 2km length overlap for a 120 km long strip and 16 km long scene.
    if min([ratioref ratiotar])> 1 && mean([ratioref ratiotar])> 5 %good case
    else
        fprintf(['\n Overlapping area of this pair is too small:',num2str([ratioref ratiotar]),' %.',]);
        continue;
    end
    
    %overlapping boundary
    rangref=[min(data0r.x) max(data0r.x) min(data0r.y) max(data0r.y)];
    rangtar=[min(XYbj(:,1)) max(XYbj(:,1)) min(XYbj(:,2)) max(XYbj(:,2))];
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
    rangeov=round(rangeov/resr)*resr;

    metafile=[fdir{j},'/',f{j}];
    clear datar
    [datar,nptsubrt]=readdem(XYbj,metafile,rangeov,resr);
    
    %in case datar miss one pixel.
    rangtar=[min(datar.x) max(datar.x) min(datar.y) max(datar.y)];
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
    rangeov=round(rangeov/resr)*resr;

    sat=filename2sat(metafile); ymd=filename2ymd(metafile);
    texttar=[sat,'_',ymd];
    
    %all DEMs are checked in dem1; here we only read the overlapping part of dem2;
    %not representing the data quality of whole strip. -> so don't store it in idd2.
    if nptsubrt<0.5 
%         idd2=[idd2;idem2];
        continue
    end    
    
    %crop the overlapping data
    refdem=[];tardem=[];
    idx=find(data0r.x>=rangeov(1) & data0r.x<=rangeov(2));
    idy=find(data0r.y>=rangeov(3) & data0r.y<=rangeov(4));
    refdem.x=data0r.x(idx);refdem.y=data0r.y(idy);
    refdem.z=data0r.z(idy,idx);
    idx=find(datar.x>=rangeov(1) & datar.x<=rangeov(2));
    idy=find(datar.y>=rangeov(3) & datar.y<=rangeov(4));
    tardem.x=datar.x(idx);tardem.y=datar.y(idy);
    tardem.z=datar.z(idy,idx);
    [n,m]=size(tardem.z) 
    minm=min(n,m);
    size(refdem.z)
    if sum(size(tardem.z)~=size(refdem.z))||minm<3
        warning(['Wrong overlapping crop, check i:',num2str([i j])]);
        tardem=refdem;
        tardem.z = interp2(datar.x,datar.y,double(datar.z),refdem.x,refdem.y','*linear',nan); %logical

%         idd=[idd;j];
%         continue
    end
    if exist('tag','var')
        if ~isempty(tag)
	    rocktag = interp2(tag.x,tag.y,tag.z,refdem.x,refdem.y','*nearest',0); %logical
 	    rocktag = logical(rocktag);
        else
        rocktag=[];
        end
    else
        rocktag=[]; %for mp1 in coregisterdems.m
    end
    
    %use 40 meter resolution DEM for coregistration
    %reduce the resolution to 40 m
    ranget=[[min(refdem.x) max(refdem.x) min(refdem.y) max(refdem.y)]/resrc];
    ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resrc;
    tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);
    data=refdem;
    data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
    tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
    tz(isnan(tz))=-9999; %return to -9999
    refdemr= struct();
    refdemr.x=tx;refdemr.y=ty;  refdemr.z=tz;
    data=tardem;
    data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
    tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
    tz(isnan(tz))=-9999; %return to -9999
    tardemr= struct();
    tardemr.x=tx; tardemr.y=ty;  tardemr.z=tz;

    datatarz=tardem.z;datatarz(datatarz== -9999) = NaN; 
    datarefz=refdem.z;datarefz(datarefz== -9999) = NaN;

    iter= 1;%initialize
    perr=[0.01 0.01 0.01];
    %output:z2out, cpts,dz;p,perr;
    if coregflag==1 % minus rock mean
         %convert nodata values to nan for interpolation
%           refz=data0r.z;refz(refz==-9999)=NaN;%avoid -9999 minus a dem or -9999 minus -9999
          % avoid interpolation for speed
%           demapr = interp2(data0r.x,data0r.y,refz,datar.x,datar.y','*linear',NaN);
%           datatarz=tardem.z;datatarz(datatarz== -9999) = NaN; 
%           datarefz=refdem.z;datarefz(datarefz== -9999) = NaN;
          iter= 1;
          %avoid -9999 - meanhrock
          out{1,1}.z=datatarz;%initial
          [X,Y]=meshgrid(tardem.x,tardem.y);
          out{1,1}.x=X(:);
          out{1,1}.y=Y(:);
          out{1,2}=out{1,1};out{1,2}.z=ones(size(out{1,1}.z));
          out{1,3}=out{1,1};out{1,3}.z=refdem.z;        
          dz = datatarz-datarefz;
          if isempty(rocktag)
              rocktag1=false(size(refdem.z));
          else
              rocktag1=rocktag;
          end
          rockfilter=rocktag1&~isnan(dz)&refdem.z~=-9999&tardem.z~=-9999&abs(dz)<50;%&abs(dz)>1e-4;
          
          if sum(rockfilter(:))<10
            warning(['coregistration failure:',metafile]); p=zeros(3,1);
            iter=49;   
          end
          
          dzrock=dz(rockfilter);
          if isempty(dzrock)
              meanhrock=0; iter= 49;
          else
              meanhrock=mean(dzrock);
          end
          p=[meanhrock,0 0];perr=[0.01 0.01 0.01];
          out{1,1}.z=out{1,1}.z-meanhrock;dz=dz-meanhrock;
          z2out=out{1,1}.z;
          cpts=[X(rockfilter),Y(rockfilter),tardem.z(rockfilter)];

    elseif coregflag==2 % use reg.txt file
          iter= 1;
          out{1,1}.z=tardem.z;%initial
          [X,Y]=meshgrid(tardem.x,tardem.y);
          out{1,1}.x=X(:);
          out{1,1}.y=Y(:);
          out{1,2}=out{1,1};out{1,2}.z=ones(size(out{1,1}.z));
          out{1,3}=out{1,1};out{1,3}.z=refdem.z;
          dzxyt=dzxy{idregion==i};
          if regflag==0 || isempty(dzxyt)
             iter=49;
             dxyz=zeros(1,3);
          else
              dxyztar=[dzxyt(2:3) dzxyt(1)];
              dxyz=dxyztar-dxyzref; % read from reg.txt file
          out{1,4}.P{1,2}(1,1)=1; out{1,4}.P{1,2}(2:4,1)=0;
          out{1,4}.P{1,2}(5:7,1)=-dxyz;
        %in coregisterdems2, tardem+txyz=refdem;in coregisterdems and reg.txt, tardem-dxyz=refdem;
          [outx]=tarx(out,tardem,refdem,[resr resr],params);
          out=outx;
          end
          dz=out{1,1}.z-out{1,3}.z;
          cpts=[];
          
    elseif coregflag==3
        mp1 = rocktag;%interp2(tag.x,tag.y,tag.z,data0r.x,data0r.y','*nearest',0);
        mp2 = rocktag;

        iter= 1;
    %  	[dx,dy,sigma]=vmimc(infile1,infile2) 
%         [z2out,p,sigma] = coregisterdems(refdem.x,refdem.y,double(datarefz),tardem.x,tardem.y,double(datatarz));
%           [z2out,p,d0] = coregisterdems(refdem.x,refdem.y,double(datarefz),tardem.x,tardem.y,double(datatarz),mp1,mp2);
        if isempty(rocktag)||sum(rocktag(:))<10
            [z2out,p,perr,sigma]= coregisterdems(refdem.x,refdem.y,double(datarefz),tardem.x,tardem.y,double(datatarz));
        else
            [z2out,p,perr,sigma]= coregisterdems(refdem.x,refdem.y,double(datarefz),tardem.x,tardem.y,double(datatarz),mp1,mp2);
        end
%         rmsreg2(j)=sigma;
        if sum(size(z2out)~=size(refdem.z)) || sigma>10 || isnan(sigma) || max(abs(p(2:3))) >= 10 %control parameter
            warning(['coregistration failure:',metafile]); p=zeros(3,1);
            iter=49;   
        else
            dz = z2out-refdem.z;
            M=isnan(z2out)|refdem.z==-9999;dz(M)=nan;
            out{1,1}.z=z2out;
            [X,Y]=meshgrid(refdem.x,refdem.y);
            out{1,1}.x=X(:);
            out{1,1}.y=Y(:);
            out{1,2}=out{1,1};out{1,2}.z=ones(size(out{1,1}.z));
            out{1,3}=out{1,1};out{1,3}.z=refdem.z;

            dx=p(2);dy=p(3); %z, x, y
%             pg(j,1:3)=p;
            % 	dzxyt=dzxyd(idreg(iref),:); %dx in the reg.txt file
            % 	if abs(dx)+abs(dy)~=0
            %         dx3=p(2:3);
%             dX4Sg(j,1:3)=p+dx4;%p(2:3)+dx4(2:3); %
            cpts=[X(mp1),Y(mp1),refdem.z(mp1)];

            %         df=p+dx4-dzxyd(j,:); %validation z, x, y
        end
	elseif coregflag==7 %see volcano/code/volcano.m

        [out] = coregisterdems2(tardemr, refdemr,[resrc resrc],params);% best results

        [outx]=tarx(out,tardem, refdem, [resr resr],params);
        out=outx;

%         dX4Sg{j}=out;
        
        % interpolation in coregisterdems2 seems to give bad values between 
        % real height and -9999, producing ripples in height change.
        %filter out edges 
        M=(imdilate(out{1,1}.z== -9999,ones(3)));
        out{1,1}.z(M)=-9999;

        iter=out{1,4}.S{1,1};
        % control points percentage needs to > 0.3%
        [m,n]=size(out{1,2}.z);npt=m*n;
        ncs=out{1,4}.S{1,2}(1); %(npt-sum(out{1,2}.z(:)))
        RC=ncs/npt*100;
        display(['RC=',num2str(RC)])
        %           if RC < 0.0001;iter=49;end
        if ~(RC > 0.3);iter=49;end % SEE Dai and Howat, 2017 Supplementary.

        dz=out{1,1}.z-out{1,3}.z;
        M=out{1,1}.z==-9999|out{1,3}.z==-9999;dz(M)=nan;
        z2out=out{1,1}.z;
        
        %control points; x, y ,z
        M=~out{1,2}.z(:);
        cpts=[out{1,2}.x(M),out{1,2}.y(M), out{1,3}.z(M)];
        %p=-[out{1,4}.P{2}(7), out{1,4}.P{2}(5), out{1,4}.P{2}(6)];%out{1,4}.P{2}(5:7) tx ty tz;pzxy
        p=[0 0 0 ]; %cant just use the 3 parameters only without the rotational and scale parameters.
        perr=[0.01 0.01 0.01];
        
    elseif coregflag==8
        
        odircoregi=[deblank(odircoreg),'/i',num2str(i),'j',num2str(j),'/'];
        if ~exist(odircoregi,'dir')
          mkdir(odircoregi)
        end
        
        if strcmp(metafileref(end-7:end),'meta.txt')
        refimagep=strrep(metafileref,'meta.txt','dem.tif');
        elseif strcmp(metafileref(end-7:end),'_mdf.txt')
        refimagep=strrep(metafileref,'mdf.txt','dem.tif');
        end
        if strcmp(metafile(end-7:end),'meta.txt')
        tarimagep=strrep(metafile,'meta.txt','dem.tif');
        elseif strcmp(metafile(end-7:end),'_mdf.txt')
        tarimagep=strrep(metafile,'mdf.txt','dem.tif');
        end
        str=['time /fs/project/howat.4/SETSM/setsm -coreg 2 -image ',refimagep,' -image ', tarimagep, ' -outpath ', odircoregi];
        fprintf([str,'\n'])
        [status, cmdout]=system(str);
        
        %plot the animation of images and control points before and after coregistration
        %refers to /home/dai.56/arcticdemapp/river/riverwork/coregtest1/plotcontrolpts.m
        filecpt=[odircoregi,'/DEM_gcps_1.txt'];
        if ~exist(filecpt,'file')
           fprintf(['\n setsm DEM coregistration failure id: ',num2str([idem2]),' ',tarimagep,' ',filecpt,'\n'])
           iter=49;   
        else
        cpts=load(filecpt); %x y z
        
        %get coregistration parameters
        %txy=[-1.86 5.04 ];
        coregfile=[odircoregi,'DEM_coreg_result.txt'];
        c=textread(coregfile,'%s','delimiter','\n');   
        [~,name,~] =fileparts([strtrim(tarimagep)]);
        r=find(~cellfun(@isempty,strfind(c,name)));
        %Sigma0[meter]   Tx[meter]       Ty[meter]       Tz[meter]       Tx_std[meter]   Ty_std[meter]   Tz_std[meter]  
        %Dist(mean)      Dist_std(mean)  Dist(med.)      Dist_std(med.)  dh_mean         dh_med.         dh_std         
        %NumberOfCPs     processing time
        c2=c{r};
        r1=strfind(c2,'dem');c2([1:r1(1)+2])='';
        [tmp]=sscanf(c2, '%f',[1,16]);

        txy=-[tmp(2), tmp(3)];   
        p=-[tmp(4) txy]; %zxy
        perr=[tmp(7) tmp(5) tmp(6) ];
        dzmms=tmp(12:14); %dh_mean         dh_med.         dh_std  
        z2out = interp2(tardem.x'-p(2) ,tardem.y-p(3),double(datatarz)-p(1) ,refdem.x',refdem.y,'*linear');

        dz=z2out-refdem.z;
        M=isnan(z2out)|refdem.z==-9999;dz(M)=nan;
        end  % exist
    
	end % if coregflag
    	fprintf(['\n Coregistration method:',num2str(coregflag),'; pzxy=',num2str(p(:)'),'.\n'])

     %check the quality
    if iter==49
%         idd=[idd;j];
%         continue %hi
    end

    %save results to offsets
		
    %write to a reg2.txt file
    if  0
    i=idregion(j);
    infile= strrep([demdir,'/',f{i}],'meta.txt','reg2.txt');
    demfile= strrep([f{i}],'meta.txt','dem.tif');
    fid10 = fopen(infile);
    fprintf(fid10,'DEM Filename: %s \n',demfile);
    fprintf(fid10,'Registration Dataset 1 Name: %s \n',tilefile);
    fprintf(fid10,'Registration Software: coregisterdems  \n');
    fprintf(fid10,'Translation Vector (dz,dx,dy)(m)= %d \n',[dX4S(j,1:3)] );
    fclose(fid10)
    end

    if flagplot==1

        % The case when coregistration is not applied.
        z2n = interp2(tardem.x' ,tardem.y,double(datatarz) ,refdem.x',refdem.y,'*linear');
        nsr=round(resrc/resr);%1; %plot low resolution data
        
        [X,Y]=meshgrid(refdem.x,refdem.y);
        [LATa,LONa]=polarstereo_inv(X,Y,[],[],70,-45);
        if ~isempty(cpts)
            [LATc,LONc]=polarstereo_inv(cpts(:,1),cpts(:,2),[],[],70,-45);
        else
            LATc=[]; LONc=[];
        end

        figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
        % imagesc(data0r.x*1e-3,data0r.y*1e-3,z2n-data0r.z);caxis([-5 5]);colorbar
        surf(LONa(1:nsr:end,1:nsr:end),LATa(1:nsr:end,1:nsr:end),z2n(1:nsr:end,1:nsr:end)-double(datarefz(1:nsr:end,1:nsr:end))); shading interp;
        colorbar;colormap jet;view(0,90)
        hl=xlabel('Longitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        hl=ylabel('Latitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        caxis([-5 5])
        title([texttar,' - ',textref,'; NO coregistration'])
        ofile=[texttar,'m',textref,'Nocoreg'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end

        figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
%         surf(LONa(1:nsr:end,1:nsr:end),LATa(1:nsr:end,1:nsr:end),z2out(1:nsr:end,1:nsr:end)-double(datarefz(1:nsr:end,1:nsr:end))); shading interp;
        surf(LONa(1:nsr:end,1:nsr:end),LATa(1:nsr:end,1:nsr:end),dz(1:nsr:end,1:nsr:end)); shading interp;
        % imagesc(data0r.x*1e-3,data0r.y*1e-3,z2out-data0r.z);caxis([-5 5]);colorbar
        % title('2011/10/08-2013/05/26 DEM (m); After coregistration')
        title([texttar,' - ',textref,'; After coregistration'])
%         plot(X(rockfilter)*1e-3,Y(rockfilter)*1e-3,'k.')
        hold on;plot(LONc,LATc,'k.','Markersize',0.1,'Linewidth',2) %Plot control points 
        colorbar;colormap jet;view(0,90)
        hl=xlabel('Longitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        hl=ylabel('Latitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        caxis([-8 8])
        ofile=[texttar,'m',textref,'wcoreg'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end
        
        %Plot control points %see  scripts/plotsrtmvs.m
        if coregflag==7
        
            if 0 %plot hillshade
        tmpy=reshape(out{1,1}.y,size(out{1,1}.z));tmpx=reshape(out{1,1}.x,size(out{1,1}.z));
        datadsx=tmpx(1,:)';datadsy=tmpy(:,1); % to check the space
        [LAT,LON]=polarstereo_inv(tmpx,tmpy,[],[],70,-45);
        %control points
%         [LATc,LONc]=polarstereo_inv(out{1,2}.x(~out{1,2}.z(:)),out{1,2}.y(~out{1,2}.z(:)),[],[],70,-45);
        [LATc,LONc]=polarstereo_inv(cpts(:,1),cpts(:,2),[],[],70,-45);
        %interpolate to regular lon lat mesh grids for PLOTTING FIGURES
        lonmesh=linspace(min(LON(:)),max(LON(:)),length(datadsx));
        latmesh=linspace(min(LAT(:)),max(LAT(:)),length(datadsy));
        [LONmesh,LATmesh]=meshgrid(lonmesh,latmesh);
        [xm,ym]=polarstereo_fwd(LATmesh,LONmesh,[],[],70,-45);
        tz=double(out{1,3}.z);tz(tz== -9999) = NaN; 
        demrefmesh= interp2(tmpx,tmpy,tz,xm,ym,'linear',nan);
        demrefmesh(isnan(demrefmesh))=-9999;

        % plot control points in reference DEM hillshade
        hills=hillshade(double(out{1,3}.z),datadsx,datadsy,'plotit');
        set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.3 0 3.8 3]); set(gcf, 'PaperSize', [4 3]);
        axis equal
        xlabel('{\itx} (m)')
        ylabel('y (m)')
        % title(['DEM ']);
        hold on
        plot(cpts(:,1),cpts(:,2),'r.','Markersize',2,'Linewidth',2)
        % plot(out{1,1}.x(~out{1,2}.z(:))*1e-3,out{1,1}.y(~out{1,2}.z(:))*1e-3,'r>','Linewidth',2)
        % saveas(gcf,'DEM','fig')
        ofile=[texttar,'',textref,'cpts'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end

        % plot control points in reference DEM hillshade in lon lat
        % hills=hillshade(demrefmesh,lonmesh,latmesh,'plotit','altitude',1,'azimuth',10);
        hills=hillshade(demrefmesh,lonmesh,latmesh,'azimuth',10,'altitude',10,'plotit');
        set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.3 0 3.8 3]); set(gcf, 'PaperSize', [4 3]);
        axis square
        hl=xlabel('Longitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        hl=ylabel('Latitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        hold on
        % plot(LONc,LATc,'r.','Markersize',20,'Linewidth',2)
        plot(LONc,LATc,'r.','Markersize',0.1,'Linewidth',2)
        ofile=[texttar,'',textref,'cptslat'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end
            end %if 0 %plot hillshade

        %control surfaces histogram; dh is distance instead of height difference
        dh=out{1,4}.S{1,4}.height_diff;
        dhaf=out{1,4}.S{1,5}.height_diff;
        figure;
        hold all
        set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.3 0 3.8 3]); set(gcf, 'PaperSize', [4 3]);
        % edges=[-15:1:10];
        edges=[-15:0.1:10];
        % histo1=histogram(dh,edges,'Normalization','pdf');
        histo1=histogram(dh,'Normalization','pdf');
        histo=histogram(dhaf,'Normalization','pdf');
        histo.FaceColor='r';
        alpha(histo,.7)
        legend('Before Coregistration','After Coregistration','Location','NorthWest')
        % legend('After Coregistration','Before Coregistration','Location','NorthWest')
        box on
        ofile=[texttar,'',textref,'pdf'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end

        fprintf ('\n height mean std before coregistration: %f %f\n ',out{1,4}.S{1,4}.height_mean,out{1,4}.S{1,4}.height_std)
        fprintf ('\n height mean std after coregistration: %f %f\n ',out{1,4}.S{1,5}.height_mean,out{1,4}.S{1,5}.height_std)

        % out{1,4}.S{1,4}.height_mean
        % out{1,4}.S{1,4}.height_std
        % out{1,4}.S{1,5}.height_mean
        % out{1,4}.S{1,5}.height_std
        
        end

    end %if plot

% close all
    	
    end % if j %idem2

save test3.mat -v7.3

% remove strips that are not coregistered.
if 0  %
%remove strips that are not coregistered (idd3) and bad quality (idd2).
idd3=find(any(isnan(dX4Sg),2));%find the id of any dZ dX dY that are nan;
idd=unique([idd2(:);idd3(:);]);
idregion(idd)=[];%XYb(idd)=[];dzxy(idd)=[];
% dzxyd(idd,:)=[];%rmsreg(idd)=[];
dX4Sg(idd,:)=[];%rmsreg2(idd)=[];
end
data0r=[];%SAVE MEMORY.

fprintf ('\n Done with coregistration.\n ')
toc
clear data datar 

return
end
