function [Co]=Landslide(tilefile,range,XYbg,f,fdir)
% 
% Requirements: 
% Versions: v2: work on scene dems, since strip DEM are heavily filtered.
% April 2018
% Modification (November 2019): 
%  1\ using adjustOffsets.m for readjusting coregistration offsets.
%  2\ using coreg3.m for coregistration;
%  3\ do not apply the grouping of DEM files based on overlapping zones.
% 

% Chunli Dai chunlidai1@gmail.com
% April 2018
% 

% close all
% clear 
% clc

constant

params.I = [1 0 0 0 0 0 0];
params.C = [50 0.001 0.001 0.05 0.0001];
%         params.G = [3000 20];
params.M = '3D';
params.V = [10 20 10];
params.G = [9000 20]; %

width=100;%2e3;%7e3; %buffer width of the a priori coastline, e.g., 2km.

%flagplot=0;

res='2';
% demdir=dir(['/data2/ArcticDEM/region_',regionnum,'*']);
macdir='/Users/chunlidai/surge/';
% macdir='/Users/chunli/surge/';
macdir=[];
%addpath(genpath([macdir,'/data/chunli/scripts/']));
%addpath(genpath([macdir,'/data/chunli/coastline/']));
%addpath(genpath(pwd));

if 0 %for testing
tilefile=[macdir,'/data/chunli/coastline/orthorect1/54_06_2_2_5m_v2.0/54_06_2_2_5m_v2.0_reg_dem.tif'];
tilefile='44_08_1_2_5m_v2.0_reg_dem.tif'; %'44_08_1_1_5m_v2.0' 20151017Landslide
tilefile='21_38_1_1_5m_v2.0_reg_dem.tif'; % Nuugaatsiaq Greenland landslide and tsunami 2017 Jun 17
end

yr=365.25;
neq=1;
dsr=0.2;%0.04; %200m ; %0.2; % 8m to 40m
nsr=1./dsr;
%resr=2; %4;%40.; %2 %str2double(res)/dsr;
% resr=40.; 
resrc=40.; %for coregisteration
% res=mt.info.map_info.dx;
% demdir=[macdir,'/data2/ArcticDEM/region_08_canada_baffin/tif_results/8m/'];
%coregflag=3; %3;%1 parameter (vertical), 3 parameter (Ian's); 7 parameter (MJ's)

filename='boundaries_regall_strip.dat'; %'boundaries_reg31.dat';
filename='boundaries_reg31.dat';
% filename='boundaries_reg03_2m.dat'; %'boundaries_reg31.dat';
Co=[];

%Grids of ArcticDEM Tiles
% 54_06_2_2_5m_v2 yid_xid_xids_yids
% 1,2 ; 2,2
% 1,1 ; 2,1
%x=x0+(xid-1)*dx; %left edge of a box;
dx=100e3;x0=-4000e3;y0=-4000e3;%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;
%dxs=dx/2/10; %5km
dxs=dx/50; %2km
% %%% Preparation: get rang0 from tile name
if isa(tilefile,'char')
    [dir,ifile,ext] =fileparts(tilefile);
    r=1;
    xid= sscanf(ifile(r+3:(r+4)), '%g', 1);
    yid= sscanf(ifile(r:(r+1)), '%g', 1);
    xids= sscanf(ifile(r+6), '%g', 1);
    yids= sscanf(ifile(r+8), '%g', 1);
    xidss= sscanf(ifile(r+10:r+11), '%g', 1);
    yidss= sscanf(ifile(r+13:r+14), '%g', 1);
    % [xid,yid,xids,yids]
    % xid=14;yid=51;xids=2;yids=1;%51_14_2_1_5m_v2.0
    % x=x0+(xid-1)*dx+(xids-1)*dx/2;y=y0+(yid-1)*dx+(yids-1)*dx/2;
    % rang0b=[x x+dx/2 y y+dx/2]; %exact tile boundary
    % rang0=[x-width x+dx/2+width y-width y+dx/2+width]; %tile boundary with buffer width
    x=x0+(xid-1)*dx+(xids-1)*dx/2+(xidss-1)*dxs;y=y0+(yid-1)*dx+(yids-1)*dx/2+(yidss-1)*dxs;
    rang0=[x-width x+dxs+width y-width y+dxs+width]; %tile boundary with buffer width
    odir=[tilefile(1:15)];
else %given rang0
    rang0=tilefile;
    odir=['givensite']; %notice odir should be all characters;
end
    
% odir=['./',tilefile(1:15),'/'];
%if ~exist(odir,'dir')
%    mkdir(odir)
%    cd(odir)
%end

% rang0 = [  -3450000    -3400000     1350000     1400000];
if 0 %for testing
    rang0 = [-3241005,-3112012,207228,381469];rang0b=rang0;
end

xeq=(rang0(1)+rang0(2))/2;yeq=(rang0(3)+rang0(4))/2;
[lateq,loneq]=polarstereo_inv(xeq,yeq,[], [],70,-45);

if 0
loneq=-85.9418555555556; lateq=	79.988633; %Eureka
[xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
exb=5e3;rang0=[xeq-exb xeq+exb yeq-exb yeq+exb ]; %for test only

rang0=[-711.4e3-500 -710.4e3+500 -821.9e3-500 -820.9e3+500]; %2km by 2km
%rang0=[-711.4e3     -710.4e3     -821.9e3     -820.9e3    ]; %1km by 1km
xeq=(rang0(1)+rang0(2))/2;yeq=(rang0(3)+rang0(4))/2;
[lateq,loneq]=polarstereo_inv(xeq,yeq,[], [],70,-45);
end

if 0
lateq=61.093035; loneq=-140.310993; % landslide Yukon 2007 July 24
[xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
teq=datenum('20070724','yyyymmdd');
texteq={'20070724Landslide'};
eqepoch=teq;timefix=1; % 1 assume event time known; 0 assume time unknown
% eqepoch=0;timefix=0;
end

if 0
%Nuugaatsiaq Greenland landslide and tsunami 2017 Jun 17
%https://ds.iris.edu/ds/nodes/dmc/specialevents/2017/06/22/nuugaatsiaq-greenland-landslide-and-tsunami/
lateq=71.640; loneq=-52.344;
[xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
teq=datenum('20170617','yyyymmdd');
texteq={'20170617Landslide'};
eqepoch=teq;timefix=1; % 1 assume event time known; 0 assume time unknown
%eqepoch=0;timefix=0;
exb=6e3;rang0=[xeq-exb xeq+exb yeq-exb yeq+exb ];
end


%Plot for time series
if 0 %Tyndall
latq1=6.0178293e+01; lonq1= -1.4118490e+02; % Point for time series. Oct. 17 2015 , largest landslide %wrong lat lon R e
%latq1=lateq;lonq1=loneq;
[xq1,yq1]=polarstereo_fwd(latq1,lonq1,[],[],70,-45);

xeq=xq1;yeq=yq1;
[lateq,loneq]=polarstereo_inv(xeq,yeq,[], [],70,-45);
exb=6e3;
%exb=3.4e3;
rang0=[xeq-exb xeq+exb yeq-exb yeq+exb ]; %for test only

xq1(1)=xeq       ; yq1(1)=yeq       ; %-711.186      -820.982; time 2011;dist 20m
xq1(2)=-3283462.34; yq1(2)=355822.63;
xq1(3)=-3302e3    ; yq1(3)=348.5e3  ;
end

if 0
xq1=-3283462.34;yq1=355822.63; % Point for time series. Oct. 17 2015 , largest landslide
xq1=-3302e3;yq1=348.5e3; % -19.56 std
xq1= -255800;yq1=   -1988720;%greenland landslide 2017
[latq1,lonq1]=polarstereo_inv(xq1,yq1,[], [],70,-45); %[60.1768700940536 -141.184899979958]
end


if 0 %Point for time series
latq1=lateq;lonq1=loneq; 
[xq1,yq1]=polarstereo_fwd(latq1,lonq1,[],[],70,-45);
xq1(1)=-711.186e3; yq1(1)=-820.982e3; %-711.186      -820.982; time 2011;dist 20m
xq1(2)=-711.178e3  ;yq1(2)=-820.984e3 ; %-711.178      -820.984; time 2013; dist 30m
xq1(3)=-711.146e3 ; yq1(3)=-820.992e3;% -711.146      -820.992 ; time 2015; dist 58m
xq1(4)=-711.134e3 ; yq1(4)=-820.996e3; %-711.134      -820.996; time 2017; dist 74m
xq1(5)=-711.118e3 ; yq1(5)=-821.004e3; %-711.118      -821.004;time 2011; dist 94m
xq1(6)=-711.098e3; yq1(6)=-821.018e3;%-711.098      -821.018;time 2013; dist 113m
xq1(7)=-711.074e3 ; yq1(7)=-821.04e3 ; %-711.074       -821.04; time 2014; dist 150m
xq1(8)=-711.066e3 ; yq1(8)=-821.054e3; %-711.066      -821.054; time 2016; dist 168m
xq1(9)=-711.06e3  ; yq1(9)=-821.064e3;%-711.06      -821.064; time 2017;dist 176m
xq1(10)=-710607.9712; yq1(10)=-821210.022; %-710607.9712      -821210.022
xq1(11)=-711419.9829; yq1(11)=-821502.0142;

nq1=length(xq1);
end

%Read files including all points (lon lat) for time series;
%For a 5 km block, check four points on regular grids (at e.g. +(1.25, 1.25), (3.75, 1.25), (1.25, 3.75), (3.75, 3.75) km pixels ).
dxq1=abs(rang0(2)-rang0(1))/2;dyq1=abs(rang0(4)-rang0(3))/2;
xq1grid=[rang0(1)+dxq1/2; rang0(1)+dxq1/2*3; rang0(1)+dxq1/2; rang0(1)+dxq1/2*3;];
yq1grid=[rang0(3)+dyq1/2; rang0(3)+dyq1/2; rang0(3)+dyq1/2*3; rang0(3)+dyq1/2*3;];
%[latq1grid,lonq1grid]=polarstereo_inv(xq1grid,yq1grid,[], [],70,-45); %[60.1768700940536 -141.184899979958]
filenameq1='xyq1.dat'; %x y
if exist(filenameq1,'file')
fidq1 = fopen(filenameq1);
nq1 = linecount(fidq1);
fidq1 = fopen(filenameq1);
xyq1=fscanf(fidq1, '%f', [2, nq1])';
fclose(fidq1);
xq1=[xyq1(:,1);xq1grid(:);];
yq1=[xyq1(:,2);yq1grid(:);];
else
xq1=[xq1grid(:);];
yq1=[yq1grid(:);];
end
[latq1,lonq1]=polarstereo_inv(xq1,yq1,[], [],70,-45);
nq1=length(xq1);

fprintf(['\n Total of ',num2str(nq1),' points to check time series. \n'])
fprintf(['\n Input points: '])
for i=1:nq1-length(xq1grid)
    fprintf(['\n ',num2str(i), ' (lon lat x y): ',num2str([lonq1(i), latq1(i), xq1(i), yq1(i)])])
end
fprintf(['\n\n Regular grid points: '])
for i=nq1-length(xq1grid)+1:nq1
    fprintf(['\n ',num2str(i), ' (lon lat x y): ',num2str([lonq1(i), latq1(i), xq1(i), yq1(i)])])
end

% xq1=-3283e3;yq1=355e3; % 2 glacier
% xq1=-3285e3;yq1=354.4e3; % 3 river
% xq1=-3283e3;yq1=359.8e3; %4 void
% xq1=-3284.3e3;yq1=362e3; % 5

% rang0=[-3467 -3455 110 124 ]*1e3;
x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
[lat0,lon0]=polarstereo_inv(x0,y0,[], [],70,-45);
ranget=round(rang0/resr)*resr;rang0=ranget;
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
xout=tx;yout=ty;
fprintf(['\n rang0 resr= ',num2str([rang0, resr]),'.\n'])

if 0 %move to main program to save time
fprintf ('\n Step 0: geting the boundary for all files in the region.\n')
%READ INPUT PARAMETERS; getting the boundaries for all files
% filename='boundaries_reg31.dat';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
% range=fscanf(fid, '%f', [4, n))';
range=zeros(n,4);
for i=1:n
   range(i,1:4)=fscanf(fid, '%f', [4, 1])';ifile=fgetl(fid);
   [demdir,name,ext] =fileparts([macdir,strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
end
display(['demdir=',demdir])
end

x=[range(:,1) range(:,2) range(:,2) range(:,1) range(:,1) ];y=[range(:,4) range(:,4) range(:,3) range(:,3) range(:,4) ];
id=find(range(:,2)>rang0(1) & range(:,1)<rang0(2) & range(:,4)>rang0(3) & range(:,3)<rang0(4));
% id=1:length(range(:,1));

%select measurements in desired time span
mon=zeros(length(id),1);
year=zeros(length(id),1);
t1=zeros(length(id),1);
Mflag=true(size(mon)); %1 use data; 0 do not use this data; depends on the time span of the cross-track images.
Mflagc=false(size(mon)); %1 cross track; 0 in-track
for j=1:length(id)
%       ymd=filename2ymd(f{id(j)});
%       mon(j)=str2num(ymd(5:6));

        infile= [fdir{id(j)},'/',f{id(j)}];
        [ymd,flago]=strip2date(infile,3);
        mon(j)=str2num(ymd(5:6));
        year(j)=str2num(ymd(1:4));
        t1(j)=datenum(ymd,'yyyymmddHHMMSS');
	if flago==0; Mflag(j)=0;end

	% flag/skip cross-track strips W1W1
%       satname=filename2sat(infile);
%       expression='W[0-9]';
%       satstr=regexp(satname,expression,'match');
%       if ~isempty(satstr) %match W1W2; cross track
%            Mflagc(j)=1;
%       end
	expression='WV';infile1=f{id(j)};
	satstr=regexp(infile1,expression,'match');
	if ~isempty(satstr) %match in-track;
	    %do nothing
	else %cross track
	    Mflagc(j)=1;
	end


	% delete given dates of DEMs in constant.m, April 2024
        for mi=1:length(deldate)
        if (abs(deldate(mi)-t1(j))< 1)
            Mflag(j)=0;
	    fprintf(['\n Landslide.m delete data: ',infile,'\n '])
        end
        end

end
%idd=~((mon>=mons&mon<=mone)&Mflag&(year>=year_start&year<=year_end));
idd=~((mon>=mons&mon<=mone)&Mflag&(t1>=date_start&t1<=date_end));
%save testyears2.mat -v7.3
id(idd)=[];
mon(idd)=[]; Mflag(idd)=[]; year(idd)=[];

%reduce data repeats
if length(id)>=200
   fprintf (['\n Data is too many:',num2str(length(id)),' repeats; Only keep summer data and In-track data. \n'])
   mons=6;mone=9;
   %idd=~((mon>=mons&mon<=mone)&Mflag&(year>=year_start&year<=year_end)&Mflagc==0);
   idd=~((mon>=mons&mon<=mone)&Mflag&(t1>=date_start&t1<=date_end)&Mflagc==0);
   id(idd)=[];
   fprintf (['\n Summer data count:',num2str(length(id)),' repeats. \n'])
end

if isempty(id);
fprintf (['\n No data in the selected region!\n'])
return;
end

fprintf ('\n Step 1: Preparing the grouping of common overlapping pieces.')
fprintf ('\n Step 1.1: getting the real Polygon boundary for all files over the output zone.')
%load all reg.txt and meta.txt files in the coverage
idregion=id;XYb=cell(size(idregion));dzxy=cell(size(idregion));count=0;
rmsreg=zeros(size(idregion));idd=[];dzxyd=zeros(length(idregion),3);
flagcb=zeros(size(idregion));
for j=1:length(idregion)
    i=idregion(j);
    XYbi=XYbg{i};
    XYb{j}=XYbi;
end
for j=1:0 %length(idregion)
    i=idregion(j);
 % get the Polygon boundary for actual data
        demdir=fdir{i};
	infile=[demdir,'/',f{i}];
        if strcmp(infile(end-7:end),'meta.txt')
          str1='meta.txt';
        elseif strcmp(infile(end-7:end),'_mdf.txt')
          str1='mdf.txt';
        end

        infile= strrep([demdir,'/',f{i}],str1,str1);
        [XYbi,rangei]=imagebd(infile);
        XYb{j}=XYbi;
        
 % get the reg.txt
        infile= strrep([demdir,'/',f{i}],str1,'reg.txt');
        if exist(infile,'file')
        c=textread(infile,'%s','delimiter','\n');
        r=find(~cellfun(@isempty,strfind(c,'Translation')));
%       dzxys=deblank(strrep(c{r(1)},'Translation vector (dz,dx,dy)(meters)= ',''));
        c1=c{r(1)};r1=strfind(c1,'=');c1(1:r1)='';dzxys=c1;
        dzxy1=textscan(dzxys,'%f %f %f','Delimiter',',');
        dzxy{j}=cell2mat(dzxy1);
        dzxyd(j,1:3)=cell2mat(dzxy1);
        %Mean Vertical Residual
        r=find(~cellfun(@isempty,strfind(c,'Mean Vertical Residual')));
        c1=c{r(1)};r1=strfind(c1,'=');c1(1:r1)='';
        rmst=sscanf(c1, '%g', 1);
        rmsreg(j)=rmst;
        else
            count=count+1;idd=[idd;j];
            display([num2str(count),' Reg file do not exist:',infile])
        end
end
display([num2str(count),' out of ',num2str(length(idregion)),' Reg file do not exist!'])
idd=find(rmsreg>10|abs(dzxyd(:,2))>10|abs(dzxyd(:,3))>10);
dzxyd(idd,:)=0;
%End of loading
if 0 %if use scene tiff files rather than strip files.
idregion=id;XYb=cell(size(idregion));
for j=1:length(idregion)
    i=idregion(j);
    XYb{j}=[x(i,:)',y(i,:)'];
end
end

%Plot the coverage
if flagplot==1
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
plot(xeq,yeq,'r*','Markersize',12)
hold on;plot(rang0([1 2 2 1 1]), rang0([3 3 4 4 3]), 'r-')
% plot(x(id), y(id) ,'.-')
for i=1:length(id)
    plot(x(id(i),:), y(id(i),:),'b>-','linewidth',4)
%     hold off
title([f{id(i)}(1:13)])
% pause
end
hold on

figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
plot(xeq,yeq,'r*','Markersize',12)
for i=1:length(id)
    Xb=XYb{i}(:,1);
    Yb=XYb{i}(:,2);
    plot(Xb,Yb,'b>-','linewidth',4)
title([f{id(i)}(1:13)])
% pause
end
hold on;plot(rang0([1 2 2 1 1]), rang0([3 3 4 4 3]), 'r-')
saveas(gcf,'polygons','fig')


[lat,lon]=polarstereo_inv(x,y,[], [],70,-45);
% lon(lon>=0)=lon(lon>=0)-360;
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
% plot(lon', lat' ,'-')
for i=1:length(id)
    plot(lon(id(i),:), lat(id(i),:),'b>-','linewidth',4)
%     hold off
title([f{id(i)}(1:13)])
end
plot(loneq, lateq ,'r*','Markersize',12)
end %if flagplot==1

fprintf ('\n Step 1.2: Searching for overlappings of Polygons at regular big grids.')
if flagplot==1
tic
%sub-zones for searching overlappings
dx=500;%400;%3600;%7e3;%%17e3; %3600; %80; %3600; 
rang0t=round(rang0);
rangx=rang0t(1):dx:rang0t(2);nx=length(rangx)-1;
rangy=rang0t(3):dx:rang0t(4);ny=length(rangy)-1;
ns=nx*ny;rang0s=zeros(ns,4);

novlp=zeros(ny,nx);baovlp=zeros(ny,nx);idg=cell(ns,1);%idg{ns}=[];
tspan=zeros(ny,nx);
enl=0;%;0.3;
for ix=1:nx
    for iy=1:ny
        ixy=iy+(ix-1)*ny;
        rang0s(ixy,:)=[rangx(ix) rangx(ix+1) rangy(iy) rangy(iy+1) ];
        id=find(range(:,1)<=rangx(ix)-enl*dx & range(:,2)>=rangx(ix+1)+enl*dx & range(:,3)<=rangy(iy)-enl*dx  & range(:,4)>=rangy(iy+1)+enl*dx);
%         novlp(iy,ix)=length(id);
        % Refinement of the overlapping count, checking the actual data
        % coverage overlapping with the subzone.
        x0si=[rang0s(ixy,1)-enl*dx rang0s(ixy,2)+enl*dx rang0s(ixy,2)+enl*dx rang0s(ixy,1)-enl*dx rang0s(ixy,1)-enl*dx ];
        y0si=[rang0s(ixy,4)+enl*dx rang0s(ixy,4)+enl*dx rang0s(ixy,3)-enl*dx rang0s(ixy,3)-enl*dx rang0s(ixy,4)+enl*dx ];
        idd=[];str=cell(length(id),1); 
        for j=1:length(id)
            % get the Polygon boundary for actual data
            i=id(j); 
          % str{j}=f{id(j)}(1:13);
	    sat=filename2sat(f{id(j)});
            ymd=filename2ymd(f{id(j)});
            str{j}=[sat,'_',ymd];

            M=idregion==i;nt=sum(M); % only work on the strips that are selected (cover coastline and successfully coregistered).
            if nt==0 ;idd=[idd;j];continue;end
            Xb=XYb{M}(:,1);
            Yb=XYb{M}(:,2);

            % new method
            in = inpolygon(x0si,y0si,Xb,Yb); %whether subtiles are inside the polygon
            if any(in==0) %any subtile corners not in the boundary
               idd=[idd;j];
            end
        end % j=1:length(id)
        id(idd)=[];str(idd)=[];
            
        %get rid of the strips have the same date and same sensor
	if 0
        [un idx_last idx] = unique(str(:));
        id1=1:length(id);idd=id1(~ismember(id1,idx_last));
        id(idd)=[];
	end
        
        %find the maximum time span
        t=zeros(length(id),1);
        for j=1:length(id)
        %ymd=f{id(j)}(6:13);i=id(j); 
	filename=f{id(j)};
%        ymd=filename2ymd(filename);
        infile= [fdir{id(j)},'/',f{id(j)}];
        [ymd,~]=strip2date(infile,3);
	i=id(j); 
        %t(j)=datenum(ymd,'yyyymmdd');
        t(j)=datenum(ymd,'yyyymmddHHMMSS');
        end
        epoch=t/yr;  
        dt=max(epoch)-min(epoch); % in year % maybe for later version 
        if isempty(dt); dt=0;end
        tspan(iy,ix)=dt;

        %if sum(t>teq)&&sum(t<teq); baovlp(iy,ix)=1;end % if data exist both before and after earthquake;

        novlp(iy,ix)=length(id);
        idg{ixy}=sort(id);   %finding the DEMs at each zones.    
    end
end
display(['Counting overlapping...']);
toc

% [X,Y]=meshgrid(rangx(1:end-1),rangy(1:end-1));
[X,Y]=meshgrid((rangx(1:end-1)+rangx(2:end))/2,(rangy(1:end-1)+rangy(2:end))/2);
Xun=X;Yun=Y; %save X Y for idun map

end %if flagplot==1

if flagplot==1
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
surf(X*1e-3,Y*1e-3,novlp);colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot3(xeq*1e-3, yeq*1e-3,1e2 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
saveas(gcf,'OverlappingCount_Alaska','fig')

[LAT,LON]=polarstereo_inv(X,Y,[], [],70,-45);

figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
surf(LON,LAT,novlp);colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot3(loneq, lateq,1e2 ,'r*','Markersize',12)
saveas(gcf,'OverlappingCount','fig')
end

if 0
id=X>=rang0b(1)&X<=rang0b(2)&Y>=rang0b(3)&Y<=rang0b(4);
out=[LAT(id) LON(id) novlp(id) tspan(id)];
% save -ascii latlonnov_barnesstrip8m.dat out
fid = fopen('latlonnov.dat','a'); %statistics
fprintf(fid,' %f %f %f %f \n',out');
fclose(fid);
end

% return

% Initialize
nsuby=length(yout);nsubx=length(xout);
fprintf(['\n nsuby nsubx = ',num2str([nsuby, nsubx]),'.\n'])

% fid2 = fopen('trend.txt', 'w');
% fid4 = fopen('demTyndalbig.dat','w');
% fid3= fopen('trenddf.txt', 'a');

clear lenfc
% novlpf=zeros(nsuby,nsubx,'int32');%maximum repeat of the chosen group at a pixel
% jump=zeros(nsuby,nsubx);jumpstd=zeros(nsuby,nsubx);
% timec=zeros(nsuby,nsubx);timestd=zeros(nsuby,nsubx);
% oflagc=zeros(nsuby,nsubx,'int32');

% [XOUT,YOUT]=meshgrid(xout,yout);

fprintf ('\n Step 2: Coregistration and time series analysis.')
%get a priori stable surfaces e.g. rock data;
%tagfile='/home/chunli/scripts/BarnesrockRGI.mat';
if flagrock ==1
	%tagfile=[pwd,'/BarnesrockRGI.mat']; %mannually retrieved mask
	tagfile=[pwd,'/rockmask.mat']; %mannually retrieved mask
	if exist(tagfile,'file')
	   load(tagfile); %tag in logical
	   fprintf(['\n Use existing rock mask file: ',tagfile,'. \n'])
	else
	   if flaggreenland==1 %1 use GIMP masks for greenland; 0 otherwise
	   [tag]=getrockicegreenland(rang0);
	   else
	   [tag]=getrgi(rang0); 
	   end
	end
else
	tag=[];
end

listfile=['demlistprecoreg.txt'];
fidlist = fopen(listfile, 'w');
for j=1:length(idregion) 
        i=idregion(j);
        demdir=fdir{i};
        infile= [demdir,'/',f{i}];
        fprintf(fidlist,'%s \n',infile);
end %j
fclose(fidlist)

% %%%% coregister all pairs of DEM strips.
listfile=[odir,'_listused.txt'];
if ~exist('dX4Sg.mat','file')
% [dX4Sg,idregion,data0r]=coreg3([],idregion,XYbg,f,fdir,tag); %load the whole strip of data and coregister.
[dX4Sg,idregion,data0r]=coreg3(rang0,idregion,XYbg,f,fdir,tag); %only load data within rang0;
%save dX4Sg.mat dX4Sg idregion
save dX4Sg.mat dX4Sg idregion XYbg f fdir -v7.3
%j=1; i=idregion(j); f(i) corresponds to dX4Sg(j)
else
load dX4Sg.mat %hi
[n1,m1]=size(dX4Sg);
  if n1~=length(idregion)
	warning('Landslide.m: loaded dX4Sg not matching the ids!')
  end
end

%write the list of files used for change detection

    %sort idregion based on time
    idblock=1:length(idregion);
    nid=length(idregion);
    t=zeros(nid,1); 
    for j=1:nid 
	i=idregion(j);
        filename=f{i};
        infile= [fdir{i},'/',f{i}];
        [ymd,~]=strip2date(infile,3);
        t(j)=datenum(ymd,'yyyymmddHHMMSS');
    end
    [~,idsort]=sort(t);idblock=idblock(idsort);

listfile=[odir,'_listused.txt'];
fidlist = fopen(listfile, 'w');
fprintf(fidlist,'%s \n','# px     py   pz   DEMFileName');
fprintf(fidlist,'%s \n','# Strip DEM - [px py pz] =Coordinates aligned with the reference DEM (reference in bundle adjustment).');
for j=1:length(idregion) %[3,5,24,26] % 2012 Kamchatka Volcano 
%        i=idregion(j);
	i=idregion(idblock(j));
        demdir=fdir{i};
        infile= [demdir,'/',f{i}];
        p=dX4Sg(idblock(j),:) ; %pzxy
	%fprintf(fidlist,'%s \n',infile);
        px=p(2);py=p(3);pz=p(1);
        fprintf(fidlist,'%12.2f %12.2f %12.2f  %s \n',px, py, pz, infile);
%fprintf(['\n Loading DEM files in chronological order (',num2str(j),'/',num2str(nid),', isel ',num2str(isel),'): ',infile])
end %j
fclose(fidlist)

if flag_diffdems==1
[Co] = diffdems_sub(dX4Sg, idregion, XYbg, f, fdir, rang0);
end

%Divide the study zone to 1 km by 1 km blocks; Computation time 15 minutes per block.
tcpu1 = cputime;
ck0=clock;
blocksize=round(1000./resr)+1; % +1 to avoid a block with size of 1 pixel.
%blocksize=round(1000./resr); 
nblocky=length(1:blocksize:nsuby);
nblockx=length(1:blocksize:nsubx);
npc=nblocky*nblockx;
fprintf(['\n Target zone is divided to ',num2str([nblocky]),' by ',num2str([nblockx]),' blocks.\n'])
%initialize cells for output variables
novlpf_cell=cell(npc);%maximum repeat of the chosen group at a pixel
jump_cell=cell(npc);jumpstd_cell=cell(npc);
rate_cell=cell(npc);ratestd_cell=cell(npc);
timec_cell=cell(npc);timestd_cell=cell(npc);
oflagc_cell=cell(npc);
%poolobj=parpool(poolsize);
%parfor isel=1:npc %[71,85] %1:npc  %block %must be consecutive increasing integers
for isel=1:npc %[71,85] %1:npc  %block %must be consecutive increasing integers

    [jumpt2,jumpstdt2,ratet2,ratestdt2,timec2,timestd2,novlpf2,oflagc2]=changedetectionblock(isel,blocksize,xout,yout,dX4Sg,idregion,XYbg,f,fdir,xq1,yq1);

    novlpf_cell{isel}=novlpf2;%maximum repeat of the chosen group at a pixel
    jump_cell{isel}=jumpt2;jumpstd_cell{isel}=jumpstdt2;
    rate_cell{isel}=ratet2;ratestd_cell{isel}=ratestdt2;
    timec_cell{isel}=timec2;timestd_cell{isel}=timestd2;
    oflagc_cell{isel}=oflagc2;
    
end %im   or isel  
%delete(poolobj)

%assign all data to one big matrix;
%initialize
jump=zeros(nsuby,nsubx);jumpstd=zeros(nsuby,nsubx);
rate=zeros(nsuby,nsubx);ratestd=zeros(nsuby,nsubx);
timec=zeros(nsuby,nsubx);
if algorithmin==1 %ice melting
timestd=zeros(nsuby,nsubx,2);
time1=zeros(nsuby,nsubx);
time2=zeros(nsuby,nsubx);
else
timestd=zeros(nsuby,nsubx);
end
novlpf=zeros(nsuby,nsubx,'int32');%maximum repeat of the chosen group at a pixel
oflagc=zeros(nsuby,nsubx,'int32');
iselop=zeros(nsuby,nsubx,'int32'); %isel id of the chosen group/each block
for isel=1:npc
    jumpt2=jump_cell{isel};jumpstdt2=jumpstd_cell{isel};
    ratet2=rate_cell{isel};ratestdt2=ratestd_cell{isel};
    timec2=timec_cell{isel};timestd2=timestd_cell{isel};
    novlpf2=novlpf_cell{isel};oflagc2=oflagc_cell{isel};
    
    %assgin to big matrix
    iblock=ceil(isel/nblocky);jblock=isel-(iblock-1)*nblocky;
    %mx0 my0 index of xout yout within this block.
    mx0=(iblock-1)*blocksize+1:min([nsubx,iblock*blocksize]); 
    my0=(jblock-1)*blocksize+1:min([nsuby,jblock*blocksize]); 
    jump(my0,mx0)=jumpt2;jumpstd(my0,mx0)=jumpstdt2;
    rate(my0,mx0)=ratet2;ratestd(my0,mx0)=ratestdt2;
    timec(my0,mx0)=timec2;
    if algorithmin==1 %ice melting
       time1(my0,mx0)=timestd2(:,:,1);
       time2(my0,mx0)=timestd2(:,:,2);
       timestd(my0,mx0)=timestd2(:,:,2)-timestd2(:,:,1);
    else
    timestd(my0,mx0)=timestd2;
    end
    novlpf(my0,mx0)=novlpf2;oflagc(my0,mx0)=oflagc2;
    iselop(my0,mx0)=isel;
    
end %isel
%start and end time
if algorithmin~=1 %ice melting
time1=timec-timestd;time2=timec+timestd;
end

e = cputime-tcpu1
ck2=clock;
dt=ck2-ck0;tsec=dt(3)*86400+dt(4)*3600+dt(5)*60+dt(6);
fprintf (['\n Change estimation process takes: ',num2str(tsec),' sec'])
    
save sv1.mat xout yout oflagc novlpf iselop timestd timec jumpstd jump rate ratestd time1 time2 -v7.3
clear datarsv data  datamtrsv

%% write output
projstr='polar stereo north';
if algorithmin==2 ||algorithmin==3
OutName=[odir,'_jump.tif'];
%writeGeotiff(OutName,xout,yout,double(jump),4,0,projstr)
jump(jump==0)=nan; % 0 to nan
writeGeotiff(OutName,xout,yout,double(jump),4,nan,projstr); 

OutName=[odir,'_jumpstd.tif'];
%writeGeotiff(OutName,xout,yout,double(jumpstd),4,0,projstr)
jumpstd(jumpstd==0)=nan; % 0 to nan
writeGeotiff(OutName,xout,yout,double(jumpstd),4,nan,projstr); 
end

if algorithmin==1 ||algorithmin==3
OutName=[odir,'_rate.tif'];
%writeGeotiff(OutName,xout,yout,double(rate),4,0,projstr)
rate(rate==0)=nan; % 0 to nan
writeGeotiff(OutName,xout,yout,double(rate),4,nan,projstr); 

OutName=[odir,'_ratestd.tif'];
%writeGeotiff(OutName,xout,yout,double(ratestd),4,0,projstr)
ratestd(ratestd==0)=nan; % 0 to nan
writeGeotiff(OutName,xout,yout,double(ratestd),4,nan,projstr); 
end

if 0
OutName=[odir,'_eventtime.tif']; %sugggest to convert to julian date, and get year.doy/365.25
writeGeotiff(OutName,xout,yout,double(timec),4,0,projstr)

OutName=[odir,'_eventtimestd.tif']; %unit years
writeGeotiff(OutName,xout,yout,double(timestd),4,0,projstr)
end

%date format YYYYMMDD for T1 T2, where T1 is the closest date of measurement before change
%	T2 is the closest date of measurement after change.
OutName=[odir,'_eventtimeT1.tif']; %
timet1=str2num(datestr(time1*yr,'yyyymmdd')); %3; fmtstr = 'int32';
timet1=reshape(timet1,size(time1));
writeGeotiff(OutName,xout,yout,int32(timet1),3,0,projstr)

timet2=str2num(datestr(time2*yr,'yyyymmdd'));
timet2=reshape(timet2,size(time2));
OutName=[odir,'_eventtimeT2.tif']; %unit years
writeGeotiff(OutName,xout,yout,int32(timet2),3,0,projstr)
if algorithmin==2 ||algorithmin==3
elseif algorithmin==1
%Ian asks for average time for all measurements used for ice melting 

OutName=[odir,'_eventtimeTave.tif']; %
if 1 %format YYYYMMDD
timetc=str2num(datestr(timec*yr,'yyyymmdd')); %3; fmtstr = 'int32';
timetc=reshape(timetc,size(timec));
writeGeotiff(OutName,xout,yout,int32(timetc),3,0,projstr)
end
%format year.
if 0
timec(timec==0)=nan; % 0 to nan
writeGeotiff(OutName,xout,yout,double(timec),4,nan,projstr)
end % if 0

end

OutName=[odir,'_nov.tif']; %number of overlapping DEMs
writeGeotiff(OutName,xout,yout,int32(novlpf),3,0,projstr)

%algorithmin=2;%fit model for time series. 
        %e.g.  1 linear (ice melting); 2 constant (landslides) ; 3 constant + linear (Okmok volcano)


%output the time in date strings as requested by Melissa
if 0
%start and end time
%time1=zeros(nsuby,nsubx);time2=zeros(nsuby,nsubx);
fid = fopen('latlontime.dat','w'); %statistics
%fprintf(fid,' %f %f %f %f \n',out');
[X,Y]=meshgrid(xout,yout);
[LAT,LON]=polarstereo_inv(X,Y,[], [],70,-45);
%slow: 1 hour for 5km by 5 km area with 2 m resolution.
%to do: matrix computation and print
for i=1:length(xout)
for j=1:length(yout)
	eqm=timec(j,i)*yr; % year to days
	eqstd=timestd(j,i)*yr;
	%eqmc=datestr(eqm*yr,26);
	%eqmc=datestr(eqm,'yyyy-mm-ddTHH:MM:SS');
	eqmc=datestr(eqm,26);
	eqsc=datestr(eqm-eqstd,26);
	eqec=datestr(eqm+eqstd,26);
	%fprintf(fid,'%12.6f %12.6f  %s %s %s\n',LAT(j,i),LON(j,i),eqmc,eqsc,eqec);
	fprintf(fid,'%12.6f %12.6f %d  %s %s %s\n',LAT(j,i),LON(j,i),novlpf(j,i),eqmc,eqsc,eqec);
        %time1(j,i)=datenum(eqsc)/yr; time2(j,i)=datenum(eqec)/yr;
end
end
fclose(fid);

%else %not tested
eqsc=datestr(time1*yr,26);eqec=datestr(time2*yr,26);
eqm=datestr(timec*yr,26);
fid = fopen('latlontime.dat','w'); %statistics
[X,Y]=meshgrid(xout,yout);
[LAT,LON]=polarstereo_inv(X,Y,[], [],70,-45);
fprintf(fid,'%12.6f %12.6f %d  %s %s %s\n',LAT(:),LON(:),novlpf(:),eqmc{:},eqsc{:},eqec{:});
fclose(fid);
end %if 0



if algorithmin==2 ||algorithmin==3
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=jumpstd;
imagesc(xt,yt,zt);
colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
caxis([-10 10])
axis(rang0*1e-3)
title('Elevation change uncertainty')
print('-dpdf','-r300','jumpstd')   
saveas(gcf,'jumpstd','fig')

figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=jump;
imagesc(xt,yt,zt,'alphadata',(jumpstd < 10)); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
caxis([-50 50])
axis(rang0*1e-3)
title('Elevation change')
print('-dpdf','-r300','jump')   
saveas(gcf,'jump','fig')
end

if algorithmin==1 ||algorithmin==3
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=ratestd;
imagesc(xt,yt,zt);
colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
caxis([-10 10])
axis(rang0*1e-3)
title('Elevation rate uncertainty')
print('-dpdf','-r300','ratestd')   
saveas(gcf,'ratestd','fig')

figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=rate;
imagesc(xt,yt,zt,'alphadata',(ratestd < 10)); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
caxis([-50 50])
axis(rang0*1e-3)
title('Elevation rate')
print('-dpdf','-r300','rate')   
saveas(gcf,'rate','fig')
end

figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=timec;
imagesc(xt,yt,zt); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
if nanmax(zt(zt~=0)) > nanmin(zt(zt~=0))
caxis([nanmin(zt(zt~=0)) nanmax(zt(zt~=0)) ])
end
axis(rang0*1e-3)
print('-dpdf','-r300','eventtime')   
saveas(gcf,'eventtime','fig')

if algorithmin==1 % 1 linear (ice melting); 2 constant (landslides) ; 3 constant + linear (Okmok volcano)
  multi=2; %For linear trend, plot the time span of all used measurements instead of the std (half the duration between two sequential measurements that have big change).
  text1='Start date of all used measurements';
  text2='End date of all used measurements';
else
  multi=1; %
  text1='Date of the closest measurements before the event';
  text2='Date of the closest measurements after the event';
end
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=multi*timestd;
imagesc(xt,yt,zt); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
if nanmax(zt(zt~=0)) > nanmin(zt(zt~=0))
caxis([nanmin(zt(zt~=0)) nanmax(zt(zt~=0)) ])
end
axis(rang0*1e-3)
if algorithmin==1
 title('The time span of used measurements (years)')
else
 title('Estimated time uncertainty (years)')
end
print('-dpdf','-r300','eventtimestd')   
saveas(gcf,'eventtimestd','fig')

if 1
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=time1;
imagesc(xt,yt,zt); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
if nanmax(zt(zt~=0)) > nanmin(zt(zt~=0))
caxis([nanmin(zt(zt~=0)) nanmax(zt(zt~=0)) ])
end
axis(rang0*1e-3)
%title('Start date of used measurements')
title(text1)
saveas(gcf,'eventtime_start','fig')

figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=time2;
imagesc(xt,yt,zt); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
if nanmax(zt(zt~=0)) > nanmin(zt(zt~=0))
caxis([nanmin(zt(zt~=0)) nanmax(zt(zt~=0)) ])
end
axis(rang0*1e-3)
%title('End date of used measurements')
title(text2)
saveas(gcf,'eventtime_end','fig')
end

%Final count of repeats for each pixel.
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
imagesc(xt,yt,novlpf);
colormap jet;colorbar; shading interp; axis equal; view(0,90)
set(gca,'FontSize', 18);
hold on
plot(xeq*1e-3, yeq*1e-3,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
saveas(gcf,'OverlappingCount_Final','fig')

% plot the chosen piece id for each output pixel
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=iselop; %to check if the order is consistent with idun
imagesc(xt,yt,zt);
colormap jet;colorbar; shading interp; axis equal; view(0,90)
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
axis(rang0*1e-3)
title('Chosen piece isel')
saveas(gcf,'Chosenisel','fig')

%% filter the jump
[co]=Landslidefilter(xout,yout,jump,jumpstd,timec,timestd,oflagc,odir);

