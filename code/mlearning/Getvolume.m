constant

%New versoin May 2023: Use the given candidates from Change2poly6classespar.m to get test images.
%new version Jan, 2022: use constraint 1 to get candidate clusters, and then get images for all clusters.

%codedir='/home/dai.56/arcticdemapp/landslide/code1/';  %Directory of codes.
%addpath(genpath('/home/chunlidai/blue/apps/landslide/code1/'));  %Directory of codes.

flagorthotif=0; % save mosacked ortho image in Geotiff
flagdemtif=0; % save mosacked dem in Geotiff
flagrescale=0;  %1, fixed range [e.g. -20 10] to keep the information of positive/negative.
                %0, re-adjust the range for each image.
flagband=3; %3; %1 writing 1 band data for training images; 3 writing 3 band data; 2: elevation change and DEM curvature, two bands.

classtypeg={'_m','_l','_vl'}; % _m _l _vl

ncount=0;
% nall=2850; %
% 
% volumeg=cell(nall,1); 
% volumestdg=cell(nall,1);
% timespang=cell(nall,1);
% regionidg=cell(nall,1);
% tiffileg=cell(nall,1);

if 0

for iregid=3:34; %10 %3:34;

%[demdir,name,ext] =fileparts([strtrim(regiondir)]);
%regid=name(11:12); %03 to 34
regid=sprintf('%02d', iregid);
regiondir=['/blue/chunlidai/chunlidai/landslide/','arcticdem_',regid,'*']; % arcticdem_09_canada_victoria; /blue/chunlidai/chunlidai/landslide/arcticdem_09_canada_victoria

regionstr=['region',regid]; %'region09';

for k=1:3 %3  %1:3
	classtype=classtypeg{k};
classtypein=classtype;
shapefile=['finalrts/change2polyg_',regid,classtypein,'_ai_com.shp']; % finalrts/change2polyg_18_vl_ai_com.shp
%mergedfile=[regiondir,'/region',regid,'_jump.tif'];

if ~exist(shapefile)
	continue
end

S1=shaperead(shapefile);
ns=length(S1);

%imagefiles_val={mergedfile};
%changefile=mergedfile; %imagefiles_val{1};

%step: 500 m ; buffer 200 m on each size. for Peel Plateau
% try 800 pixels (2m ): 1600 m ; For Eureka
%buff=100; dx=300; %-> Chandi uses 250 pixels by 250 pixels;
%dy=-dx;

%The largest slump is roughly 600 by 600 m.
buff=600; % best
%buff=2000; %test the effect to the covariance matrix;

fprintf(['There are a total of ',num2str(ns),' slumps ',' in ',shapefile,' .\n'])

ncounts=ncount;
clear Soi
%Soi(ns)=struct();
%areathres=2*2e3*2e3;
sz = getenv('SLURM_NTASKS');
sz=str2num(sz);
fprintf(['\n ',num2str(sz),' worker(s) allocated in job.slurm.\n'])
poolobj=parpool(sz);
parfor ix=1:ns
%for ix=1:ns %1 %1:ns

%	ncount=ncount+1; %does not work with parfor
	ncount=ncounts+ix;

	%save in all related geotiff files: slumptif_arctic/slumpX/
	odir1=['slumptif_arctic/slump',num2str(ncount),'/'];
        if ~exist(odir1,'dir')
  		 mkdir(odir1)
	end

        fprintf(['\n Working on slump ',num2str(ix),' ; ',num2str(ix),'/',num2str(ns),' in ',shapefile,' .\n'])
	tic
	%run 5 classes seperately.
        %[Mpc,idsi]=dividemapsubo(ix,S1,buff,data0,flagrescale,flagband,f,fdir,regiondir,odir1);
        [Soi(ix)]=dividemapsubo(ix,S1,buff,[],flagrescale,flagband,[],[],regiondir,odir1);
	%So(ncount)=Soi;
	
	fprintf(['\n ',num2str(ix),' takes ',num2str(toc),' seconds.\n']);
end %ix
delete(poolobj)

%collect results
for ix=1:ns %1 %1:ns
        ncount=ncounts+ix; %does not work with parfor
	So(ncount)=Soi(ix);
end

end %for k

end %for  iregid 

end 
%save test2.mat -v7.3
load test2_sv4.mat

nall=length(So);

%% plot volume per year time series for selected RTS
%2012 to 2022; final year;
yearf=2012:2022;

%volumegv same as vol
volumeg=cell(nall,1);
for i=1:nall
    volumeg{i}=So(i).volume; 
end
volumegv=cell2mat(volumeg);
[~,idk]=sort(volumegv,'descend');
for i=1:50; regtop10{i,1}=So(idk(i)).regionid;end

Mgood=false(nall,1); %select RTS with >= 9 obs.
Mbad=false(nall,1);
Rcor=nan(nall,3); %correlation with temperature, precipitation 07 08, precipitation 05 06
pvalg=nan(nall,3);
% for k=1:50
%    i=idk(k);
for i=1:nall % 6 hours to read temperature/precipitation data
   %i =  2110; maximum volume and area; region 10
    % i=2100; % covariance figure; area, 114,400 m^2,  very large.
    % i=1893; FM3; area 35,500, median
    %i =2831; % second largest area; region 30
    %i=2221; %second laregest volume in region 12, north of mountain range. too short
    %i=2210; %third largest volume in region 12.

    year_s(i)=min(So(i).year);
    year_e(i)=max(So(i).year);
    vol(i)=So(i).volume;
    volstd(i)=So(i).volume_uncertainty;
    
    binyear=So(i).year;
    vpy=So(i).volume_per_year;
    vpystd=So(i).volume_per_year_uncertainty;
    
    %interpolate the volume per year
    M1=~isnan(vpy);M2=~isnan(vpystd);
    if sum(M1(:))>=9;Mgood(i)=1;  
    elseif sum(M1(:))<=2 %Mbad is not counted toward total volume per year.
        fprintf(['\n slump ',num2str(i),' has only ',num2str(sum(M1(:))),' valid points. \n'])
        Mbad(i)=1;
        continue
    end %
    volperyear_int(:,i)=interp1(binyear(M1),vpy(M1),yearf,'nearest','extrap');
%     volperyear_int(:,i)=interp1(binyear(M1),vpy(M1),yearf,'linear','extrap');   %risky extrapolation
    volperyearstd_int(:,i)=interp1(binyear(M2),vpystd(M2),yearf,'nearest','extrap');

    ctype{i}=class(So(i).year);
    if ~strcmp(ctype{i},'double')
        fprintf(['Class type ',ctype{i},', i=', num2str(i)])
    end
    
    % add temperature time series
    if 1
    flagmethod=2; %1, regular year interval; 2, use the actual data interval
    [yearv,tempe_yr]=latlon2tempe(So(i).lat,So(i).lon,flagmethod,So(i));
    shift=0;  M=yearv>=(yearf(1)-shift)&yearv<=(yearf(end)-shift);
    if flagmethod ==1
        fprintf(['\n Temperature uses the summer mean of each year.\n '])
        [R,pval] = corr(volperyear_int(:,i),tempe_yr(M));
        tempe_yr1=interp1(yearv,tempe_yr,binyear(M1),'nearest');
        [R,pval]=corr(vpy(M1),tempe_yr1); %correlation of the data, not the interpolated data
        
    elseif flagmethod==2
        fprintf(['\n Temperature uses the summer mean of actual data time span.\n '])
        if length(tempe_yr)~=length(vpy); fprintf([' \n  Temperature vector from method 2 seems to have the wrong size. \n ']); end
        [R,pval]=corr(vpy(M1),tempe_yr(M1));
    end
    Rcor(i,1)=R;
    pvalg(i,1)=pval;

    %precipitation;
    monthsel=[7,8];[yearv,precip_yr]=latlon2precip(So(i).lat,So(i).lon,flagmethod,So(i),monthsel);
    if length(precip_yr)~=length(vpy); fprintf([' \n  precipitation vector from method 2 seems to have the wrong size. \n ']); end
    [R2,pval]=corr(vpy(M1),precip_yr(M1));
    Rcor(i,2)=R2; pvalg(i,2)=pval;
    %spring
    monthsel=[4,5,6];[yearv,precip_yrb]=latlon2precip(So(i).lat,So(i).lon,flagmethod,So(i),monthsel);
    [R3,pval]=nancorr(vpy(M1),precip_yrb(M1));
    Rcor(i,3)=R3;pvalg(i,3)=pval;
    
    end

    if 1

    right_color = [0 0 0]; %black
    left_color = [255,0,0]/255; %red
    fig=figure; %figure(5);
    clf(fig) %clear preview plots
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    set(gcf,'Color','white')
%     set(gca,'FontSize', 14);
%     fontsize(gcf,14,"points")
    set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4
    subplot(2,1,1)
    yyaxis left
    hold on
    %     area(binyear,volperyr)
%     errorbar(yearf, volperyear_int(:,i),volperyearstd_int(:,i),'k.-','Markersize',14,'linewidth',1)
    errorbar(binyear,vpy,vpystd,'r.','Markersize',14,'linewidth',1)
    xlabel('Year')
    box on
    ylabel('Each Volume (m^3/yr)')
%     legend('Interpolated','Observed')
    legend('Observed')
    % add temperature time series
    yyaxis right
    plot(yearv,tempe_yr,'>-k')
    ylabel('Temperature (°C)')
%     plot(yearv,precip_yrb*1e3*24*30,'>-b')
%     ylabel('Percipitation (mm monthly)') %    %Unit: The units of this parameter are depth in metres of water equivalent within an hour. 
    title([So(i).regionid, num2str([i,So(i).volume, k, R, R2, R3 ])])
    subplot(2,1,2)
    yyaxis left
    hold on
    errorbar(binyear,vpy,vpystd,'r.','Markersize',14,'linewidth',1)
    xlabel('Year')
    box on
    ylabel('Each Volume (m^3/yr)')
    legend('Observed')
    % add temperature time series
    yyaxis right
    plot(yearv,precip_yr*1e3*24*30,'>-b')
    ylabel('Percipitation (mm monthly)') % 

%     figure;mapshow(So(i))
%     title(num2str([i,So(i).volume, k]))
     pause
    end
end
figure;plot(1:nall,Rcor(:,1),'r.-',1:nall,Rcor(:,2),'b.-',1:nall,Rcor(:,3),'c.-')
figure;histogram(Rcor(:,1));title('Correlation volume per year vs July August temperature')
figure;histogram(Rcor(:,2));title('Correlation volume per year vs July August precipitation')
figure;histogram(Rcor(:,3));title('Correlation volume per year vs April to June precipitation')
[sum(Rcor(:,1)>=0.5)/nall] %1333

save test2.mat -v7.3

fprintf(['\n ',num2str(sum(Mbad)),' rts has <=2 valid obs. \n']); %38; 54 bad ones, say median volume loss per year is 6e5 -> 32e6

%Mbad (38) is not counted to the total volume per year.
totalvolperyear=nan(size(yearf));totalvolperyearstd=nan(size(yearf));
for j=1:length(yearf)
    totalvolperyear(j)=sum(volperyear_int(j,1:nall)); %Mgood
    totalvolperyearstd(j)=sqrt(sum(volperyearstd_int(j,1:nall).^2));
end
num2str([median(totalvolperyear) median(totalvolperyearstd) ]*1e-6) % 32084692.5568      46823.6412962'

totalarea=sum(areav);
% total volume for all slumps
totalvol=sum(vol(:));
totalvolstd=sqrt(sum(volstd.^2));
fprintf(['\n Total volume is , ',num2str([totalvol totalvolstd]),' m^3, for ',num2str(nall),' slumps. \n'])
%Mass unit Gt, gigatonne, 1 Gt = 1 Pg = 10^12 kg
density=917; %ice 917 kg/m3; %kg/m^3
fprintf(['\n Total mass is , ',num2str([totalvol totalvolstd]*density*1e-12),' Gt, for ',num2str(nall),' slumps. \n'])
fprintf(['\n Total mass is , ',num2str([median(totalvolperyear) median(totalvolperyearstd)]*density*1e-12),' Gt/year, for ',num2str(nall),' slumps. \n'])
%to Pg carbon (petagrams of carbon); 1Pg =1 e15g = 1 billion (1e9) tons (1e6 grams) = 1Gt;
% assuming carbon content of 11–14 kg C /m3, which accounts for ground ice
cdensity=12;
fprintf(['\n Total carbon release is , ',num2str([totalvol totalvolstd]*cdensity*1e-12),' Pg carbon, for ',num2str(nall),' slumps. \n'])
fprintf(['\n Total carbon release is , ',num2str([median(totalvolperyear) median(totalvolperyearstd) ]*cdensity*1e-12),' Pg carbon per year. \n'])
% 0.00038502  5.6188e-07 Pg. 
fprintf(['\n Total carbon release is , ',num2str([median(totalvolperyear) median(totalvolperyearstd) ]*cdensity*1e-12*1e15/365/totalarea),' g C m−2 day−1. \n'])

figure;
plot(1:nall,year_s)
hold on;plot(1:nall,year_e)
set(gcf,'Color','white')
set(gca,'FontSize', 14);
set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4
legend('Earlest detected active year','Latest detected active year')

% figure;histogram(year_s)
% figure;histogram(year_e)

%get 0.5 degree grid of the final RTS.
% consistent with grids in gettemperature.m
if 0
lat=60.25:0.5:85.25;
lon=-179.75:0.5:179.75;
nlat=length(lat);nlon=length(lon);
Mrts=false(nlat,nlon);
for i=1:nall
    lati=So(i).lat;loni=So(i).lon;
    ilat=round((lati-lat(1))/0.5)+1;
    ilon=round((loni-lon(1))/0.5)+1;
    Mrts(ilat,ilon)=1;
end
figure;imagesc(lon(:),lat(:),Mrts);colorbar
view(0,-90)
save('Mrts.mat','lon','lat','Mrts')
end

latgv=nan(nall,1);longv=nan(nall,1);
for i=1:nall
    latgv(i)=So(i).lat;
    longv(i)=So(i).lon; 
end
%precipitation; Computation time long, ~14 hours
flagmethod=1;
tic;[yearv,tempe_yr]=latlon2tempe(latgv,longv,flagmethod);toc %4 hours
tempe_yr1=interp1(yearv,tempe_yr,yearf,'nearest');
[R1,pval]=corr(totalvolperyear(:),tempe_yr1(:));
tic; monthsel=[7,8];[yearvp,precip_yr]=latlon2precip(latgv,longv,flagmethod,[],monthsel);toc
precip_yr1=interp1(yearvp,precip_yr,yearf,'nearest');
[R2,pval]=corr(totalvolperyear(:),precip_yr1(:));
%spring not correlated
tic;monthsel=[4,5,6];[yearvpb,precip_yrb]=latlon2precip(latgv,longv,flagmethod,[],monthsel);toc
precip_yrb1=interp1(yearvpb,precip_yrb,yearf,'nearest');
[R3,pval]=corr(totalvolperyear(:),precip_yrb1(:));

%CRU TS (Climatic Research Unit gridded Time Series, v. 4.07) 
t1=load('tempe_yr_mean.mat');
hold all;plot(t1.years:t1.yeare,t1.tempe_yr_med,'.-')
tempe_yr_cru=interp1(t1.years:t1.yeare,t1.tempe_yr_med,yearf,'nearest');
[R1_cru, pval]=corr(totalvolperyear(:),tempe_yr_cru(:));
[R1_cru_era5, pval]=corr(tempe_yr1(:),tempe_yr_cru(:));
% If pval(a,b) is small (less than 0.05), then the correlation rho(a,b) is significantly different from zero / correlated.

save test2.mat -v7.3

% plot volume per year time series
right_color = [0 0 0]; %black
left_color = [255,0,0]/255; %red
fig=figure;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
%     area(binyear,volperyr)
errorbar(yearf,totalvolperyear*1e-6,totalvolperyearstd*1e-6,'r.-','Markersize',14,'linewidth',1)
xlabel('Year')
box on
ylabel('Volume loss rate (10^6 m^3/yr)')
% axis([min(binyear)-1 max(binyear)+1 0 max(volperyr)*1.2])
yyaxis right
% load tempe_yr_med.mat
% plot(years:yeare,tempe_yr_med,'>-k')
plot(yearv,tempe_yr,'>-k')
ylabel('Temperature (°C)')
plot(yearvp,precip_yr*1e3*24*30,'>-b')
ylabel('Percipitation (mm monthly)') %    %Unit: The units of this parameter are depth in metres of water equivalent within an hour. 

set(gcf,'Color','white')
set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4
fontsize(gcf,14,"points")
legend('Volume loss','ERA5','CRU')

saveas(gcf,'totalvolperyear','fig')

%% Scatter plots of Temperature and Volume

figure;
errorbar(tempe_yr1,totalvolperyear*1e-6,totalvolperyearstd*1e-6,'r.','Markersize',14,'linewidth',1)
xlabel('Temperature (°C)')
ylabel('Volume loss rate (10^6 m^3/yr)')
% y =kx+b
lenf=length(totalvolperyear);mp=2;
AM=[tempe_yr1(:),ones(lenf,1)];
yobs=[totalvolperyear(:)*1e-6];
yobserr=totalvolperyearstd(:)*1e-6;
P=diag(yobserr(:).^-2); 
var=inv(AM'*P*AM);
est=var*AM'*P*yobs;

etilde=yobs-AM*est;
sigma02hat=etilde'*P*etilde/(lenf-mp);
var=var*sigma02hat; %estimated variance
fit=AM*est;
fitstdall=AM*var*AM';
fitstd=zeros(length(fitstdall(:,1)),1);
for j=1:length(fitstdall(:,1))
    fitstd(j)=sqrt(fitstdall(j,j));
end
hold all;
[tempe_yr1s,ids]=sort(tempe_yr1);
plot(tempe_yr1,fit,'k-','linewidth',2,'markersize',14)
% shadedErrorBar(tempe_yr1s,fit(ids),fitstd(ids),'k-',1)
set(gcf,'Color','white')
set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4
fontsize(gcf,14,"points")

%%
%calcualte correlation
% yeart=years:yeare;
% for shift=0:4
%     2012-shift
%     M=yeart>=2012-shift&yeart<=2022-shift;
%     R = corrcoef(totalvolperyear(:),tempe_yr_med(M)')
% end
% median temperature of Arctic, correlation 0.175, shift to 2010:2020
% mean, 0.3, at 2009:+10
% mean all rts, 0.4, at 2009:+10;
% median all rts, 0.25, at 2011:+10;

figure;
set(gcf,'Color','white')
set(gca,'FontSize', 14);
set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
errorbar(1:nall,vol,volstd,'r.','Markersize',14,'linewidth',1)
xlabel('Year')
box on
ylabel('Total Volume (m^3) for each RTS')

%% Histogram of area, volume, elevation changes.
areav=nan(nall,1); minh=nan(nall,1);dhpixel=[];
%vol ,
for i=1:nall
    areav(i)=So(i).SHAPE_Area;
    minh(i)=min(So(i).dh);
    dhpixel=[dhpixel(:);So(i).dh(:)];
end
[minminh,id]=min(minh);
num2str([median(areav),median(vol),minminh,id,median(minh),median(dhpixel)])

figure;histogram(areav(:),'BinWidth',1e3);
% figure;histogram(vol(:),'BinWidth',5e3);
% figure;histogram(dhpixel,'Normalization','pdf')
figure;histogram(minh); title('minh')
set(gcf,'Color','white')
set(gca,'FontSize', 14);
set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3

%%
load('test2_sv4Rcor_all.mat')
%vol areav volstd
[alpha,beta]=volumearea(vol,areav,volstd);
% [alpha,beta]=volumearea(vol,areav);
[alpha,beta]=volumearea(vol(1:2),areav(1:2),volstd(1:2));

%% write output

%total volume, median volume per year
output=[];
latg=cell(nall,1);long=cell(nall,1);
volumeg=cell(nall,1);medianvolperyr=cell(nall,1);
Rcorg1=cell(nall,1);Rcorg2=cell(nall,1);Rcorg3=cell(nall,1);
for i=1:nall
%     output=[output;So(i).lon, So(i).lat, So(i).volume, So(i).median_volume_per_year];
    latg{i}=So(i).lat;
    long{i}=So(i).lon; 
    volumeg{i}=So(i).volume; 
    medianvolperyr{i}=So(i).median_volume_per_year;
    Rcorg1{i}=Rcor(i,1);Rcorg2{i}=Rcor(i,2);Rcorg3{i}=Rcor(i,3);
end
% save('vollatlon.gmt','output','-ascii')
Sopt=struct('Geometry', 'Point','X',long,'Y',latg,'volume',volumeg, ...
    'median_volume_per_year',medianvolperyr,...
    'correlation_temperature',Rcorg1,...
    'correlation_precipitation_JulyAugust',Rcorg2,'correlation_precipitation_ApriltoJune',Rcorg3);
shapewrite(Sopt, 'vollatlon.shp');

shpname='arcticrts.shp';
% shapewrite(So, shpname);
% Error using makedbfspec>validateAttributes
% Attribute field year of S contains at least one value that is not a scalar.
% So1=So;
% So1 = rmfield(So1, 'year');
% So1 = rmfield(So1, 'volume_per_year');
% So1 = rmfield(So1, 'volume_per_year_uncertainty');
% shapewrite(So1, shpname);

% recontruct the struct 
BoundingBox=cell(nall,1);xo=cell(nall,1);yo=cell(nall,1);SHAPE_Area=cell(nall,1);
latg=cell(nall,1);long=cell(nall,1);regionid=cell(nall,1);
volumeg=cell(nall,1);volumestdg=cell(nall,1);
binyear=cell(nall,1);volperyr=cell(nall,1);volperyrstd=cell(nall,1);
medianvolperyr=cell(nall,1);medianvolperyrstd=cell(nall,1);fileg=cell(nall,1);scarg=cell(nall,1);
for i=1:nall
    BoundingBox{i}=So(i).BoundingBox;
    xo{i}=So(i).X; yo{i}=So(i).Y;
    SHAPE_Area{i}=So(i).SHAPE_Area; latg{i}=So(i).lat;
    long{i}=So(i).lon; regionid{i}=So(i).regionid;
    volumeg{i}=So(i).volume; volumestdg{i}=So(i).volume_uncertainty;
    binyear{i}=So(i).year; volperyr{i}=So(i).volume_per_year;
    volperyrstd{i}=So(i).volume_per_year_uncertainty; 
    medianvolperyr{i}=So(i).median_volume_per_year;
    medianvolperyrstd{i}=So(i).median_volume_per_year_uncertainty;
    fileg{i}=So(i).file;
    scarg{i}='scar';
end
% array with different size from xo, yo is not okay.
%     'year',binyear,'volume_per_year',volperyr,'volume_per_year_uncertainty',volperyrstd,...
So2=struct('Geometry', 'PolyGon','BoundingBox',BoundingBox, 'X', xo, 'Y', yo ...
    ,'SHAPE_Area',SHAPE_Area,'zone',scarg,'lat',latg,'lon',long,'regionid',regionid,'volume',volumeg, ...
    'volume_uncertainty',volumestdg,...
    'median_volume_per_year',medianvolperyr,'median_volume_per_year_uncertainty',medianvolperyrstd,'file',fileg);
shapewrite(So2, shpname);


