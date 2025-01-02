
function [yearv,tempe_yr_med]=latlon2precip(latin,lonin,flagmethod,Soi,monthsel)
%given lat lon, output precipitation time series.
flagdata=2; %1, CRU TS (Climatic Research Unit gridded Time Series, v. 4.07) surface air temperature dataset (Harris et al., 2020)
    %       2,ECMWF Reanalysis version 5 (ERA5).

%flagmethod=2; %1, regular year interval; 2, use the actual data interval
shift=365; %0, shift in days for correlation lag.

if flagdata==1 %not updated
    %Gettemperature.m years yeare tempe
load test_temp.mat

% consistent with grids in gettemperature.m

lat=60.25:0.5:85.25;
lon=-179.75:0.5:179.75;
nlat=length(lat);nlon=length(lon);
Mrts=false(nlat,nlon);

nall=length(latin);
for i=1:nall
    lati=latin(i);loni=lonin(i);
    ilat=round((lati-lat(1))/0.5)+1;
    ilon=round((loni-lon(1))/0.5)+1;
    Mrts(ilat,ilon)=1;
end

year_s=2008;
yearv=year_s:yeare;
yearv=yearv';
tempe_yr_med=nan(size(yearv));
count=0;
 for k=1:length(yearv) %iyr=years:yeare
    count=count+1;
    iyr=yearv(k);
    iym=(iyr-year_s)*12+monthsel;
    %tempe_yr=tempe(iym,:,:);
    tempe_yr=tempe(iym,Mrts);
%     tempe_yr_med(k)=nanmedian(tempe_yr(:));
    tempe_yr_med(k)=nanmean(tempe_yr(:));
 end
 
elseif flagdata==2
    %
    %2m temperature;
    ncFile='./era5/ERA5precip.nc'; %2008 to 2023, June July August September, hourly
    %Precipitation:2008 to 2023, April May June July August September, hourly
        flagmon=0; % 1, input is montlhy data; 0 hourly data
    
    if flagmethod==3
%         ncFile='./era5/ERA5preciGreenland2010to22hr.nc';
        if flagmon==1
            ncFile='/Users/chunlidai/test/era5/monthly_temppreci.nc'; % lon 0 to 360
            if lonin<0; lonin=lonin+360;end
        elseif flagmon==0
            ncFile='/Users/chunlidai/test/era5/ERA5preciGreenland2009to13hr.nc';
%             ncFile='/Users/chunlidai/test/era5/ERA5preciGreenland2014to19hr.nc';
        end
    end
    fileInfo = ncinfo(ncFile);

    %0.25 degrees
    lat = ncread(ncFile, 'latitude');
    lon = ncread(ncFile, 'longitude');
    time = ncread(ncFile, 'time'); % hours since  elapsed since January 1, 1900.
    timed=datenum('19000101','yyyymmdd')+double(time)/24; %days
    %hour to calendar year 
    timeyr=str2num(datestr(timed,'yyyy')); %1900+double(time)/(365.25*24); 
    timemon=str2num(datestr(timed,'mm')); 

    nlat=length(lat);nlon=length(lon);
%     Mrts2=false(nlon,nlat);

    dlat=median(diff(lat));dlon=median(diff(lon));
    nall=length(latin);
    
    %output 
    if flagmethod==1
    year_s=min(timeyr);yeare=max(timeyr);
    yearv=year_s:yeare;
    yearv=yearv';
    elseif flagmethod==2
        yearv=Soi.year;% T1bin T2bin (datenum days).
            
    elseif flagmethod==3
          yearv=(timed)/365.25; %years from datenum %Original interval; hourly or monthly
%           yearv=(floor(min(timed)):1:floor(max(timed)))/365.25; % daily interval, in unit of years
    end
    tempe_yr_med_n=nan(length(yearv),nall);

%     monthsel=[7,8]; %selected month

    for i=1:nall
        lati=latin(i);loni=lonin(i);
        ilat=round((lati-lat(1))/dlat)+1;
        ilon=round((loni-lon(1))/dlon)+1;
%         Mrts2(ilon,ilat)=1;
    
        %Temperature measured in kelvin can be converted to degrees Celsius (°C) by subtracting 273.15.
        %1440x121x46848, lon lat time
        if 1
        start  = [ilon ilat 1]; 
        count  = [1 1 Inf]; 
        stride = [1 1 1]; %Space between variable indices
        if flagmon==1 % flagmethod==3
            % 1440 721 2 199        
            start  = [ilon ilat 1 1]; 
            count  = [1 1 1 Inf]; 
            stride = [1 1 1 1]; %Space between variable indices

        end
        tempe = ncread(ncFile, 'tp',start,count,stride);
        % Sep. 2024; scale and offset; already applied in ncread.m

%          if flagmethod ==3
%              tempe_yr_med_n(:,i)=tempe(1,1,:);
%              continue
%         end

        %get the mean of July August
        
        for k=1:length(yearv) %iyr=years:yeare
            iyr=yearv(k);
            Mtm=false(size(timeyr));
            for kj=1:length(monthsel)
                Mtm=Mtm|(timemon==monthsel(kj));
            end
            
            if flagmethod==1
                Mty=timeyr==iyr;
    
                Mt=Mty&Mtm;
            elseif flagmethod ==2 
                T1k=Soi.T1bin(k);T2k=Soi.T2bin(k);
                Mt12=timed>=T1k-shift&timed<=T2k-shift;
                Mt=Mt12&Mtm;
            elseif flagmethod==3
%                 Mt=floor(timed)==yearv(k)*365.25; % find all data within the selected day.
                Mt=(timed)==yearv(k)*365.25; % original interval
            end
            if flagmon==1 %4D %flagmethod==3 %4D
                tempe_yr=tempe(1,1,1,Mt);
            else %3D
                tempe_yr=tempe(1,1,Mt);
            end
            %should replace NaN with zeros.%Sep. 2024
            M1=isnan(tempe_yr);tempe_yr(M1)=0;
            tempe_yr_med_n(k,i)=nanmean(tempe_yr(:));
        end
    end
    
    %Unit: The units of this parameter are depth in metres of water equivalent within an hour per 0.25 degrees grid. 
% https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

    tempe_yr_med=nanmean(tempe_yr_med_n,2);
    
    %get 2D map
    if 0
        start  = [1 1 1]; 
        count  = [nlon nlat 1]; 
        stride = [1 1 1]; %Space between variable indices
        tempe = ncread(ncFile, 'tp',start,count,stride);
        figure;imagesc(lon,lat,tempe');hold on;plot(lonin,latin,'k>')
        colorbar;view(0,-90)
    end
    
    
end %flagdata

return


end
