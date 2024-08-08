

str=['wget https://crudata.uea.ac.uk/cru/data//hrg/cru_ts_4.07/ge/grid05/cells/N77.5W42.5/tmp/N75.25W44.75.tmp.txt'];
%name convention: 5 degrees, 72.5 for 70 to 75, 177.5 for 175 to 180
% 0.5 degrees: 179.75 for 179.5 to 180; 71.25 for 71 to 71.5;

cd tempe/

years=2008; yeare=2022;
nm=(yeare-years+1)*12; %180, 2008 Jan to 2022 Dec.
lat=60.25:0.5:85.25;
lon=-179.75:0.5:179.75;
nlat=length(lat);nlon=length(lon);
tempe=nan(nm,nlat,nlon);
monthsel=[7,8]; %selected month
%monthsel=[6,7,8,9]; %selected month
Mk=false(nm,1);
datev=nan(nm,1);

count=0;

for j=1:length(lat) %lat=60.25:0.5:85.25
	latj=lat(j);
	lat5=floor(latj/5)*5+5/2;
	if latj<0
		latc=['S',num2str(latj)];
		lat5c=['S',num2str(lat5)];
	else
		latc=['N',num2str(latj)];
		lat5c=['N',num2str(lat5)];
	end

	for i = 1:length(lon) % -179.75:0.5:179.75
		loni=lon(i);
		lon5=floor(loni/5)*5+5/2;
		if loni<0
			lonc=['W',num2str(abs(loni))];
			lon5c=['W',num2str(abs(lon5))];
		else
			lonc=['E',num2str(abs(loni))];
			lon5c=['E',num2str(abs(lon5))];
        end
        
        str5=[lat5c,lon5c];
        str0d5=[latc,lonc]; %N75.25W44.75
        
        if strcmp(str0d5,'N70.75E178.75') %'N75.25W44.75'
            %double check
%             str0d5
        end
        
        str1=['wget https://crudata.uea.ac.uk/cru/data//hrg/cru_ts_4.07/ge/grid05/cells/',str5,'/tmp/',str0d5,'.tmp.txt'];
%       [status, cmdout]=system(str1);
        
        file1=[str0d5,'.tmp.txt'];
        if exist(file1,'file')
            count=count+1;
            latsv(count)=latj; lonsv(count)=loni;
            %read data
            
            fileID = fopen(file1, 'r');

            % Read the header lines
            headerLines = textscan(fileID, '%s', 7, 'Delimiter', '\n');
            headerLines = headerLines{1};

            % Read the data using importdata
            data = importdata(file1, ' ', 7);

            % Close the file
            fclose(fileID);
            
            %
            for iyr=years:yeare
                for imonth=1:12
                    irow=find(data.data(:,1)==iyr&data.data(:,2)==imonth);
                    iym=(iyr-years)*12+imonth;
                    if isempty(irow)
                        fprintf(['\n ', file1,', No data found for year month :',num2str([iyr,imonth]),'. \n '])
                    else
                        tempe(iym,j,i)=data.data(irow,3);
                    end
                end 
            end

        end

    end %for i = 1:length(lon)
end %j=1:length(lat) 

%collect date vector
for iyr=years:yeare
    for imonth=1:12
        iym=(iyr-years)*12+imonth;
        datestr1=[num2str(iyr,'%.4d'),num2str(imonth,'%.2d'),'01'];
        datev(iym)=datenum(datestr1,'yyyymmdd');

    end
end


save test_temp.mat -v7.3

figure;
plot(lonsv,latsv,'ks')
saveas(gcf,'templatlon','fig')

%monthsel 

load('../Mrts.mat')

count=0;
 for iyr=years:yeare
%                 Mk=false(nm,1);
%                 for imonth=1:12
%                     if ismember(imonth,monthsel)
%                          iym=(iyr-years)*12+imonth;
%                          Mk(iym)=1;
%                     end
%                 end
    count=count+1;
    iym=(iyr-years)*12+monthsel;
    %tempe_yr=tempe(iym,:,:);
    tempe_yr=tempe(iym,Mrts);
%   tempe_yr_med(count)=nanmedian(tempe_yr(:));
    tempe_yr_med(count)=nanmean(tempe_yr(:));
 end
 
 ofile=['tempyear',num2str([monthsel],'%.2d%.2d')];
 figure;
 plot(years:yeare,tempe_yr_med,'>-k')
 title(ofile)
 saveas(gcf,ofile,'fig')
save('tempe_yr_med.mat','years','yeare','tempe_yr_med')

% lat lon temperature map for a selected month
iym=(2015-years)*12+monthsel
for j=1:nlat;for i=1:nlon;temp2(j,i)=tempe(91,j,i);end;end
figure;imagesc(lon(:),lat(:),temp2);colorbar
view(0,-90);title('2015/07')

%temperature time series for a selected points
i=130;j=25;
figure;plot(datev(:),tempe(:,j,i));title(['lat lon:',num2str([lat(j), lon(i)])])
datetick('x','mm/yy')

