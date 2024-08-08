function [Soi]=dividemapsubo(ix,S1,buff,data0,flagrescale,flagband,f,fdir,regiondir,odir1)
% Given shapefiles, to get the volume estimates
% The elevation change maps can be either generated from mosaicking or copied from known folders.

%   odirt1=['slumps_ortho/'];
%     odirt1=odir1; %strrep(odir1,'slumpsub_','slumptif_');

    flagfiles=2; %1 default, generate tif files from mosaicing; 2 copy from existing folders and rename it.
    
    resr=10; %

    ixy=ix; i=ix;

    %% read files
    if flagfiles==1

    % see prepareDetectron2.m changemap2polygon.m
    %get the box in polar coordinates
    if isfield(S1, 'BoundingBox') 
    ll=S1(ix).BoundingBox(1,:);ru=S1(ix).BoundingBox(2,:);
    [llx,lly]=polarstereo_fwd(ll(2),ll(1),[],[],70,-45); %lat lon
    [rux,ruy]=polarstereo_fwd(ru(2),ru(1),[],[],70,-45); %lat lon
    [rlx,rly]=polarstereo_fwd(ll(2),ru(1),[],[],70,-45); %right low
    [lux,luy]=polarstereo_fwd(ru(2),ll(1),[],[],70,-45); %left up
    x=[llx,rlx,rux,lux];y=[lly,rly,ruy,luy];
    elseif isfield(S1,'X')
    	[s1x,s1y]=polarstereo_fwd(S1(ix).Y,S1(ix).X,[],[],70,-45);
	x=s1x;y=s1y;
    end

    rang2=[min(x)-buff max(x)+buff min(y)-buff max(y)+buff];
    rang2=round(rang2/resr)*resr;

           try
           [datao]=box2mosaic_output(rang2,flagband,f,fdir,ix,regiondir,odir1); %test hi
           catch e
           fprintf('\n dividemapsubo.m box2mosaic_output: There was an error! The message was:\n%s',e.message);
           % initialize datao
           end

    elseif flagfiles==2
        % copy files to target
        %'/Users/chunlidai/UFwork/ChangeDetection/finalpolygons/slumptif_arctic/slump1287/slump1287_jump.tif';
        %infiles=['/Users/chunlidai/UFwork/ChangeDetection/finalpolygons/',S1(ix).file]; %hi
        infiles=[S1(ix).file]; %hi
        [indir,name,~]=fileparts(infiles);
        str1=['cp ',indir,'/*txt ',odir1];
        [status, cmdout]=system(str1);

        %read datao;

        [~,slumpname,~]=fileparts(odir1(1:end-1));
        filenameo=[odir1,slumpname,'_jump.tif']; %slumptif_arctic/slumpk/'

        filenamet1i=infiles;
        datao=readGeotiff(filenamet1i);
        str1=['cp ',infiles,'  ',filenameo];
        [status, cmdout]=system(str1);

	filenamet1i_org=filenamet1i;
	filenamet1i=strrep(filenamet1i,'_jumpfilter1filled.tif','_jump.tif')  %Find the correct bitmask for givensite_jumpfilter1filled.tif
	filenamet1=strrep(filenamet1i,'_jump.tif','_jumpstd.tif');
    data=readGeotiff(filenamet1);
    datao.z(:,:,2)=data.z(:,:);
    filenamet1o=strrep(filenameo,'_jump.tif','_jumpstd.tif');
        str1=['cp ',filenamet1,' ',filenamet1o];
        [status, cmdout]=system(str1);

	filenamet1=strrep(filenamet1i,'_jump.tif','_eventtimeT1.tif');
    data=readGeotiff(filenamet1);
    datao.z(:,:,3)=data.z(:,:);
    filenamet1o=strrep(filenameo,'_jump.tif','_eventtimeT1.tif');
        str1=['cp ',filenamet1,' ',filenamet1o];
        [status, cmdout]=system(str1);

	filenamet1=strrep(filenamet1i,'_jump.tif','_eventtimeT2.tif');
    data=readGeotiff(filenamet1);
    datao.z(:,:,4)=data.z(:,:);
        filenamet1o=strrep(filenameo,'_jump.tif','_eventtimeT2.tif');
        str1=['cp ',filenamet1,' ',filenamet1o];
        [status, cmdout]=system(str1);

	filenamet1=strrep(filenamet1i,'_jump.tif','_bitmask.tif');
    data=readGeotiff(filenamet1);
    datao.z(:,:,5)=data.z(:,:);
        filenamet1o=strrep(filenameo,'_jump.tif','_bitmask.tif');
        str1=['cp ',filenamet1,' ',filenamet1o];
        [status, cmdout]=system(str1);

	filenamet1=strrep(filenamet1i,'_jump.tif','_nov.tif');
    data=readGeotiff(filenamet1);
    datao.z(:,:,6)=data.z(:,:);
        filenamet1o=strrep(filenameo,'_jump.tif','_nov.tif');
        str1=['cp ',filenamet1,' ',filenamet1o];
        [status, cmdout]=system(str1);

	try
	filenamet1=strrep(filenamet1i,'_jump.tif','_curvature.tif');
    data=readGeotiff(filenamet1);
    datao.z(:,:,7)=data.z(:,:);
        filenamet1o=strrep(filenameo,'_jump.tif','_curvature.tif');
	
        str1=['cp ',filenamet1,' ',filenamet1o];
        [status, cmdout]=system(str1);
	catch e
            fprintf('\n dividemapsubo.m.m There was an error! The message was:\n%s',e.message);
	end


        %update filenamet1i for output
        filename = filenameo;  %'
        % Find the index of the second occurrence of the file separator '/'
        idx = strfind(filename, '/') ;
	if length(idx)>=3
        filenamet1i=filename(idx(end-2)+1:end) ;
	else
		filenamet1i=filename;
	end

    end

           %% Get mask

	   resrx=mean(diff(datao.x));resry=mean(diff(datao.y));
    Mpc=[];
    Mpc.x=datao.x;Mpc.y=datao.y;

    %store the candidate mask to data0.z;
    % poly2mask may not fully recover the original mask, and unintentionally separate clusters.
    % hence artificially increase the number of clusters.
    [s1x,s1y]=polarstereo_fwd(S1(ix).Y,S1(ix).X,[],[],70,-45);
    Xb=round((s1x-Mpc.x(1))/resrx)+1; 
    Yb=round((s1y-Mpc.y(1))/(resry))+1; %descending        
    XYbi=[Xb(:),Yb(:)]; %n by 2
    id1=find(any(isnan(XYbi),2));
    XYbi(id1,:)=[];
    %[nyi,nxi,nzi]=size(Mpcpre.z);
    nyi=length(Mpc.y);nxi=length(Mpc.x);
    Mpc.z=poly2mask(XYbi(:,1),XYbi(:,2), nyi,nxi); %not accurate

%%  write files


    if flagfiles==1

    	M1=~isnan(datao.z(:,:,1));
	
        %save in all related geotiff files: slumptif_arctic/slumpX/
		%        % jump jumpstd T1 T2 bitmask nov dem/curvature
        [~,slumpname,~]=fileparts(odir1(1:end-1));
        filenamet1i=[odir1,slumpname,'_jump.tif']; %slumptif_arctic/slumpk/'
	filenamet1i_org=filenamet1i;
	if sum(M1(:))~=0 %if there are some valid points, save output files.
        fprintf(['\n Write data for slump',filenamet1i,'.\n'])
        projstr='polar stereo north';
        writeGeotiff(filenamet1i,datao.x,datao.y,double(datao.z(:,:,1)),4,nan,projstr)
	
	filenamet1=strrep(filenamet1i,'_jump.tif','_jumpstd.tif');
	writeGeotiff(filenamet1,datao.x,datao.y,double(datao.z(:,:,2)),4,nan,projstr)

	filenamet1=strrep(filenamet1i,'_jump.tif','_eventtimeT1.tif');
	writeGeotiff(filenamet1,datao.x,datao.y,int32(datao.z(:,:,3)),3,0,projstr)

	filenamet1=strrep(filenamet1i,'_jump.tif','_eventtimeT2.tif');
	writeGeotiff(filenamet1,datao.x,datao.y,int32(datao.z(:,:,4)),3,0,projstr)

	filenamet1=strrep(filenamet1i,'_jump.tif','_bitmask.tif');
	writeGeotiff(filenamet1,datao.x,datao.y,int32(datao.z(:,:,5)),3,0,projstr)

	filenamet1=strrep(filenamet1i,'_jump.tif','_nov.tif');
	writeGeotiff(filenamet1,datao.x,datao.y,int32(datao.z(:,:,6)),3,0,projstr)

	filenamet1=strrep(filenamet1i,'_jump.tif','_curvature.tif');
	writeGeotiff(filenamet1,datao.x,datao.y,double(datao.z(:,:,7)),4,nan,projstr)
	else
%       fprintf(['\n No valid data for slump',filenamet1i,'.\n'])
	end %if 
    end %flagfiles

    %% calculate volumes
    % see also /home/chunlidai/blue/apps/landslide/code1/volcanovolume/jump2vol.m
    ptx=nanmedian(s1x);pty=nanmedian(s1y); %point
    [latpt,lonpt]=polarstereo_inv(ptx,pty,[],[],70,-45); %lat lon
    
    jump=datao.z(:,:,1); jumpstd=datao.z(:,:,2);M=Mpc.z;
    %fill the nans for std
    meanstd = nanmedian(jumpstd(M));
    jumpstd(jump==0|isnan(jumpstd))=meanstd;

    vol=-abs(resrx*resry)*sum(jump(M));
    fprintf(['\n Total volume: ',num2str(vol),'m^3, ',num2str(vol*1e-9),' km^3.  \n'])
    
    %box of four zones for calculating the covariance of elevation errors based on the periodogram method.
    
    rangeg=[min(datao.x) min(datao.x)+buff min(datao.y) min(datao.y)+buff;
    max(datao.x)-buff max(datao.x) min(datao.y) min(datao.y)+buff;
    min(datao.x) min(datao.x)+buff max(datao.y)-buff max(datao.y);
    max(datao.x)-buff max(datao.x) max(datao.y)-buff max(datao.y);
        ]; 
    rangeg=rangeg';
    
    [volstd,cor]= corrv2(datao.x,datao.y,jump,jumpstd,M,rangeg); % 2.6764e6 m^3
    fprintf(['\n The uncertainty of total volume considering correlation: ',num2str(volstd) , ' m^3.\n']);
    close all
    
    yr =365.25; 
    T1=datao.z(:,:,3);T2=datao.z(:,:,4);
    minT1=nanmin(T1(M));maxT2=nanmax(T2(M));
    timespan=[num2str(minT1),' to ',num2str(maxT2)];
    dt=(datenum(num2str(maxT2),'yyyymmdd')-datenum(num2str(minT1),'yyyymmdd'))/yr; %in years
    if ~(dt<=100)
        fprintf(['\n Time span is wrong, longer than 100 years, dt=',num2str(dt),'. Check slump',filenamet1i,'.\n'])
    end

    vol
    volstd
    dt

       try
    totalvolperyear=[vol volstd]/dt;
           catch e
           fprintf('\n dividemapsubo.m There was an error! The message was:\n%s',e.message);
	   totalvolperyear=[nan nan];
           % initialize datao
       end

%   save test1.mat -v7.3

    %plot volume per year as a function of year

    % get the delta t in matrix
%     tic %slow 103 seconds
%     dtmat1=( datenum(num2str(T2(:)),'yyyymmdd')-datenum(num2str(T1(:)),'yyyymmdd') )/yr;
%     dtmat=reshape(dtmat1,size(T1));
%     toc
    %get the time tag in terms of year
% If a mass loss were mapped before the summer (<June 20) of a given year, it’s counted as the thawing of the previous year.
    tic
    T2year=floor(T2*1e-4); %20160619 -> 2016 
    %20160619 -> 2015. 20160620 -> 2016
    T1year=floor(T1*1e-4);lastFourDigits=mod(T1,10000);Msm=lastFourDigits<620;T1year(Msm)=T1year(Msm)-1;
    lastFourDigits = mod(T2, 10000);Msm=lastFourDigits<620;T2year(Msm)=T2year(Msm)-1;
    binyear=min(T2year(M)):max(T2year(M));
    binyear=binyear(:);
    nyr=length(binyear);
    volperyr=nan(nyr,1); volperyrstd=nan(nyr,1);
    T1k=nan(nyr,1);T2k=nan(nyr,1);
    idd=[];
    for k=1:nyr
        yeark=binyear(k);
        Myr=T2year==yeark;
        M1=M&Myr;
        if sum(M1(:))==0
            fprintf(['\n No pixels for year T2 =',num2str(yeark),'.\n'])
            continue
        end
	T2k(k)=median(datenum(num2str(T2(M1)),'yyyymmdd')); %
	T1k(k)=median(datenum(num2str(T1(M1)),'yyyymmdd'));
%       dtk=(median(datenum(num2str(T2(M1)),'yyyymmdd'))-median(datenum(num2str(T1(M1)),'yyyymmdd')))/yr; %median time span  %outvolume_sv2
%	dtk=median(T2year(M1))-median(T1year(M1));  % outvolume_sv3; bad for 20210802 to 20220416

	% count the number of summers (June 20 to September 22) within the date span; 
	%like 20200728 to 20200802, 20190604 to 20200308 , 20190301 to 20200502
	cyrT1=str2num(datestr(T1k(k),'yyyy')); cyrT2=str2num(datestr(T2k(k),'yyyy')); %calendar year
	str1=[num2str(cyrT1),'0922'];dayinsum=(datenum(str1,'yyyymmdd')-T1k(k))/(4*30);
	if dayinsum>=1/4; sumyrT1=cyrT1;else;sumyrT1=cyrT1+1;end %days in summer should be > 1/4 of the summer
	str1=[num2str(cyrT2),'0620']; dayinsum=(T2k(k)-datenum(str1,'yyyymmdd'))/(4*30);
	if dayinsum>=1/4; sumyrT2=cyrT2;else sumyrT2=cyrT2-1;end
	dtk=sumyrT2-sumyrT1+1; %outvolume_sv4

        fprintf(['\n Time span for year T2 =',num2str(yeark),' is ',num2str(dtk),' years: ',datestr(T1k(k)),' to ', datestr(T2k(k)),'.\n'])
	if dtk <=0
		idd=[idd;k];
		continue
	end
        volperyr(k)=-abs(resrx*resry)*sum(jump(M1))/dtk;
        volstdk= corrv2(datao.x,datao.y,jump,jumpstd,M1,rangeg,cor);
        volperyrstd(k)= volstdk/dtk;
    end
    toc
    close all

    medianvolperyr=[nanmedian(volperyr),nanmedian(volperyrstd)]; %better than total volume per year estimate, i think.
    fprintf(['\n ',filenamet1i,num2str(ix),' Volume per year total (m^3/yr) =',num2str(totalvolperyear),', median (m^3/yr) = ',num2str(medianvolperyr),'.\n'])

    %Plot time series
    figure;
    set(gcf,'Color','white')
    set(gca,'FontSize', 14);
    set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
    hold all
%     area(binyear,volperyr)
    errorbar(binyear,volperyr,volperyrstd,'r.','Markersize',14,'linewidth',2)
    xlabel('Year')
    box on
    ylabel('Volume (m^3/yr)')
%   axis([min(binyear)-1 max(binyear)+1 0 max(volperyr)*1.2])
%     text(max(epochp)-(max(epochp)-min(epochp))*0.3,max(T6)-(max(T6)-min(T6))*0.4 ,['Trend=',num2str(rate,formatSpec),'\pm',num2str(ratestd,formatSpec),'m/yr'],'FontSize',12)
    close all

%     figure;histogram(T1(M)) ;title('T1');figure;histogram(T2(M)) ;title('T2');
    
	if 0 % elseif flagband ==3
	  A1(:,:,2)=uint8((datao.z(:,:,2))); % %curvature
	  A1(:,:,3)=uint8((datao.z(:,:,3))); % %time
    end

    %collect outputs
    %Soi=S1(ix);
    %delete some fields
    %Soi = rmfield(Soi, 'Feature_Ty'); 
    %Soi = rmfield(Soi, 'layer');    
    %Soi = rmfield(Soi, 'path');
    Soi(1).Geometry=S1(ix).Geometry; %Soi(1).BoundingBox=S1(ix).BoundingBox;
    Soi(1).X=S1(ix).X; Soi(1).Y=S1(ix).Y;
    if isfield(S1, 'BoundingBox')
    Soi(1).BoundingBox=S1(ix).BoundingBox;
    end
    Soi(1).SHAPE_Area=abs(resrx*resry)*sum((M(:)));
    Soi(1).zone=S1(ix).zone;
    

    Soi(1).lat=latpt; Soi(1).lon=lonpt; %point cooridnates

    if flagfiles==1
    matchingFiles =dir(regiondir);
    if ~isempty(matchingFiles) %file found
                regionid=matchingFiles(1).name;
    else
                file1b=strrep(regiondir,'/blue/chunlidai/chunlidai/landslide/','/orange/chunlidai/results/landslide/');
                matchingFilesb =dir(file1b);
                if ~isempty(matchingFilesb)
                        regionid=matchingFilesb(1).name;
		else
			regionid='';
                end
    end
    elseif flagfiles==2
	   if isfield(S1, 'regionid')
	        regionid=S1(ix).regionid;
	   else
	   regionid='';
           end
    end
    Soi(1).regionid=regionid;
    Soi(1).volume=double(vol);
    Soi(1).volume_uncertainty=double(volstd);
    Soi(1).year=double(binyear);
    Soi(1).volume_per_year=volperyr;
    Soi(1).volume_per_year_uncertainty=volperyrstd;
    Soi(1).median_volume_per_year=medianvolperyr(1);
    Soi(1).median_volume_per_year_uncertainty=medianvolperyr(2);
%     Soi(1).timespan=timespan;
    Soi(1).dh=double(jump(M)); %elevation change within the scar area;
    Soi(1).T1bin=double(T1k);
    Soi(1).T2bin=double(T2k);

    %Soi(1).file=filenamet1i;
    Soi(1).file=filenamet1i_org; %hi

    if flagfiles==2 && isfield(S1, 'Comment')
        Soi(1).Comment=S1(ix).Comment;
    end

%   save test1.mat -v7.3

%     exit

return
end
