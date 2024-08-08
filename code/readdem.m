function [datar,nptsubrt]=readdem(XYbi,metafile,rang0,resr)
%read DEM for a given index i.
% nptsubrt: ratio of good DEMs for whole scene based on filtering mask.
% 1\ get ratio of good DEMs; 2\ apply bitmask; 3\ apply geoid; 4\ resampling

%     demdir=fdir{i};
%   infile= strrep([demdir,'/',f{i}],'meta.txt','dem.tif')
    constant

    if strcmp(metafile(end-7:end),'meta.txt')
    str1='meta.txt';
    elseif strcmp(metafile(end-7:end),'_mdf.txt')
    str1='mdf.txt';
    else
        fprintf(['\n warning: metafile name is not as expected: ',metafile,'.\n'])
        nptsubrt=0;datar.x=[];datar.y=[];datar.z=[];
        return;        
    end

    filegeoid='geoidegm96grid15mpolar.tif';
    flaggeoid =0; %1 apply geoid undulation; 0 do not apply geoid;
    %SPOT 5 DEM, ALOS, SRTM, ASTER; to add (AST)
    if contains(metafile,'SPI')||contains(metafile,'GES')||contains(metafile,'ALPSML')||contains(metafile,'SRTM')
       flaggeoid=1;
	fprintf('\n Apply Geoid correction. \n')
    end

    infile= strrep(metafile,str1,'demfill.tif');
    if ~exist(infile,'file')
      %infile= strrep(metafile,str1,'dem.tif');
      infile= strrep(metafile,str1,demext);
      %if dem.tif file does not exist, try dem_10m.tif.
      if strcmp(demext,'dem.tif') && ~exist(infile,'file')
         infile= strrep(metafile,str1,'dem_10m.tif');
	 fprintf(['\n dem.tif does not exist, use dem_10m.tif instead for ',infile,'\n.'])
      end
    end


    try
    if ~isempty(rang0)
        data=readGeotiff(infile,'map_subset',rang0);
    else
        data=readGeotiff(infile); %% to keep long profiles
    end
    catch e
        fprintf('readdem.m There was an error! The message was:\n%s',e.message);
	nptsubrt=0;datar.x=[];datar.y=[];datar.z=[];
        return;
    end % try

    data.z=double(data.z); %data.z(data.z == -9999) = NaN; nan is 0 in int16

    [m,n]=size(data.z);
    if min([m,n])<=2
        fprintf(['\n data has only two rows/columns: ',metafile,'.\n'])
        nptsubrt=0;datar.x=[];datar.y=[];datar.z=[];
        return;        
    end
    
    
    if flaggeoid==1
        geoid=readGeotiff(filegeoid);
    end
    
%get the mask
    infile= strrep(metafile,str1,'bitmask.tif');
    if exist(infile,'file') 
        if ~isempty(rang0)
            mask=readGeotiff(infile,'map_subset',rang0); %0 is good; 1 is edge; > 1 bad data.
        else
            mask=readGeotiff(infile); %0 is good; 1 is edge; > 1 bad data.
        end
 	% 2m res to 10 m resolution
        M1= interp2(mask.x,mask.y,double(mask.z),data.x,data.y','*nearest',1); %logical
        mask.z= (M1); %compatible with older version

        data.z(mask.z>0)=-9999;
    end

% use matchtag to filter bad points
    if flagmatchtag==1 %use 

        %matchtag file 
        infile1= strrep(metafile,str1,'matchtag.tif');
        if exist(infile1,'file')
	   if ~isempty(rang0)
             mt=readGeotiff(infile1,'map_subset', rang0);
	   else
	     mt=readGeotiff(infile1);     
	   end
	   %dilate matchtag mask by 6 10 pixels
	   M=(imdilate(mt.z==0,ones(6))); %mt: 0 bad, 1 good; M: 1 is void
 	   % 2m res to 10 m resolution
           M1= interp2(mt.x,mt.y,double(M),data.x,data.y','*nearest',1); %logical
           M= logical(M1);

	   data.z(M)=-9999;
        else
           fprintf(['\n Step using matchtag as filter; Matchtag file not exist:',infile1,'\n'])
        end

    end

    % ratio of good DEMs for whole scene
    res=data.info.map_info.dx;
%     XYbi=XYbg{i};
    Xb=XYbi(:,1);Yb=XYbi(:,2);
    if 0 %super slow
        [X,Y]=meshgrid(data.x,data.y);
        in = inpolygon(X,Y,Xb,Yb); % 1 inside the polygon
        nptsubrt=sum(sum(data.z~=-9999))/(sum(in(:))) 
        clear X Y in
    else %fast
        idx=round((Xb-data.x(1))/res)+1;
        idy=round(-(Yb-data.y(1))/res)+1;
        in= poly2mask(idx,idy, length(data.y),length(data.x)); % faster than inpolygon
    	nptsubrt=sum(sum(data.z~=-9999))/(sum(in(:)));
    end
	
%     if nptsubrt<0.5 
%         idd2=[idd2;j];
% %         continue
%     end

%   reduce resolution -> resample data in given resolution
    ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resr];
    ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
    tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
    res=data.info.map_info.dx;
    nsr=resr/res;
    if 1 %nsr <1 
     data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
     tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
     tz(isnan(tz))=-9999; %return to -9999
    %the following get the exactly the same as from interp2
    else %nsr >1 and nsr is integer
    idrxs=find(abs(data.x-ranget(1))<1e-3);idrxe=find(abs(data.x-ranget(2))<1e-3);
    idrys=find(abs(data.y-ranget(4))<1e-3);idrye=find(abs(data.y-ranget(3))<1e-3);
    idrx=idrxs:nsr:idrxe;
    idry=idrys:nsr:idrye;
    dd=[idrx(2:end)-idrx(1:end-1),idry(2:end)-idry(1:end-1)];
    if isempty(dd)||any(abs(dd-nsr)>1e-3);warning('Resized grid is not constant spacing.');end
    if length(tx)~=length(idrx)||length(ty)~=length(idry);warning('Wrong downsizing of DEM.');end
    tz=data.z(idry,idrx); %0.09s
    end
    datar= struct();
    datar.x=tx;datar.y=ty;  datar.z=double(tz);

    if flaggeoid==1 %appy geoid undulation to convert geoid height to ellipsoid height.
    geoidz =interp2(geoid.x,geoid.y,double(geoid.z),tx,ty','*linear');% 65.5m ~ 67m
    datar.z(datar.z == -9999) = NaN;
    datar.z=datar.z+geoidz;
    datar.z(isnan(datar.z))=-9999; %return to -9999
    end
        
    return
end
