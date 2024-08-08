function [datao]=box2mosaic(rangi,flagband,f,fdir,ix);
%Given a box (rangi), to get the mosaicked files (elevation change file, mosaicked DEM and orthoimages)
% flagband, 1 only output change file; 3, output three files in matrix data.z(:,:,1:3).
%  4: 2 band: elevation change and time of change. 
% f, list of filenames for finding the coregistration offsets.

	constant
    	resr2m=2;
    	%rangi=[minx maxx miny maxy];
	rang0=rangi;
	f_org=f;fdir_org=fdir;

	maxlen=max([abs(rangi(2)-rangi(1))*1e-3, abs(rangi(4)-rangi(3))*1e-3]);% maximum length in km
        if maxlen <= 1 % maximum side of the box is <=4 km; small tiles
               resr2m=2;
        else  % very large slump
               resr2m=10;
	       fprintf(['\n This image ',num2str(ix),' is written at 10 m resolution. rangi =',num2str(rangi),'.\n'])
        end

	tx=rang0(1):resr2m:rang0(2);ty=rang0(4):-resr2m:rang0(3);
	xout=tx;yout=ty;
	nx0=length(xout);ny0=length(yout);
	data0.x=xout;data0.y=yout;data0.z=nan(ny0,nx0); %initial matrix
	data0T2.x=xout;data0T2.y=yout;data0T2.z=nan(ny0,nx0); %initial matrix
	M1=ones(size(data0.z));
	if flagband==1
	datao=data0; %initialize
	elseif flagband==3
	  datao=data0;datao.z=nan(ny0,nx0,3); 
	elseif flagband==2||flagband==4
	  datao=data0;datao.z=nan(ny0,nx0,2); 
	end

	xr0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];yr0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];

	buff=100; %extra buffer zone for reading DEMs
	rangwide=[min(xr0)-buff max(xr0)+buff min(yr0)-buff max(yr0)+buff];

    	%Find the subtiles that may cover this given box.
	%refers to ~/arcticdemapp/landslide/code1/getmosaic.m 
	dx=100e3;x0=-4000e3;y0=-4000e3;	dxs=dx/50; %2km Change detection 2 km tiles;
	
	icount=0;clear xc yc
	for i=1:4 %[1,3]
	    icount=icount+1;
	    x=xr0(i);y=yr0(i);

	    xid=floor((x-x0)/dx)+1; yid=floor((y-y0)/dx)+1;
	    xids=floor((x-(x0+(xid-1)*dx))/(dx/2))+1;
	    yids=floor((y-(y0+(yid-1)*dx))/(dx/2))+1;
	    xidss=floor((x-(x0+(xid-1)*dx+(xids-1)*dx/2))/dxs)+1;
	    yidss=floor((y-(y0+(yid-1)*dx+(yids-1)*dx/2))/dxs)+1;
	    tilefile=sprintf('%02d_%02d_%01d_%01d_%02d_%02d.tif',yid,xid,xids,yids,xidss,yidss) ;

	    xc(icount)=x0+(xid-1)*dx+(xids-1)*dx/2+(xidss-1)*dxs; yc(icount)=y0+(yid-1)*dx+(yids-1)*dx/2+(yidss-1)*dxs; %bottom left
	end

	if 0
    	[tilenamell]=gettilename(xr0(1),yr0(1),1);
    	[tilenameru]=gettilename(xr0(3),yr0(3),1);
    	if ~strcmp(tilenamell,tilenameru)
    	  fprintf(['This polygon encompasses two different subtiles. tilenamell tilenameru',' ',tilenamell,' ',tilenameru,'.\n']);
    	end
	end % if 0

	icount=0;
	for xi=min(xc):dxs:max(xc)
	    for yi=min(yc):dxs:max(yc)
	        icount=icount+1;
		[tilename]=gettilename(xi,yi,1);
		
		% recover the tile range in change detecion
		width=100;%buffer width m
		ifile=tilename;	
		r=1;
		xid= sscanf(ifile(r+3:(r+4)), '%g', 1);
		yid= sscanf(ifile(r:(r+1)), '%g', 1);
		xids= sscanf(ifile(r+6), '%g', 1);
		yids= sscanf(ifile(r+8), '%g', 1);
		xidss= sscanf(ifile(r+10:r+11), '%g', 1);
		yidss= sscanf(ifile(r+13:r+14), '%g', 1);
		x=x0+(xid-1)*dx+(xids-1)*dx/2+(xidss-1)*dxs;y=y0+(yid-1)*dx+(yids-1)*dx/2+(yidss-1)*dxs;
		rang0ti=[x-width x+dxs+width y-width y+dxs+width]; %tile boundary with buffer width
		xr0i=[rang0ti(1) rang0ti(2) rang0ti(2) rang0ti(1) rang0ti(1) ];
		yr0i=[rang0ti(4) rang0ti(4) rang0ti(3) rang0ti(3) rang0ti(4) ];
		
		%coverage of this tile over area of interest
		idx=round((xr0i-data0.x(1))/resr2m)+1;idy=round((yr0i-data0.y(1))/(-resr2m))+1;
		sm=poly2mask(idx,idy,ny0,nx0); % fast, apply to each polygon one by one.
		%sm:  1 has data; 0 has no dem data.
		sm=double(sm);
		ratio=sum(sum(sm))/sum(M1(:))*100; %percentage of tile pixels covering the input box.
		smg{icount}=sm;

	        ratiog(icount)=ratio;	
		tilenameg{icount}=tilename;
	    end
	end % for xi
	
	%sort tiles in descending order of ratio
	[~,idsort]=sort(ratiog,'descend'); % sorted from high to low
	ratiog=ratiog(idsort); tilenameg=tilenameg(idsort);smg=smg(idsort);

	smac=false(size(M1)); %accumulative coverage
	ortho.z=nan(ny0,nx0,4);dem.z=nan(ny0,nx0,4);

	fprintf(['\n This box ',num2str(ix),' encompasses ',num2str(length(ratiog)),' different subtiles.','\n'])
	% 84 subtiles take 0.6 days; 252 subtiles take 5.9 days.
	if length(ratiog) > 100
	   fprintf(['\n This image ',num2str(ix),' is too big, and has',num2str(length(ratiog)),' different subtiles, which may take too long; Skip .\n'])
	   return
	end

	for k1=1:length(ratiog) %all tiles
	tilenamell=tilenameg{k1};
	fprintf(['\n Checking tile ',tilenamell,', which covers ',num2str(ratiog(k1)),'%% of the input box area.'])
    	% 36_22_2_1/36_22_2_1_01_01/36_22_2_1_01_01_jump.tif
       %	num2str(rangi)
    	%file1=[tilenamell(1:9),'/',tilenamell,'/',tilenamell,'_jump.tif'];
        file1=[regiondir,'/',tilenamell(1:9),'/',tilenamell(1:12),'/',tilenamell,'/',tilenamell,'_jump.tif']
	if ~exist(file1,'file')
	  fprintf(['\n ',file1,' does not exist. \n'])
	  if k1>=2
	     fprintf(['\n ',file1,' does not exist, work on next tile. \n'])
	  end
	continue
	end
    	data1=readGeotiff(file1);
%   	idx1=find(data1.x>=rangi(1)&data1.x<=rangi(2));
%   	idy1=find(data1.y>=rangi(3)&data1.y<=rangi(4));
%   	data.x=data1.x(idx1);data.y=data1.y(idy1);
%   	data.z=data1.z(idy1,idx1);
	data=data1;
        %change data to 2 m resolution;-> 2m yields better accuracy for Detecton2
	
	[ny,nx]=size(data.z);
	if min([nx,ny])>=2
	    data.z(data.z <= -9999) = NaN; % %convert nodata values to nan for interpolation
	    tz =interp2(data.x,data.y,double(data.z),data0.x,data0.y','*linear',nan);%toc;%4s
%	    M=~isnan(tz);
%	    data0.z(M)=tz(M);
            tzpre=data0.z;
            M=isnan(tzpre)&~isnan(tz); %Only add the part where previous mask is nan, and this mask provides value.
            fprintf(['\nAdding extra ',num2str([sum(M(:))]),' pixels (',num2str((sum(M(:))/sum(M1(:)))*100),'%%) for elevation change.\n']);
            tzpre(M)=tz(M);
	    data0.z=tzpre;
	   
	    %add time of change
	    T2file=strrep(file1,'_jump.tif','_eventtimeT2.tif');

	    if ~exist(T2file,'file')
              fprintf(['\n file does not exist:',T2file,'; Skip this image ',num2str(ix),'. \n'])
	      return
   	    end

	    dataT2=readGeotiff(T2file);
	    dataT2.z(dataT2.z<=20000101)=NaN;
	    tz=interp2(dataT2.x,dataT2.y,double(dataT2.z),data0.x,data0.y','*nearest',nan);
	    tzpre=data0T2.z;
	    M=isnan(tzpre)&~isnan(tz);
	    tzpre(M)=tz(M);
	    data0T2.z=tzpre;
	else
	    fprintf('Elevation change file is empty.\n')
	end

	% time map	
           T2m=data0T2.z;
           T2c=nan(size(T2m));
           for ii=1:length(T2m(:));
	     if ~isnan(T2m(ii))
		     %MM taken as minutes instead of month.
             %T2c(ii)=(datenum(num2str(T2m(ii)),'YYYYMMDD')-datenum('20000101','YYYYMMDD'))/365.25*10; % scale it to 0 to 255.
             T2c(ii)=(datenum(num2str(T2m(ii)),'yyyymmdd')-datenum('20000101','yyyymmdd'))/365.25*10; % scale it to 0 to 255.
	     end
           end
           T2c=reshape(T2c,size(T2m));
           %end of time map

	if flagband==1 
	   datao=data0;
	   %use the actual datao;
	   smac=~isnan(datao.z);
	elseif flagband==4
           datao.z(:,:,1)=data0.z;
           datao.z(:,:,2)=T2c;
           %use the actual datao;
           smac=~isnan(datao.z);
	elseif flagband==3||flagband==2
	   %datao=data0;
	   datao.z(:,:,1)=data0.z;
	%dir1=[tilenamell(1:9),'/',tilenamell,'/'];
	dir1=[regiondir,'/',tilenamell(1:9),'/',tilenamell(1:12),'/',tilenamell,'/'];
        %61_25_1_1/61_25_1_1_01/61_25_1_1_01_10/61_25_1_1_01_10_listused.txt
	file1=[dir1,'dX4Sg.mat'];	
	if ~exist(file1,'file')
	  listfile=[dir1,tilenamell,'_listused.txt']; %34_23_2_2/34_23_2_2_25_25/34_23_2_2_25_25_listused.txt
          fprintf(['\n file does not exist:',file1,'; Read pxyz in ',listfile,' file instead. \n'])
	  [f, fdir,idregion,dX4Sg,oflag]=getdX4Sg(listfile);
	  if oflag==0 %bad
             %fprintf('Failed retrieving pxyz from list.txt, skip this tile.\n')
            % continue
             fprintf('\n Failed retrieving pxyz from list.txt, download the mosaic tile.\n')
	  end
        else
          oflag=1;
  	  dX4Sgmat=load([dir1,'dX4Sg.mat']); %get idregion dX4Sg  %save dX4Sg.mat dX4Sg idregion
				%save dX4Sg.mat dX4Sg idregion XYbg f fdir -v7.3
	  if isfield(dX4Sgmat, 'fdir') %dX4Sg.mat include f fdir
	    fprintf('Use f fdir in dX4Sg.mat.\n')
	    fdir=dX4Sgmat.fdir;f=dX4Sgmat.f;
	  else
	    fprintf('Use f fdir in mat0.mat.\n')
	    fdir=fdir_org;f=f_org;
	  end
	    dX4Sg=dX4Sgmat.dX4Sg;idregion=dX4Sgmat.idregion;
	end

	%update the filepath in fdir
	fdir=strrep(fdir,'/home/chunli/chunliproject/ArcticDEMdata','/blue/chunlidai/data/ArcticDEMdata');
	fdir=strrep(fdir,'/home/chunli/chunliwork/ArcticDEMdata','/blue/chunlidai/data/ArcticDEMdata');

        if oflag==0
        fprintf('\n Get DEM from the downloaded DEM mosaic.\n')
        [dem1]=getmosaic(rang0);
        dem1.z(dem1.z==-9999)=nan;
        tz=interp2(dem1.x,dem1.y,double(dem1.z),data0.x,data0.y','*linear',nan);
%         dem.x=data0.x;dem.y=data0.y;
        
        k=1;
        tzpre=dem.z(:,:,k);
        M=isnan(tzpre)&~isnan(tz); %Only add the part where previous mask is nan, and this mask provides value.
        fprintf(['\nAdding extra ',num2str([sum(M(:))]),' pixels (',num2str((sum(M(:))/sum(M1(:)))*100),'%%) from DEM mosaic.\n']);
        tzpre(M)=tz(M);
        dem.z(:,:,k)=tzpre;

        else %oflag==1

        fprintf('\n Mosaic the orthoimage and DEM from the latest four strip files.\n')
        %j=1; i=idregion(j); f(i) corresponds to dX4Sg(j)

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
		
	js=max([1,length(idregion)-3]);
	n3=length(idregion)-js+1; %maximum 4
	%ortho.z=nan(ny0,nx0,4);dem.z=nan(ny0,nx0,4);
	k=0;
	for j=js:length(idregion) %[3,7,8,10]; %js:length(idregion) %
		%k=j-js+1;
		k=k+1; %good for picking any j.
	        i=idregion(idblock(j));
	        demdir=fdir{i};
	        infile= [demdir,'/',f{i}];

		p=dX4Sg(idblock(j),:) ; %pzxy
		fprintf(['\n Loading its translational offsets pzxy ',num2str(p),' ',infile])

	metafile=infile;
    if strcmp(metafile(end-7:end),'meta.txt')
    str1='meta.txt';
    elseif strcmp(metafile(end-7:end),'_mdf.txt')
    str1='mdf.txt';
    else
        fprintf(['\n warning: metafile name is not as expected: ',metafile,'.\n'])
    end

	
                demfile= strrep(infile,str1,demext);
		clear demk
                if ~exist(demfile,'file');
	            fprintf(['\n file does not exist:', demfile])
                    continue
                end
	        demk=readGeotiff(demfile,'map_subset',rangwide);
	        [n1,m1]=size(demk.z);
	        if min([n1,m1])<=2
		fprintf('\n No overlapping.\n')
		continue
	        end

     		demk.z(demk.z <= -9999) = NaN; % %convert nodata values to nan for interpolation
        	tz=interp2(demk.x-p(2),demk.y-p(3),double(demk.z)-p(1),data0.x,data0.y','*linear',nan);
	        tzpre=dem.z(:,:,k);
	        M=isnan(tzpre)&~isnan(tz); %Only add the part where previous mask is nan, and this mask provides value.
                fprintf(['\nAdding extra ',num2str([sum(M(:))]),' pixels (',num2str((sum(M(:))/sum(M1(:)))*100),'%%) for DEM strip.\n']);
	        tzpre(M)=tz(M);
	        dem.z(:,:,k)=tzpre;

		orthofile= strrep(infile,str1,'ortho.tif');
                if ~exist(orthofile,'file');
            fprintf(['\n file does not exist:', orthofile])
                        continue
                end
		clear orthok
	        orthok=readGeotiff(orthofile,'map_subset',rangwide);
		orthok.z=double(orthok.z); %int16(nan)=0;
	     	orthok.z(orthok.z==0)=NaN; % %convert nodata values to nan for interpolation
        	tz=interp2(orthok.x-p(2),orthok.y-p(3),double(orthok.z),data0.x,data0.y','*linear',nan);
		% orthoimage has to have coverage >98% to avoid incoherent mosaic.	
		%hi Beware: this may cause problem for very large slump. -solved with thres.
		Mt1=~isnan(tz);
		ratio1=sum(Mt1(:))/length(tz(:))*100; %ratio of valid pixels

	        if maxlen <= 4 % maximum side of the box is <=4 km; small tiles
		   thres=98;
	        else  % very large slump
	 	   thres=10;
	        end

	       	if ratio1 <thres
		    fprintf(['\n skip this orthoimage, valid pixels ratio ',num2str(ratio1),' < ',num2str(thres),'%%. \n'])
		    continue
		end 
		tzpre=ortho.z(:,:,k);
                M=isnan(tzpre)&~isnan(tz); %Only add the part where previous mask is nan, and this mask provides value.
		%fprintf(['\nAdding extra pixels:',num2str([sum(M(:))]),' ',num2str((sum(M(:))/sum(M1(:)))*100),'%% for orthoimage.\n']);
                fprintf(['\nAdding extra ',num2str([sum(M(:))]),' pixels (',num2str((sum(M(:))/sum(M1(:)))*100),'%%) for orthoimage.\n']);
                tzpre(M)=tz(M);
	    	ortho.z(:,:,k)=tzpre;
		%Beware: For different subtiles, the orhoimages may have different dates in different pixels! But this is the best we can do.
	end %j

        end %Done with mosaicking DEMs and Orthoimages.

%	save test1.mat -v7.3
        dem_median.z=median(dem.z,3,'omitnan'); %DEM;
        dem_median.x=datao.x;dem_median.y=datao.y;
        [~,~,~,totalc] = curvature2(dem_median);

	if flagband==3
%	datao.z(:,:,2)=median(dem.z,3,'omitnan'); %DEM;
%	datao.z(:,:,3)=median(ortho.z,3,'omitnan'); %image %bad mosaic; inconsistent colorbar scale
%	%datao.z(:,:,3)=ortho.z(:,:,end); %pick the latest image; risk of clouds
	   datao.z(:,:,2)=totalc*1e6;   
           datao.z(:,:,3)=T2c;
	elseif  flagband==2
%       dem_median.z=median(dem.z,3,'omitnan'); %DEM;
%	dem_median.x=datao.x;dem_median.y=datao.y;
%	[~,~,~,totalc] = curvature2(dem_median); 
	datao.z(:,:,2)=totalc*1e6; %%scale it by 1e6 to convert the data to the range of [0 255].
	end

%	fprintf('hi1 sub 1')
	
	%use the actual datao;
	smac=~isnan(datao.z(:,:,1))&~isnan(datao.z(:,:,2))&~isnan(datao.z(:,:,end));
	end %if flagband

	%save datao.mat datao ortho

	%if the accumulated coverage from 1: k1 tiles covers all the polygon, stop
	%smi=logical(smg{k1});
	%smac=smac|smi;
	if sum(smac(:))==sum(M1(:))  
	  fprintf(['\n All areas of box are covered from 1 to ',num2str(k1),' tiles. No need to use the rest ',num2str(length(ratiog)-k1),' tiles. \n'])
	  if k1>=2
%     fprintf(['\n This box encompasses ',num2str(k1),' different subtiles.','\n'])
	  end
	  return
	end

	end %for k1

	return
end
