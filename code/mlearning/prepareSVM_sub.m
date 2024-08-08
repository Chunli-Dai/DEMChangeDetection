function [featlabels,featgt,category_id_feat]=prepareSVM_sub(S1,i,buff,shapefiles_k,data0,flagdata,flagband,f,fdir,scenario)
	
	featlabels=[];featgt=[];category_id_feat=[];

    if isfield(S1, 'zone')
    %if contains(S1(i).zone,'debris');continue;end %only detect scars
    %if ~(contains(S1(i).zone,'scar')|| contains(S1(i).zone,'nois'));continue;end   % two classes
    %if ~(contains(S1(i).zone,'scar'));continue;end 
    %if ~(contains(S1(i).zone,'nois') || contains(S1(i).zone,'gull')  );return;end  %four classes: stable, slumps, noise, gullies
    if ~(contains(S1(i).zone,'scar')||contains(S1(i).zone,'nois') || contains(S1(i).zone,'gull')  );return;end  %four classes: stable, slumps, noise, gullies; %run_allseasons_pgcnew/test_svm2_t1_sv2.mat, pgcnew/test_svm2_t1_banks_v1.mat
    end
    ll=S1(i).BoundingBox(1,:);ru=S1(i).BoundingBox(2,:);
    [llx,lly]=polarstereo_fwd(ll(2),ll(1),[],[],70,-45); %lat lon
    [rux,ruy]=polarstereo_fwd(ru(2),ru(1),[],[],70,-45); %lat lon
    x=[llx,rux];y=[lly,ruy];
    rangi=[min(x)-buff max(x)+buff min(y)-buff max(y)+buff];

    %if  rux >=-2485300 %for peel, right 70% (this is the training area) %hi
    if rux <-2495300; % for peel, left 10% (this is the test area).
    else; return;
    end
    
    %Get multiple the polygons within the buffer zone
    if 1 % 
         x0=[rangi(1) rangi(2) rangi(2) rangi(1) rangi(1) ];
         y0=[rangi(4) rangi(4) rangi(3) rangi(3) rangi(4) ];
         [lat,lon]=polarstereo_inv(x0,y0,[],[],70,-45);
         bbox=[min(lon) min(lat); max(lon) max(lat)]; %left bottom ; right top
         tic;S1sub=shaperead(shapefiles_k,'BoundingBox',bbox);toc %10 minutes
    else %single polygon
	  S1sub=S1(i);
    end

    %gather category ids. 1-4 means stable, slumps, noise, gullies.
    [category_idsub]=getcateg(S1sub);
    if length(S1sub) ==0; return;end
    [category_idsub_i]=getcateg(S1(i));

   % % Get the elevation change data, DEM, orthoimage in the target zone.
    %data=readGeotiff(changefile,'map_subset',rangi);
    idx=find(data0.x>=rangi(1)&data0.x<=rangi(2));
    idy=find(data0.y>=rangi(3)&data0.y<=rangi(4));
    if flagdata==1 %1 use the merged elevation change file;
      data.x=data0.x(idx);data.y=data0.y(idy);
      data.z=data0.z(idy,idx);
    elseif flagdata==2
     %fprintf(['\n Working on i polygon:',num2str(i),'/',num2str(n),'; category:',num2str(category_idsub_i),'.\n'])
     fprintf(['\n Working on i polygon:',num2str(i),'; category:',num2str(category_idsub_i),'.\n'])
     %Given a box (rangi), to get the mosaicked files
tic
     [datao]=box2mosaic(rangi,flagband,f,fdir);
%    [datao]=box2mosaic(rangi,3,f,fdir); %test hi
toc

     data=datao;data.z(:,:,2:end)=[]; %to be compatible with older version
     %get orthoimage and DEM	
     dem.x=datao.x;dem.y=datao.y;ortho=dem;delev=dem; %initialize
     dem.z=datao.z(:,:,2);ortho.z=datao.z(:,:,3);delev.z=datao.z(:,:,1);

     %make sure to set void data to NaNs;
     M1=delev.z==0;delev.z(M1)=nan;
     M1=dem.z==-9999; dem.z(M1)=nan;
     M1=ortho.z==0; ortho.z(M1)=nan;
     data.z=delev.z;

     M1=isnan(data.z);
     ratio=sum(M1(:))/length(data.z(:))*100; %ratio of void pixels
     if ratio>10 % more than 10% of the box has no data.
	%fprintf(['\n There is no elevation change data for this slump! \n '])
	fprintf(['\n There are more than 10% of the box has no elevation change data. Skip this slump! \n '])
	return	
     end

    end %flagdata==1

    %Get the matrix of classes of all pixels. 1-4 means stable, slumps, noise, gullies
    [mi,ni]=size(data.z); minm=min([mi,ni]);
    if (minm<3);
	return	
    else
        %count=count+1;
	%countaug=countaug+1;

	%save orthoimage in geotiff 
	if 0 % flagorthotif==1
%	filenamet1i=[odirt1,'/','slump',num2str(count),'.tif'];
%	projstr='polar stereo north';
%	writeGeotiff(filenamet1i,datao.x,datao.y,uint16(ortho.z),2,0,projstr)
	end

        %Get the matrix of classes of all pixels. 1-4 means stable, slumps, noise, gullies
        Mclass=zeros(size(data.z));
    
        resrx=mean(diff(data.x));resry=mean(diff(data.y));
        if 0 %single polygon
            [s1x,s1y]=polarstereo_fwd(S1(i).Y,S1(i).X,[],[],70,-45);
            Xb=round((s1x-data.x(1))/resrx)+1; 
            Yb=round((s1y-data.y(1))/(resry))+1; %descending        
            XYbi1=[Xb(:),Yb(:)]; %n by 2
            id1=find(any(isnan(XYbi1),2));
            XYbi1(id1,:)=[];
        else %multiple polygons
            XYbi=[]; %cell
	    category_idg=[];% array
            nsub=length(S1sub);idd=[];
	    countmp=0; %count of multiple polygons

            for i1=1:nsub
            [s1x,s1y]=polarstereo_fwd(S1sub(i1).Y,S1sub(i1).X,[],[],70,-45);
            Xb=round((s1x-data.x(1))/resrx)+1; 
            Yb=round((s1y-data.y(1))/(resry))+1; %descending        
            XYbi1=[Xb(:),Yb(:)]; %n by 2
            id1=find(any(isnan(XYbi1),2));
            XYbi1(id1,:)=[]; 
            %XYbi{i1}=XYbi1;
         
            %get rid of the out of boundary points.
            [nyi,nxi,nzi]=size(data.z);
            maskorg=false(size(nyi,nxi));
        
            maskorg_k=poly2mask(XYbi1(:,1),XYbi1(:,2), nyi,nxi); 
            maskorg=maskorg|maskorg_k;

            Mclass(maskorg)=category_idsub(i1);
	    
	    end %i1
        end
	mask_stable=(delev.z>=-1&delev.z<=0.5)&~(Mclass>=2); %>= -1 m & <=0.5m. 
	smlarea=1e3*2; %1000 m^2
	mask_stable = bwareaopen(mask_stable, round(abs(smlarea/resrx/resry))); %remove small clusters
	Mclass(mask_stable)=1;

	%Get feature table
	if scenario ==2
	   %get the 100 m buffer zone of stable area around the targets;

            %find out the dominating class of this window.
            classi=category_idsub_i;

	    %get mask of the target cluster
            [s1x,s1y]=polarstereo_fwd(S1(i).Y,S1(i).X,[],[],70,-45);
            Xb=round((s1x-data.x(1))/resrx)+1; 
            Yb=round((s1y-data.y(1))/(resry))+1; %descending        
            XYbi1=[Xb(:),Yb(:)]; %n by 2
            id1=find(any(isnan(XYbi1),2));
            XYbi1(id1,:)=[];
            %get rid of the out of boundary points.
            [nyi,nxi,nzi]=size(data.z);
            maskorg=false(size(nyi,nxi));
            maskorg_k=poly2mask(XYbi1(:,1),XYbi1(:,2), nyi,nxi); 
            maskorg=maskorg|maskorg_k;

	    % add 50m buffer zone to get the stable background zone;retreat rate of around 10 m /yr (Melissa Jones, 2016);
	    widpix=round(50/resrx); %
 	    ptmask= imdilate(maskorg, ones(widpix*2)); 
	    widpix2=widpix+9+2; %9 pixels band is the background;
 	    ptmask2= imdilate(maskorg, ones(widpix2*2)); 
	    maskbuff=ptmask2&~ptmask;
	    maskback=maskbuff&mask_stable; %background
	    
	    if (sum(maskback(:))<100);return;end %not enough background pixels; skip this polygon.

	    %jump.mask if -1, the pixel is the target cluster; if 1, the pixe is the background zone, 0 means void.
	    delev.mask=zeros(size(delev.z));
	    delev.mask(maskorg)=-1;delev.mask(maskback)=1;

               %countf=countf+1;
		countf=1;
               [feati,featlabels]=extractfeature(dem,ortho,delev,scenario);
               featgt(:,countf)=feati;
               category_id_feat(countf)=classi;

	elseif scenario==1
	%divide the grid to 9 by 9 small windows;
tic
	nw=9; 
	II=1:nxi;JJ=1:nyi;
	ny1=length([1:nw:nyi]);nx1=length([1:nw:nxi]);
	featgt_par=cell(ny1,nx1); category_id_feat_par=zeros(ny1,nx1);
%	poolsize=20;
%	poolobj=parpool(20);
	iic=1:nw:nxi;jjc=1:nw:nyi;
	%parfor ii1=1:nx1 %1:nw:nxi;
	for ii1=1:nx1 %1:nw:nxi;
	for jj1=1:ny1  %1:nw:nyi;
	    ii=iic(ii1);jj=jjc(jj1);
	    Mx=II>=ii&II<=ii+nw-1;My=JJ>=jj&JJ<=jj+nw-1;
	    dem_win=[];delev_win=[];ortho_win=[];Mclass_win=[];
	    dem_win.x=dem.x(Mx);dem_win.y=dem.y(My);dem_win.z=dem.z(My,Mx);
	    delev_win=dem_win; ortho_win=dem_win;
	    delev_win.z=delev.z(My,Mx); ortho_win.z=ortho.z(My,Mx);
	    Mclass_win=dem_win;Mclass_win.z=Mclass(My,Mx);
	
	    nw_total=length(dem_win.z(:));
	    if (nw_total<0.8*nw*nw); continue;end % if the window size less than 9by 9, skip

	    %find out the dominating class of this window.
	    classi=mode(Mclass_win.z(:));
	    Mi=Mclass_win.z==classi;

	    n_good=sum(sum(Mi&(~isnan(delev_win.z))));

	    if n_good/nw_total>0.8 
	       [feati,featlabels]=extractfeature(dem_win,ortho_win,delev_win,scenario);	
	        
	       if 0 %not parallel
  	     %   countf=countf+1;
	%	featgt(:,countf)=feati;
	%	category_id_feat(countf)=classi;
	       else %parallel
		featgt_par{jj1,ii1}=feati;
	        category_id_feat_par(jj1,ii1)=classi;
	       end
	    end
	    
	end %jj
	end %ii
%	delete(poolobj)

%	collect data
	countf=0;
	for ii1=1:nx1 %1:nw:nxi;
	for jj1=1:ny1  %1:nw:nyi;
        if ~isempty(featgt_par{jj1,ii1})
        countf=countf+1;
        featgt(:,countf)=featgt_par{jj1,ii1};
	category_id_feat(countf)=category_id_feat_par(jj1,ii1);
        end
	end %jj
	end %ii
toc

	end %scenario

	%save test_svm2_t1.mat -v7.3
        

    end


return
end
