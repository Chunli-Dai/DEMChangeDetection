function [Mpc,idsi]=dividemapsub(ix,S1,buff,data0,flagrescale,flagband,f,fdir,odir1)

%   odirt1=['slumps_ortho/'];
    odirt1=strrep(odir1,'slumpsub_','slumptif_');

    resrx=mean(diff(data0.x));resry=mean(diff(data0.y));

    ixy=ix; i=ix;

    % see prepareDetectron2.m changemap2polygon.m
    %get the box in polar coordinates
    ll=S1(ix).BoundingBox(1,:);ru=S1(ix).BoundingBox(2,:);
    [llx,lly]=polarstereo_fwd(ll(2),ll(1),[],[],70,-45); %lat lon
    [rux,ruy]=polarstereo_fwd(ru(2),ru(1),[],[],70,-45); %lat lon
    [rlx,rly]=polarstereo_fwd(ll(2),ru(1),[],[],70,-45); %right low
    [lux,luy]=polarstereo_fwd(ru(2),ll(1),[],[],70,-45); %left up
    x=[llx,rlx,rux,lux];y=[lly,rly,ruy,luy];
    rang2=[min(x)-buff max(x)+buff min(y)-buff max(y)+buff];

    Mx=data0.x>=rang2(1)&data0.x<=rang2(2);My=data0.y>=rang2(3)&data0.y<=rang2(4);
    
    %cropped clusters
    Mpc=[];
    Mpc.x=data0.x(Mx);Mpc.y=data0.y(My);
%   Mpcpre.z=data0.z(My,Mx); %just to get the size.

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
%   data0.z(My,Mx)=Mpc.z|Mpcpre.z;
    Mpc.Mx=Mx;Mpc.My=My;

    %get the test images
    % resolution of data0 will affect rang2 and the matrix dimension of Mpc.
    rang2=[min(Mpc.x) max(Mpc.x) min(Mpc.y) max(Mpc.y)]; %make sure the boundary is consistent with ids.

%idxs=[]; idxe=[]; idys=[]; idye=[];
	idxs=find(Mx,1,'first');idxe=find(Mx,1,'last');
	idys=find(My,1,'first');idye=find(My,1,'last');

        %ids{ix}=[idxs idxe idys idye];
	idsi=[idxs idxe idys idye];

        rangi=rang2;

%fprintf('hi1')
        filenamet1i=[odirt1,'/','slump',num2str(ixy),'.tif'];
	flagfile=0; %file does not exist;
	if exist(filenamet1i,'file');
	   fprintf(['\n File exists, read the data :', filenamet1i])
	   try
	     datao=readGeotiff(filenamet1i);
	     flagfile=1; %file exist;
	   catch e
   	     fprintf('\n dividemapsub.m There was an error! The message was:\n%s',e.message);
	     flagfile=0; %file does not exist;
	   end
	end

	if flagfile==0;
	   try
 	   [datao]=box2mosaic(rangi,flagband,f,fdir,ix); %test hi
	   catch e
           fprintf('\n dividemapsub.m There was an error! The message was:\n%s',e.message);
	   % initialize datao
	   resr2m=2;
           xout=rangi(1):resr2m:rangi(2);yout=rangi(4):-resr2m:rangi(3);
           nx0=length(xout);ny0=length(yout);
           datao.x=xout;datao.y=yout;datao.z=nan(ny0,nx0,3);
           end
	end
	%datao change, curvature, time
	mask1=datao.z(:,:,1) <= -1;
	mask1r=interp2(datao.x,datao.y,double(mask1),Mpc.x,Mpc.y','nearest',0);
	%recover the cluster that's selected originally in changemap2polygon.m
	BW=mask1r;
	Mf4cleankeep=false(size(BW));
	CC = bwconncomp(BW);
	for k=1:CC.NumObjects
	    BW3=false(size(BW));
	    BW3(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});
	    
	    ratio=sum(sum(BW3.*Mpc.z))./sum(sum(Mpc.z))*100;
	    
	    if ratio>50
	        Mf4cleankeep=Mf4cleankeep|BW3;
	        break
	    end
	end
	Mpc.z=Mf4cleankeep; %more accurately reflect the selected clusters in changemap2polygon.m
	if flagfile ==1 %all outputs exist.
		return;
	end
%[datao]=box2mosaic(rangi,1,f,fdir); %test hi
%	fprintf('hi2')
%	save testhi1.mat -v7.3

        %save orthoimage in geotiff 
        if 0 %flagorthotif==1
        filenamet1i=[odirt1,'/','slump',num2str(ixy),'.tif'];
        projstr='polar stereo north';
        writeGeotiff(filenamet1i,datao.x,datao.y,uint16(datao.z(:,:,3)),2,0,projstr)
        end
        if 1 %  flagdemtif==1 %hi
        filenamet1i=[odirt1,'/','slump',num2str(ixy),'.tif'];
        projstr='polar stereo north';
        %writeGeotiff(filenamet1i,datao.x,datao.y,int16(datao.z(:,:,2)),2,0,projstr)
        writeGeotiff(filenamet1i,datao.x,datao.y,int16(datao.z),2,0,projstr)
        end

	%A1=uint16(rescale(datao.z(:,:,1),0,65535)); %Detectron2 only works with uint8!
	A1=[];T2=[];
	if flagrescale==0
	    A1=uint8(rescale(datao.z(:,:,1),0,255));
	elseif flagrescale==1
	    inmin=-20;inmax=10; %truncate data to this range [inmin inmax], then rescale.
	    A1=uint8(rescale(datao.z(:,:,1),0,255,'InputMin',inmin,'InputMax',inmax));
	end

	if flagband==2 
          %get curvature
          totalc=datao.z(:,:,2);
          A1(:,:,2)=uint8((totalc(:,:))); %DEM;
          A1(:,:,3)= A1(:,:,1); %  %Data with 2 components not supported for JPEG files. So add 3rd band.
	elseif flagband ==4
	  T2=datao.z(:,:,2);
          A1(:,:,2)=uint8((T2(:,:))); %
          A1(:,:,3)= A1(:,:,1); %
	elseif flagband ==3
	  A1(:,:,2)=uint8((datao.z(:,:,2))); % %curvature
	  A1(:,:,3)=uint8((datao.z(:,:,3))); % %time
	end

        %ofile=['slumpsub_band1changeuint8jpg/slumppeelsubi',num2str(ixy),'.jpg']; 
        ofile=[odir1,'/slumppeelsubi',num2str(ixy),'.jpg']; 
        imwrite(A1,ofile);

return
end
