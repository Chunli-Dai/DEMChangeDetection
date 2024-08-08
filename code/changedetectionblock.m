function [jumpt2,jumpstdt2,ratet2,ratestdt2,timec2,timestd2,novlpf2,oflagc2]=changedetectionblock(isel,blocksize,xout,yout,dX4Sg,idregion,XYbg,f,fdir,xq1,yq1)
    %processing change detection (load files and time series analysis) on each block.
    %small block needs small memory size; big block requires large memory size;
    constant
    
    if flagoutput == 1
    %fid3= fopen('trenddf.txt', 'a');
    end
    yr=365.25;
    nq1=length(xq1);

    nsuby=length(yout);nsubx=length(xout);
    nblocky=length(1:blocksize:nsuby);
%     nblockx=length(1:blocksize:nsubx);

    iblock=ceil(isel/nblocky);jblock=isel-(iblock-1)*nblocky;
    %mx0 my0 index of xout yout within this block.
    mx0=(iblock-1)*blocksize+1:min([nsubx,iblock*blocksize]); 
    my0=(jblock-1)*blocksize+1:min([nsuby,jblock*blocksize]); 
    rang0st=[min(xout(mx0)) max(xout(mx0)) min(yout(my0)) max(yout(my0))]; %range of this block
    [MX0,MY0]=meshgrid(mx0,my0);
    idop=(MX0-1)*nsuby+MY0; %find(iselop==isel); % id of output pixels for this block.
    
    fprintf(['\n\n Processing block ',num2str(isel),'; [jblock,iblock]:',num2str([jblock,iblock]),'.'])

    nyj=length(my0);nxj=length(mx0);
    
    %initialize output vairables
    jumpt2=zeros(nyj,nxj);jumpstdt2=zeros(nyj,nxj);timec2=zeros(nyj,nxj);
    ratet2=zeros(nyj,nxj);ratestdt2=zeros(nyj,nxj);
if algorithmin==1 %ice melting
    timestd2=zeros(nyj,nxj,2);
else
    timestd2=zeros(nyj,nxj);
end
    novlpf2=zeros(nyj,nxj,'int32');%maximum repeat of the chosen group at a pixel
    oflagc2=zeros(nyj,nxj,'int32');
    
    %find the strip DEM within this block.  
    %coordinates of the block;
    wm=[];%temporary variable
    wm.x=xout(mx0);wm.y=yout(my0);nwy=length(wm.y);nwx=length(wm.x);
    Mcb=true(nwy,nwx);idblock=[];
    for j=1:length(idregion)
        i=idregion(j);
 % get the Polygon boundary for actual data
        demdir=fdir{i};
        infile= [demdir,'/',f{i}];
        %[XYbi,rangei]=imagebd(infile);
        XYbi=XYbg{i};
        Xb=XYbi(:,1);Yb=XYbi(:,2);
        
        % check whether this polygon intersect with the block
        idx=round((Xb-wm.x(1))/resr)+1;
        idy=round((Yb-wm.y(1))/(-resr))+1;
        Mb = poly2mask(idx,idy, nwy,nwx); % build polygon mask       
        overl=Mb&Mcb;
        %it's possible to work on one pixel block.
        if(sum(sum(overl))>0);idblock=[idblock;j]; %idblock: the index of idregion, dX4Sg
    	    fprintf(['\n Overlapping Files : ',infile])
        end
    end
    id=idregion(idblock); %selected files within this block

    %sort idblock based on time
    nid=length(id);
    t=zeros(nid,1); 
    for j=1:nid %0%length(id)
        filename=f{id(j)};
%       ymd=filename2ymd(filename);
        infile= [fdir{id(j)},'/',f{id(j)}];
        [ymd,~]=strip2date(infile,3);
        i=id(j); 
        %t(j)=datenum(ymd,'yyyymmdd');
        t(j)=datenum(ymd,'yyyymmddHHMMSS');
    end
    [~,idsort]=sort(t);idblock=idblock(idsort); %sort the id based on time
    if nid<2 || length(idop(:))<1; return;end %

%    rang0st=; %range of this block

    x0st=[rang0st(1) rang0st(2) rang0st(2) rang0st(1) rang0st(1) ];y0st=[rang0st(4) rang0st(4) rang0st(3) rang0st(3) rang0st(4) ];
    exb=max([100.,2*resr]); %expand the block by 100 m for reading DEMs.
    rang0exb=[rang0st(1)-exb rang0st(2)+exb rang0st(3)-exb rang0st(4)+exb ]  % boundary for reading DEM

%find the overlapping zone for this piece
% mx0=find(xout>=rang0st(1) & xout<=rang0st(2) );
% my0=find(yout>=rang0st(3) & yout<=rang0st(4) ); %should be the same as above.
% nyj=length(my0);nxj=length(mx0);
% demg=nan*ones(nyj,nxj,length(id));
demg=zeros(nyj,nxj,nid);
datadsx=xout(mx0);datadsy=yout(my0);

[X,Y]=meshgrid(datadsx,datadsy);
[LAT,LON]=polarstereo_inv(X,Y,[], [],70,-45);

tic
t=zeros(nid,1);%clear demg demgmt
demg=-9999*ones(nyj,nxj,nid);demgmt=zeros(nyj,nxj,nid);
idd=[];

%loading DEM and apply offsets.
for j=1:nid %[3,5,24,26] % 2012 Kamchatka Volcano 
        i=idregion(idblock(j));
        demdir=fdir{i};
        infile= [demdir,'/',f{i}];
    	fprintf(['\n Loading DEM files in chronological order (',num2str(j),'/',num2str(nid),', isel ',num2str(isel),'): ',infile])
	%ymd=filename2ymd(f{i});
        [ymd,~]=strip2date(infile,3);
	%t(j)=datenum(ymd,'yyyymmdd');
        t(j)=datenum(ymd,'yyyymmddHHMMSS');
        XYbi=XYbg{i};

        if strcmp(infile(end-7:end),'meta.txt')
          str1='meta.txt';
        elseif strcmp(infile(end-7:end),'_mdf.txt')
          str1='mdf.txt';
        end
	infile1= strrep(infile,str1,demext);
%       clear data %clear violates workspace transparency for parallel computation.

	if flagfilter==1 %apply filter of water and clouds
           [data,~]=readdem(XYbi,infile,rang0exb,resr);
	else
            data=readGeotiff(infile1,'map_subset',rang0exb);
	end

        %Interpolation requires at least two sample points in each dimension.
	[m1,n1]=size(data.z);
	if min([m1,n1])<=1;continue;end

	%refers to ~/arcticdemapp/river/rivergithub2/rivermainserial/riverprofsub.m
%filtering out bad edges, e.g. 20110804
% Me=~(imdilate((data.z==-9999),ones(round(30*8)))); %Bug 23, remove too much river area
%data.z(~Me)=-9999; %April 2019, pzxy may change the value of -9999.
%edge based on boundary files.
	dx=data.x(2)-data.x(1); dy=data.y(2)-data.y(1);[nwy,nwx]=size(data.z);
	sx=XYbg{i}(:,1); sy=XYbg{i}(:,2); 
	idx=round((sx-data.x(1))/dx)+1;idy=round((sy-data.y(1))/(dy))+1;
	Med=~poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
	Me=~(imdilate(Med,ones(round(30*8*2/resr)))); 
	data.z(~Me|data.z==-9999)=NaN;

        p=dX4Sg(idblock(j),:); %pzxy

    	fprintf(['\n Loading its translational offsets pzxy (',num2str(j),'/',num2str(nid),', isel ',num2str(isel),'): ',num2str(p)])
	data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
    	z2inp=interp2(data.x-p(2),data.y-p(3),double(data.z)-p(1),datadsx,datadsy','*linear');
    	z2inp(isnan(z2inp))=-9999; %return to -9999
        demg(1:nyj,1:nxj,j)=z2inp;

        %matchtag file as weight
	infile1= strrep(infile,str1,'matchtag.tif');
        if exist(infile1,'file')
        data=readGeotiff(infile1,'map_subset', rang0exb);
        else
        data.z=ones(size(data.z));
        end

        % interpolate to match the coregistered DEM
        mt = interp2(data.x-p(2),data.y-p(3),data.z,datadsx,datadsy','*nearest',0);    
        demgmt(1:nyj,1:nxj,j)=mt;%mt(idy,idx);
end %

toc  % loading 34 files take 8 minutes, too slow.
% coregistering 17 files take 29 minutes.

fprintf ('\n Step 2.3: Time series analysis.')
    %Preparing the time series analysis
    epochorg=zeros(length(t),1);
    epochorg(:)=t(:);
    algorithm=2; % 1 linear; 2 constant
    epoch=epochorg/yr;  % Change the unit of time from day to year, to avoid singular matrix
    tm=mean(epoch);
    
    [~,~,ni]=size(demg);
    demg(isnan(demg))=-9999;
    % Time series analysis
    tic
%   demp=zeros(ni,1);dempmt=ones(ni,1);
    AMa=ones(ni,1);
    stdmax=40;
%     jx=find(abs(datadsx(idx)-xeq)<1);jy=find(abs(datadsy(idy)-yeq)<1);
%   jxq1=find(abs(datadsx-xq1)<resrc/1.5);jyq1=find(abs(datadsy-yq1)<resrc/1.5);
%     jumpt2=zeros(nyj,nxj);jumpstdt2=zeros(nyj,nxj);timec2=zeros(nyj,nxj);
    
%   poolobj=parpool(poolsize); %need to change the following to a function for parallel.
    for jxy=1:length(idop(:))
            mxp=ceil(idop(jxy)/nsuby);myp=idop(jxy)-(mxp-1)*nsuby; %mxp, myp are index of xout, yout; idop is the 1D index of jump.
            jy=find(my0==myp);jx=find(mx0==mxp); %jy jx are index of demg.
            if length(jy)~=1 ||length(jx)~=1 ; warning(['\n Wrong index at jxy=',num2str(jxy),'!']); end
            xp=xout(mxp);yp=yout(myp); %meter

            demp=zeros(ni,1);dempmt=ones(ni,1);
            demp(:)=demg(jy,jx,:); % DEM at a point; p denotes point; g denotes group;
            dempmt(:)=demgmt(jy,jx,:);    %match tag

            %outlier dtection
%             [idkp,flagcond]=outlier(demp,dempmt,epochorg);
            
%             idkp= idkp(~ismember(idkp,[6,13,21,23,25]));
            
	    %only one requested time series
            if 0
            if ismember(jx,jxq1)&&ismember(jy,jyq1); 
                flagplotsv=1;
                jxy
            else
                flagplotsv=0;
            end %0 %plotting
            end

	    %an array of requested locations for plotting time series
	    flagxp=0;
            for jq=1:nq1
                jxq1=find(abs(datadsx-xq1(jq))<=resr/2);jyq1=find(abs(datadsy-yq1(jq))<=resr/2);
                if ismember(jx,jxq1)&&ismember(jy,jyq1); %jx=jxq1(1), jy=jyq1(2); not good
                   fprintf(['\n xq1(jq) with jq:',num2str(jq),'; [lenx leny]=',num2str([length(jxq1) length(jyq1)]),'; isel:',num2str(isel)])
                   flagxp=1;break;
                end
            end

            if flagxp==1
                flagplotsv=1;
                jxy
            else
                flagplotsv=0;
                jq=0;
            end %0 %plotting


            [oflag,trend,trest,rate,ratestd,eqm,eqs,eqe]=timeseries2(demp,dempmt,epochorg,timefix,eqepoch,flagplotsv,jq);

            
     if flagoutput == 1
        %   fprintf(fid3,'%12.6f %12.6f  %23.15e  %23.15e  %12.6f %s %s\n',LAT(jy,jx),LON(jy,jx),trend,trest,eqm,eqs,eqe);

	%e.g.  1 linear (ice melting); 2 constant (landslides) ; 3 constant + linear (Okmok volcano)
 	   if algorithmin==1
           % fprintf(fid3,'%12.6f %12.6f  %23.15e  %23.15e  %12.6f %s %s\n',LAT(jy,jx),LON(jy,jx),rate,ratestd,eqm,eqs,eqe);
 	   % to save it to jump;
 	    % trend=rate;trest=ratestd;
 	   elseif algorithmin==2
           % fprintf(fid3,'%12.6f %12.6f  %23.15e  %23.15e  %12.6f %s %s\n',LAT(jy,jx),LON(jy,jx),trend,trest,eqm,eqs,eqe);
 	   elseif algorithmin==3
           % fprintf(fid3,'%12.6f %12.6f  %23.15e %23.15e   %23.15e %23.15e\n',LAT(jy,jx),LON(jy,jx),trend,trest,rate,ratestd);
 	   end

     end 


%         if jumpflag==1
	      % m/yr
%         if (trest< 0.8)
%           fprintf(fid2,'%12.6f %12.6f  %23.15e %23.15e\n',LAT(jy,jx),LON(jy,jx),trend,trest);
%         end %if
            eqstd=(datenum(eqe,'yyyy/mm/dd')-datenum(eqs,'yyyy/mm/dd'))/2./yr;
            
            j1=idop(jxy);
%           jump(idop(jxy))=trend;
%             jump(j1)=trend;
%             jumpstd(idop(jxy))=trest;  
%             timec(idop(jxy))=eqm;  
%             oflagc(idop(jxy))=int32(oflag);
%             timestd(idop(jxy))=eqstd;
%             novlpf(idop(jxy))=sum((demp~=-9999));
            
            jumpt2(jy,jx)=trend;
            jumpstdt2(jy,jx)=trest; 
            ratet2(jy,jx)=rate; 
            ratestdt2(jy,jx)=ratestd;
            timec2(jy,jx)=eqm;
if algorithmin==1 %ice melting
	   timestd2(jy,jx,1)=datenum(eqs,'yyyy/mm/dd')/yr;
	   timestd2(jy,jx,2)=datenum(eqe,'yyyy/mm/dd')/yr;
else
            timestd2(jy,jx)=eqstd; %yr
	end
            novlpf2(jy,jx)=sum((demp~=-9999));
            oflagc2(jy,jx)=int32(oflag);

%         end

    end %jxy
    %delete(poolobj)

     if flagoutput == 1
      %    fclose(fid3)
     end 

    toc
    %End of Time series analysis
%           jump(:,:,isel)=jumpt2;
%           fprintf(fid2,'%12.6f %12.6f   %23.15e %23.15e\n',LAT(jy,jx),LON(jy,jx),trend,trest);
     if 0
        figure %(100)
        hold all
        xt=xout(mx0)*1e-3;yt=yout(my0)*1e-3;zt=jumpt2;
    %     zt(mp1)=NaN;
        imagesc(xt,yt,zt); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
        set(gca,'FontSize', 18);
        plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
        xlabel('{\itx} (km)')
        ylabel('y (km)')
        caxis([-50 50])
        plot(x0*1e-3, y0*1e-3,'g-','linewidth',4)
        plot(x0st*1e-3, y0st*1e-3,'r-','linewidth',6)
        title(['Jump (m); isel=',num2str(isel)])
%       x0sov=[rang0sov(1) rang0sov(2) rang0sov(2) rang0sov(1) rang0sov(1) ];
%       y0sov=[rang0sov(4) rang0sov(4) rang0sov(3) rang0sov(3) rang0sov(4) ];
%       plot(x0sov*1e-3, y0sov*1e-3,'m-','linewidth',6)
        title('Elevation change')

% saveas(gcf,'jumpisel7good9','fig')

        figure %(100)
        hold all
        xt=xout(mx0)*1e-3;yt=yout(my0)*1e-3;zt=jumpt2;
    %     zt(mp1)=NaN;
        imagesc(xt,yt,timec2); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
        set(gca,'FontSize', 18);
        plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
        xlabel('{\itx} (km)')
        ylabel('y (km)')
        caxis([0 18])
        plot(x0*1e-3, y0*1e-3,'g-','linewidth',4)
        plot(x0st*1e-3, y0st*1e-3,'r-','linewidth',6)
        title(['Event time (yr) after 2000; isel=',num2str(isel)])
%       x0sov=[rang0sov(1) rang0sov(2) rang0sov(2) rang0sov(1) rang0sov(1) ];
%       y0sov=[rang0sov(4) rang0sov(4) rang0sov(3) rang0sov(3) rang0sov(4) ];
%       plot(x0sov*1e-3, y0sov*1e-3,'m-','linewidth',6)

     end %if plot
%    save test5.mat -v7.3

end %function
