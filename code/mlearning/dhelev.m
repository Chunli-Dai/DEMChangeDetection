function [co]=dhelev(jump,dem,Mouts,S2,Mn,Mp)
% input: jump: elevation change data;
%        dem: DEM data
%        S2: shapefile of truth slumps;
%        Mouts: mask of the clusters to be checked; 
%               it is a structure with different matrix size for each cluster.
% reference to /Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/poly2precision.m
% ex: [co]=dhelev(jump,demr,Mouts,S2);
% Results: flagout=1: dh vs elev plot
%          flagout=2: the mask of good clusters of both negative and positive clusters.    

IoUthre=0.4; %0.4;  0.5 0.6  
buff=1e3; %assume the ground truth polygon that overlaps with the target is no wider than 1 km
resr=10; %2m grid size to creat mask from given polygon
resrmask=2; % 66sec (10m) vs 64 sec (2m)
flaglat=1; %1, shapefile cooridnates are in longitude/latitude; 0, in polar stereographic coordinates
co=[];
flagout=2; %1, calculate  for S2 ; 2 calculate for Mn. 

if flaglat==1
    fprintf('\n Assume groud truth shapefiles is in longitude/latitude. \n')
elseif flaglat==0
    fprintf('\n Assume groud truth shapefiles is in polar stereographic coordinates. \n')
end

M=Mn.z;
mask1=jump;
mask1.z=M;

% BW=M;
% CC = bwconncomp(BW);
% TPFP=CC.NumObjects;

if flagout==1
    
lenS2=length(S2);
for j=1:lenS2 %ground truth
        if flaglat==1
        [s2x,s2y]=polarstereo_fwd(S2(j).Y,S2(j).X,[],[],70,-45);
        else
        s2x=S2(j).X;s2y=S2(j).Y;
        end
        s2x(isnan(s2x))=[];s2y(isnan(s2y))=[];
        s2xg{j}=s2x;s2yg{j}=s2y;
end
   
FP=0;TP=0;FN=0;nout=0;
idp=[]; %ids of all ground truth slumps that paired with target slumps.
s2xp=[];s2yp=[]; %for plot
maskj_all=false(size(jump.z));
    for j=1:lenS2 %ground truth
        s2x=s2xg{j};s2y=s2yg{j};
        resr1m=1; %for small temporary areas; better resolution yields better mask from polygon 
        rang0=round([min(s2x)-buff max(s2x)+buff min(s2y)-buff max(s2y)+buff]/resr1m)*resr1m;
        maskx=rang0(1):resr1m:rang0(2);masky=rang0(4):-resr1m:rang0(3);
        ny=length(masky);nx=length(maskx);
        x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
        
        Xb1=round((s2x-maskx(1))/resr1m)+1; 
        Yb1=round((s2y-masky(1))/(-resr1m))+1; %descending    
        maskj=poly2mask(Xb1,Yb1, ny,nx);  %ground truth  
               
        %crop mask1 for faster interpolation
        rang2=round([min(s2x)-2*buff max(s2x)+2*buff min(s2y)-2*buff max(s2y)+2*buff]/resr1m)*resr1m;
        Mx=mask1.x>=rang2(1)&mask1.x<=rang2(2);My=mask1.y>=rang2(3)&mask1.y<=rang2(4);
        mask1r.x=mask1.x(Mx);mask1r.y=mask1.y(My);
        mask1r.z=mask1.z(My,Mx); % %input mapped area
        %crop the matrix of positive clusters 
        Mpc=mask1r;Mpc.z=Mp.z(My,Mx);
        
        minm=min(size(mask1r.z));
        if minm<3 %ground truth range is out of the boundary of mask1
            nout=nout+1;
            continue
        end
                
        jumpb=mask1r;jumpb.z=jump.z(My,Mx); 
        demb=mask1r;demb.z=dem.z(My,Mx); 
        jumpr.x=maskx;jumpr.y=masky;
        jumpr.z=interp2(jumpb.x,jumpb.y,jumpb.z,maskx,masky','*nearest',0);
        demr=jumpr;
        demr.z=interp2(demb.x,demb.y,demb.z,maskx,masky','*nearest',0);
        
        if 0
            figure;imagesc(jumpr.x*1e-3,jumpr.y*1e-3,jumpr.z,'alphadata',Mout.z);colorbar;colormap jet;caxis([-5 5])
            hold on; plot(s2x*1e-3,s2y*1e-3,'k.-')
            set(gcf,'Color','white');
            xlabel('{\itx} (km)');ylabel('{\ity} (km)');title(['Ground truth slump ',num2str(j)])
%             
%             figure;hillshade(demr.z,demr.x,demr.y,'plotit');colorbar;
%             hold on; plot(s2x,s2y,'k.-')
%             set(gcf,'Color','white');
%             xlabel('{\itx} (m)');ylabel('{\ity} (m)');
        end
        
        M1=interp2(Mpc.x,Mpc.y,single(Mpc.z),maskx,masky','*nearest',0);
        Mpc2.x=maskx;Mpc2.y=masky;Mpc2.z=logical(M1);
        maskjm.x=maskx;maskjm.y=masky;maskjm.z=maskj;
        %find the matching positive clusters.
        [Mouts, Mout]=matchclusters(maskjm,Mpc2);
        
        [co1]=dhelev_sub(jumpr,demr,logical(Mout.z));   
        
        %store ground truth mask
%         tic; maskj_all1=interp2(maskx,masky,int8(maskj),jump.x,jump.y','*nearest',0);  toc   %%input mapped area
        tic;maskjc=interp2(maskx,masky,int8(maskj),mask1r.x,mask1r.y','*nearest',0);  
        maskj_all1=false(size(jump.z));maskj_all1(My,Mx)=maskjc; toc
        maskj_all=logical(maskj_all)|logical(maskj_all1);
                
    end %for j
    
    j=1;
    
        
elseif flagout==2 
    % loop for all False positives
%     MFP=M&~(maskj_all);
    MFP=M; %all clusters
    maskcri1to3=false(size(M));
    
    smlarea=1e3*2; %1000 m^2
resx=mean(jump.x(2:end)-jump.x(1:end-1));resy=mean(jump.y(2:end)-jump.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);
    MFP = bwareaopen(MFP, smlarea/resr/resr); %remove small clusters

    [X,Y]=meshgrid(jump.x,jump.y);
    
    CC = bwconncomp(MFP);
    
    for ki=1:1:CC.NumObjects
        fprintf(['\n Working on ki :',num2str(ki),'/',num2str(CC.NumObjects)])
        BW3=false(size(MFP));
        BW3(CC.PixelIdxList{ki})=MFP(CC.PixelIdxList{ki});
        xmin=min(X(BW3));xmax=max(X(BW3));ymin=min(Y(BW3));ymax=max(Y(BW3));
        rang2=[xmin-buff xmax+buff ymin-buff ymax+buff];
        Mx=mask1.x>=rang2(1)&mask1.x<=rang2(2);My=mask1.y>=rang2(3)&mask1.y<=rang2(4);
        
        %crop the area around the target cluster with negative signal 
        mask1r.x=mask1.x(Mx);mask1r.y=mask1.y(My);
        mask1r.z=BW3(My,Mx); % %input mapped area
         %crop the matrix of positive clusters 
        Mpc=mask1r;Mpc.z=Mp.z(My,Mx);
       
        %find the matching positive clusters. Mout is the matched pairs of negative and positive clusters
        [Mouts, Mout]=matchclusters(mask1r,Mpc);
        
        jumpb=mask1r;jumpb.z=jump.z(My,Mx); %elevation change in the subzone.
        demb=mask1r;demb.z=dem.z(My,Mx); demb.z(demb.z==-9999|demb.z==0)=nan;
        [flagi]=dhelev_sub(jumpb,demb,Mout.z);
        
        if flagi==1 %good
            maskgood=interp2(Mout.x,Mout.y,double(Mout.z),jump.x,jump.y','nearest',0);
            maskcri1to3=maskcri1to3|logical(maskgood);
            
        else %don't save data
            %bad
        end
        
        if 0
        figure;imagesc(jumpb.x*1e-3,jumpb.y*1e-3,jumpb.z);colorbar;colormap jet; caxis([-3 3]);title('Elevation change (m)')
        set(gcf,'Color','white');xlabel('{\itx} (km)');ylabel('{\ity} (km)');
        figure;imagesc(jumpb.x*1e-3,jumpb.y*1e-3,jumpb.z,'alphadata',Mout.z );colorbar;colormap jet; caxis([-3 3]);title('Elevation change at the matched clusters')
        figure;imagesc(jumpb.x*1e-3,jumpb.y*1e-3,demb.z);colorbar;colormap jet;title('DEM (m)')
%         close all
        end
        
%         pause
    end
    co=maskcri1to3;
end

    
return
end

