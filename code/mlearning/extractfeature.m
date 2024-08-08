% extract features of an input image, DEM
function [feat,featlabels]=extractfeature(DEM,orthoimage,jump,scenario)
% input: DEM in 9 by 9 window; orthoimage and elevation change
% no data are NaNs.
% scenario 1: work on 9 by 9 pixels
% scenario 2: work on a negative elevation change cluster and the surrounding 9 pixels (or 100 m) buffer zone (stable background);
% input: jump.mask if -1, the pixel is the target cluster; if 1, the pixel is the buffer zone, 0 means void.

% feat2=array2table(feat,'Variablenames', featlabels);

% all possible features
featlabels={'slopeangle_mad','aspect_mad','Roughness_median','Roughness_mad','DCER1_median','DCER1_mad', ...
    'DCER2_median','DCER2_mad','hills_mad','RLOV_median','RLOV_mad','TPI_median','TPI_mad', ...
    'profc_median','profc_mad','planc_median','planc_mad','Kt_median','Kt_mad','totalc_median','totalc_mad', ... % 21
    'dh_slopeangle_mad','dh_aspect_mad','dh_Roughness_median','dh_Roughness_mad','dh_DCER1_median','dh_DCER1_mad', ...
    'dh_DCER2_median','dh_DCER2_mad',    'dh_RLOV_median','dh_RLOV_mad','dh_TPI_median','dh_TPI_mad', ...
    'dh_profc_median','dh_profc_mad','dh_planc_median','dh_planc_mad','dh_Kt_median','dh_Kt_mad','dh_totalc_median','dh_totalc_mad',... %20
    'dh_median', 'dh_mad', ...    %extra factors, 2
    'ortho_median_rd', 'ortho_mad_rd',...  %orthoimage of the scar area over background area, 2
    'cri1','cri2','cri3','slope','intercept','fitstd'}; % scenario 2. need negative/positive cluster pairs
m=length(featlabels);
feat=nan(1,m); %array matrix, n by m.

resx=mean(DEM.x(2:end)-DEM.x(1:end-1));resy=mean(DEM.y(2:end)-DEM.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);

%Must smooth DEM before calculating curvature!
DEMs.z = imgaussfilt(DEM.z,2); 
% DEMs.x=DEM.x;DEMs.y=DEM.y;
DEM.z=DEMs.z;
jumps.z = imgaussfilt(jump.z,2); jump.z=jumps.z;

if scenario ==1
    
    dem_feat=extractfeature_sub(DEM);
    jump_feat=extractfeature_sub(jump);
    
    dh_median=nanmedian(jump.z(:));
    dh_mad=mad(jump.z(:)); 
    
% feat=[dem_feat.slopeangle_mad dem_feat.aspect_mad dem_feat.Roughness_median dem_feat.Roughness_mad dem_feat.DCER1_median dem_feat.DCER1_mad  ...
%     dem_feat.DCER2_median dem_feat.DCER2_mad dem_feat.hills_mad dem_feat.RLOV_median dem_feat.RLOV_mad dem_feat.TPI_median dem_feat.TPI_mad ...
%     dem_feat.profc_median dem_feat.profc_mad dem_feat.planc_median dem_feat.planc_mad dem_feat.Kt_median dem_feat.Kt_mad dem_feat.totalc_median dem_feat.totalc_mad ...
%     jump_feat.slopeangle_mad jump_feat.aspect_mad jump_feat.Roughness_median jump_feat.Roughness_mad jump_feat.DCER1_median jump_feat.DCER1_mad ...
%     jump_feat.DCER2_median jump_feat.DCER2_mad   jump_feat.RLOV_median jump_feat.RLOV_mad jump_feat.TPI_median jump_feat.TPI_mad ...
%     jump_feat.profc_median jump_feat.profc_mad jump_feat.planc_median jump_feat.planc_mad jump_feat.Kt_median jump_feat.Kt_mad jump_feat.totalc_median jump_feat.totalc_mad ...
%     ];  %NaN NaN NaN NaN NaN NaN];

    feat=[dem_feat(1:21), ...
    jump_feat([1:8,10:21]), ...
    dh_median,dh_mad,...
    ];
    %Extra factors by Chunli

    %DEMdh: Elevation change map direction vector vs DEM direction vector
    if 0 % not working well.
    [Nx_dh,Ny_dh,Nz_dh] =surfnorm(X,Y,jump.z);
    %matrix to n by 3 series
    u=[Nx(:),Ny(:),Nz(:)]';
    v=[Nx_dh(:),Ny_dh(:),Nz_dh(:)]';

    DEMdh_ang=angleofvectors(u,v);
    DEMdh_ang=reshape(DEMdh_ang,size(DEM.z));%vector to matrix
    end

% 
elseif scenario == 2  %target cluster
    %DEMdh2.cri1 DEMdh2.cri2 DEMdh2.cri3 DEMdh2.slope DEMdh2.intercept DEMdh2.fitstd];
    
    %get masks
    %jump.mask if -1, the pixel is the target cluster; if 1, the pixe is the background zone, 0 means void.
    M_tar=jump.mask==-1; %target
    M_bac=jump.mask==1; %background
    
    % check the effect of NaNs on the edge for the calculation of curvature and other factors. 
    jump_tar=jump;DEM_tar=DEM;jump_tar.z(~M_tar)=nan;DEM_tar.z(~M_tar)=nan;
    dem_feat_tar=extractfeature_sub(DEM_tar);
    jump_feat_tar=extractfeature_sub(jump_tar);

    jump_bac=jump;DEM_bac=DEM;jump_bac.z(~M_bac)=nan;DEM_bac.z(~M_bac)=nan;
    dem_feat_bac=extractfeature_sub(DEM_bac);
    jump_feat_bac=extractfeature_sub(jump_bac);
    
    %get relative differences;
    dem_feat_rd=rela_diff(dem_feat_bac,dem_feat_tar); %x1 stable area; x2 target;
    jump_feat_rd=rela_diff(jump_feat_bac,jump_feat_tar);
    
    dh_median=nanmedian(jump_tar.z(:));
    dh_mad=mad(jump_tar.z(:));   %extra factors, 2
    
    ortho_median_rd=rela_diff(nanmedian(orthoimage.z(M_bac)),nanmedian(orthoimage.z(M_tar)));
    ortho_mad_rd= rela_diff(mad(orthoimage.z(M_bac)),mad(orthoimage.z(M_tar)));
    
    [flagi,DEMdh2]=dhelev_sub(jump,DEM,M_tar);
    
    feat=[dem_feat_rd(1:21), ...
    jump_feat_rd([1:8,10:21]), ...
    ...
    dh_median,dh_mad, ...    %extra factors, 2
    ortho_median_rd,ortho_mad_rd, ...  %orthoimage of the scar area over background area.
    ...
    DEMdh2.cri1 DEMdh2.cri2 DEMdh2.cri3 DEMdh2.slope DEMdh2.intercept DEMdh2.fitstd];
    

end

return
end

function [ThetaInDegrees]=angleofvectors(u,v)
%calculate the angles between two vectors
% u can be 3 by 1 or 3 by n;
% v should be the same size as u;

norm_u=(dot(u,u)).^0.5;norm_v=(dot(v,v)).^0.5;

CosTheta = max(min(dot(u,v)./(norm_u.*norm_v),1),-1);
ThetaInDegrees = (acosd(CosTheta)); % [0 180];

return
end



% extract features of an input image, DEM
function [out]=extractfeature_sub(DEM)
% input: 3D surface
% output: some statistics.

% no data are NaNs.
out=nan(1,21);

resx=mean(DEM.x(2:end)-DEM.x(1:end-1));resy=mean(DEM.y(2:end)-DEM.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);

%Features are the median absolute deviation of each factor. 

%Mora 2015, 
%Terrain Slope; see /Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/watermask.m
% get the slope
% [px,py] = gradient(DEM.z,resr); % second parameter is the image resolution 
% py=-py; % minus sign in front of py, since data.z has the y axis positive downwards.
[X,Y]=meshgrid(DEM.x,DEM.y);
[Nx,Ny,Nz] =surfnorm(X,Y,DEM.z);px=Nx;py=Ny; %same as above two commented lines, 
           % because unit surface normal of an implicit surface f(x, y, z)=0 definitely is the normalized gradient of f.
M=isnan(Nx)|isnan(Ny)|isnan(Nz);
Nxt=Nx;Nxt(M)=nan;Nyt=Ny;Nyt(M)=nan;Nzt=Nz;Nzt(M)=nan;%in case of different elements have NaNs, change the same elements to NaN.
[theta,rho] = cart2pol(px,py);% th = atan2(y,x);
%slope angle
slopeangle=atan(rho)*180/pi;% 0 to 90
% out.slopeangle_mad=mad(slopeangle(:),1); %median(abs(X – median(X))).
out(1)=mad(slopeangle(:),1); 

% Aspect
theta=theta*180/pi; %radius to degree
% out.aspect_mad=mad(theta(:),1);
out(2)=mad(theta(:),1);

% Roughness % R = Max(abs(Zij – Z11)) in Mora, 2015
Roughness=nan(size(DEM.z)); TPI=nan(size(DEM.z));
DCER1=nan(size(DEM.z));DCER2=nan(size(DEM.z));RLOV=nan(size(DEM.z));
II=1:length(DEM.x);JJ=1:length(DEM.y);

%Sliding window size of 3 by 3
stepi=1; %10; %1
for i=1:stepi:length(DEM.x)
for j=1:stepi:length(DEM.y)
    Mx=II>=i-1&II<=i+1;My=JJ>=j-1&JJ<=j+1;
    DEMt=DEM.z;
    Dsub=DEMt(My,Mx);
    %convert matrix to vector.
    Nxi=Nxt(My,Mx); Nxi=Nxi(:);Nyi=Nyt(My,Mx); Nyi=Nyi(:);Nzi=Nzt(My,Mx); Nzi=Nzi(:);

    % all 9 pixels must be valid.
    M=isnan(Dsub)|isnan(Nxt(My,Mx));
    if sum(M(:))>=1; continue;end
    
    % The Direction Cosine Of Eigenvalue Ratios DCER
    % Resultant Length of Orientation Vectors RLOV; 
%     T_eigen=[nansum(Nx(~M).^2) nansum(Nx(~M).*Ny(~M)) nansum(Nx(~M).*Nz(~M));
%          nansum(Nx(~M).*Ny(~M)) nansum(Ny(~M).^2) nansum(Ny(~M).*Nz(~M));
%         nansum(Nx(~M).*Nz(~M)) nansum(Ny(~M).*Nz(~M)) nansum(Nz(~M).^2) ;];
    T_eigen=[nansum(Nxi.^2) nansum(Nxi.*Nyi) nansum(Nxi.*Nzi);
         nansum(Nxi.*Nyi) nansum(Nyi.^2) nansum(Nyi.*Nzi);
        nansum(Nxi.*Nzi) nansum(Nyi.*Nzi) nansum(Nzi.^2) ;];
    eigen = eig(T_eigen); 
    % (λ1/ λ 2) and (λ 1/ λ 3) chapter 4.1.2.1
    DCER1(j,i)=(eigen(1)/eigen(2)); %
    DCER2(j,i)=(eigen(1)/eigen(3));
%     RLOV=sqrt((nansum(Nx(~M))).^2+(nansum(Ny(~M))).^2+(nansum(Nz(~M))).^2);
    RLOV(j,i)=sqrt((nansum(Nxi)).^2+(nansum(Nyi)).^2+(nansum(Nzi)).^2);
    
    DEMt(j,i)=nan; %exclude the central pixel for Roughness and TPI.
    Dsub=DEMt(My,Mx);
    Roughness(j,i)=nanmax(nanmax(abs(Dsub-DEM.z(j,i))));
    TPI(j,i)=DEM.z(j,i)-nanmedian(Dsub(:));
end
end

% out.Roughness_mad=mad(Roughness(:),1);
% out.Roughness_median=median(Roughness(:));
out(4)=mad(Roughness(:),1);
out(3)=nanmedian(Roughness(:));

% The Direction Cosine Of Eigenvalue Ratios DCER
% (equation 2 in Kasai et al., 2009)
% Lower eigenvalue ratios indicate that the unit orientation vector of the cells will have higher degrees of surface roughness

% out.DCER1_median=median(DCER1(:));
% out.DCER1_mad=mad(DCER1(:),1);
% out.DCER2_median=median(DCER2(:));
% out.DCER2_mad=mad(DCER2(:),1);
out(5)=nanmedian(DCER1(:));
out(6)=mad(DCER1(:),1);
out(7)=nanmedian(DCER2(:));
out(8)=mad(DCER2(:),1);

% Hillshade
hills=hillshade(double(DEM.z),DEM.x,DEM.y);
% out.hills_mad=mad(hills(:),1);
out(9)=mad(hills(:),1);

% Resultant Length of Orientation Vectors RLOV; 
%   Greater variations will be displayed for rough topography
% out.RLOV_median=median(RLOV(:));
% out.RLOV_mad=mad(RLOV(:),1);
out(10)=nanmedian(RLOV(:));
out(11)=mad(RLOV(:),1);

%Chang et al., 2019
% TRI, terrain roughness index; ignore this, use Roughness instead.

% TPI, terrain position index; De Reu et al., 2013, Application of the topographic position index to heterogeneous landscapes
% TPI measures the difference between elevation at the central point (z0) and the average elevation (z") around it.
out(12)=nanmedian(TPI(:));
out(13)=mad(TPI(:),1);

% Total curvature,  
% Profile curvature, 
% Plan curvature
%Must smooth DEM before calculating curvature!
if 0 %did this in the main function
DEMs.z = imgaussfilt(DEM.z,2); 
% DEMs.x=DEM.x;DEMs.y=DEM.y;
DEM.z=DEMs.z;
end
[profc,planc,Kt,totalc] = curvature2(DEM); %get NaNs
[profc,planc] = curvature(DEM.z,resr,'both');%avoid NaNs. Pick this one.

% out.profc_median=median(profc(:));
% out.profc_mad=mad(profc(:),1);
% 
% out.planc_median=median(planc(:));
% out.planc_mad=mad(planc(:),1);
out(14)=nanmedian(profc(:));
out(15)=mad(profc(:),1);

out(16)=nanmedian(planc(:));
out(17)=mad(planc(:),1);

if 1 %
% out.Kt_median=median(Kt(:));
% out.Kt_mad=mad(Kt(:),1);
out(18)=nanmedian(Kt(:));
out(19)=mad(Kt(:),1);
end

% out.totalc_median=median(totalc(:));
% out.totalc_mad=mad(totalc(:),1);

out(20)=nanmedian(totalc(:));
out(21)=mad(totalc(:),1);


%plot
if 0
%     set(gcf,'units','points','position',[x0,y0,width,height]);
%     x0=0:300:900; y0=600 290 90
    figure;imagesc(DEM.x*1e-3,DEM.y*1e-3,slopeangle);colorbar; colormap jet;title('slopeangle');view(90,-90)
    set(gcf,'units','points','position',[300,600,300,200])
    figure;imagesc(DEM.x*1e-3,DEM.y*1e-3,theta);colorbar; colormap jet;title('aspect');view(90,-90)
    set(gcf,'units','points','position',[600,600,300,200])

    figure;imagesc(DEM.x(1:stepi:end)*1e-3,DEM.y(1:stepi:end)*1e-3,Roughness(1:stepi:end,1:stepi:end));colorbar; colormap jet;title('Roughness');view(90,-90)
                set(gcf,'units','points','position',[600,90,300,200])
    figure;imagesc(DEM.x(1:stepi:end)*1e-3,DEM.y(1:stepi:end)*1e-3,TPI(1:stepi:end,1:stepi:end));colorbar; colormap jet;title('TPI');view(90,-90)
                set(gcf,'units','points','position',[900,90,300,200])
    figure;imagesc(DEM.x(1:stepi:end)*1e-3,DEM.y(1:stepi:end)*1e-3,DCER1(1:stepi:end,1:stepi:end));colorbar; colormap jet;title('DCER1');view(90,-90)
                set(gcf,'units','points','position',[0,90,300,200])
    figure;imagesc(DEM.x(1:stepi:end)*1e-3,DEM.y(1:stepi:end)*1e-3,DCER2(1:stepi:end,1:stepi:end));colorbar; colormap jet;title('DCER2');view(90,-90)
                set(gcf,'units','points','position',[300,90,300,200])
    figure;imagesc(DEM.x(1:stepi:end)*1e-3,DEM.y(1:stepi:end)*1e-3,RLOV(1:stepi:end,1:stepi:end));colorbar; colormap jet;title('RLOV');view(90,-90)
                set(gcf,'units','points','position',[900,290,300,200])

    hillshade(double(DEM.z),DEM.x,DEM.y,'plotit');view(-90,90);title('hillshade')

    figure;imagesc(DEM.x*1e-3,DEM.y*1e-3,profc);colorbar; colormap jet;title('profile curvature');view(90,-90)
            set(gcf,'units','points','position',[300,290,300,200])
    figure;imagesc(DEM.x*1e-3,DEM.y*1e-3,planc);colorbar; colormap jet;title('plan curvature');view(90,-90)
        set(gcf,'units','points','position',[0,290,300,200])
    figure;imagesc(DEM.x*1e-3,DEM.y*1e-3,Kt);colorbar; colormap jet;title('tangential curvature');view(90,-90)
                set(gcf,'units','points','position',[600,290,300,200])
    figure;imagesc(DEM.x*1e-3,DEM.y*1e-3,totalc);colorbar; colormap jet;title('total curvature');view(90,-90)
        set(gcf,'units','points','position',[900,600,300,200])

end

return
end
