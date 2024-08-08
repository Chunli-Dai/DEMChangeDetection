
function [co]=dhelev_sub(jump,dem,M)
% input: jump: elevation change data;
%        dem: DEM data
%        M: mask matrix of the target area (same grid as jump) to be checked;
% results: dh vs elev plots
    flagplot=0;
    co=[];
    
    dem.z(dem.z==-9999|dem.z==0)=nan;
    %profiles change vs elevation
    elev=dem.z(M);
    if length(elev)<=2;fprintf('\n Not enough DEM data.\n');co=0;return;end
    
    elevchange=jump.z(M);
    
%     elevstep=round(nanmin(elev)):1:round(nanmax(elev));
% %     elevstep=linspace(round(nanmin(dem.z(M))),round(nanmax(dem.z(M))),100);            
%     for k=1:length(elevstep)-1
%         elevdiff_step(k)=nanmedian(elevchange(elev>=elevstep(k)&elev<=elevstep(k+1) ));
%     end
    
    elevdiff_step2=round(nanmin(elevchange)):1:round(nanmax(elevchange));
    elevdiff_step2b=elevdiff_step2(1:end-1); %,elevstep2
    elevstep2=nan*ones(size(elevdiff_step2b)); %nan*ones(length(elevdiff_step2)-1,1);
    for k=1:length(elevdiff_step2)-1
        elevstep2(k)=nanmedian(elev(elevchange>=elevdiff_step2(k)&elevchange<=elevdiff_step2(k+1) ));
    end   
    
    %Criterion 1: Elevations at the negative area are higher than the elevations at the deposits area.
    cri1=nanmedian(elevstep2(elevdiff_step2b<0))>nanmedian(elevstep2(elevdiff_step2b>0));
    %Criterion 2: The volume loss >= volume gain
    cri2=sum(abs(elevchange(elevchange<0)))>sum(elevchange(elevchange>0));
    %Criterion 3: location of positive clusters is at the outside of the negative area;
      %fill the holes of negative clusters;
      BWn=jump.z<0&M;
      BW2 = imfill( BWn ,'holes');
      %see the overlapping area of the above mask with the positive clusters;
      BWp=(jump.z>0&M);
      BWov=BW2&BWp;
      
      % get the ratio of the overlapping area over total positve area;
      % 1 all positive pixels inside the negative cluster (bad); 
      % 0, all pixels outside (good);
      % 0.3, only 30% pixels are inside.
      ratio=sum(BWov(:))/sum(BWp(:));
      
    cri3=(ratio<0.3);
    
    if cri1&&cri2&&cri3 
        co=1;
    else; co=0;
    end
    
    if flagplot==1
%         close (figure(1))
%         figure(1);hold on;plot(elevstep(1:end-1),elevdiff_step,'o','linewidth',3)
%         box on
%         set(gcf,'Color','white')
%         set(gca,'FontSize', 12);
%         set(gcf, 'PaperPosition', [0 0 3 1.5]); 
%         set(gcf, 'PaperSize', [ 3 1.5]); 
%         xlabel('Surface elevation (m)')
%         ylabel('Elevation change (m)')
%         print('-dpng','-r400','dhelevi') 
%         saveas(gcf,'dhelevi','fig')
        
        close (figure(1))
        figure(1);hold on;plot(elevdiff_step2(1:end-1),elevstep2,'o','linewidth',3)
        box on
        set(gcf,'Color','white')
        set(gca,'FontSize', 12);
        set(gcf, 'PaperPosition', [0 0 3 1.5]); 
        set(gcf, 'PaperSize', [ 3 1.5]); 
        ylabel('Surface elevation (m)')
        xlabel('Elevation change (m)')
        print('-dpng','-r400','dhelevi') 
        saveas(gcf,'dhelevi','fig')
    end
    
return
end