function [Mouts, Mout]=matchclusters(Mnc0,Mpc0,flagcase)
% Match the negative (100 m buffer) and positive clusters.
% Mnc0: Input negative clusters. This can be a small area with a single
%       cluster (make sure to make buffer zone (> buff1km) around the cluster)
% Mpc0: Input positive clusters. This can be a large area with all clusters.
% 
% Mout: output matrix of all matched pairs of negative and positive clusters, same size as Mnc0.
        %connect the matched pair.
% Mouts: each pair of negative and positive clusters stored in structure.
% flagcase=0 : default; if no match, still output the negative cluster.
% flagcase=1 : Find the clusters in Mnc0 that include Mpc0;if no match, ignore/remove the negative cluster. no buffer zone, buff=1;
if nargin==2; flagcase=0;end

buff=20; %max distance to find the closest positive clusters. e.g. 100m or 20 m
buff1km=1e3/2/2; %max distance of the extent of the positive clusters.
if flagcase==1;buff=1;end
% Mouts=[];
Mout=[];

resx=mean(Mpc0.x(2:end)-Mpc0.x(1:end-1));resy=mean(Mpc0.y(2:end)-Mpc0.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);

[n1,m1]=size(Mnc0.z);[n2,m2]=size(Mpc0.z);
if n1==n2&&m1==m2
	flag=1 ; %two input matrix are the same size; They are all big matrix.
else
	flag=2; %two input matrix are different sizes; positive clusters larger than the negative clusters
    %make the grids same;
    Mpc0_org=Mpc0;
    Mpc1.z=interp2(Mpc0.x,Mpc0.y,int8(Mpc0.z),Mnc0.x,Mnc0.y','*nearest',0);
    Mpc0=Mnc0;Mpc0.z=logical(Mpc1.z); clear Mpc1
end

CC = bwconncomp(Mnc0.z);

[X,Y]=meshgrid(Mpc0.x,Mpc0.y);
Mout=Mnc0; Mout.z=false(size(Mnc0.z));
Mouts(CC.NumObjects)=struct('x',[],'y',[],'z',[]);

for ki=1:CC.NumObjects
    %get the individual cluster.
    BW3=false(size(Mnc0.z));
    BW3(CC.PixelIdxList{ki})=Mnc0.z(CC.PixelIdxList{ki});
    
%     if flag==2
%     %match the size to Mpc0
%     BW3r=interp2(Mnc0.x,Mnc0.y,int8(BW3),Mpc0.x,Mpc0.y','*nearest',0);
%     BW3=logical(BW3r);
%     end
    
    xmin=min(X(BW3));xmax=max(X(BW3));ymin=min(Y(BW3));ymax=max(Y(BW3));
    rang2=[xmin-buff1km xmax+buff1km ymin-buff1km ymax+buff1km];
    Mx=Mpc0.x>=rang2(1)&Mpc0.x<=rang2(2);My=Mpc0.y>=rang2(3)&Mpc0.y<=rang2(4);
    
    %cropped positive clusters
    Mpc.x=Mpc0.x(Mx);Mpc.y=Mpc0.y(My);
    Mpc.z=Mpc0.z(My,Mx);
    
    Mnc=Mpc;Mnc.z=BW3(My,Mx);
    
    %1\ Dilate around the negative clusters. 
    %2\ Dilate the positive clusters;
    %3\ Find the closest cluster (the maximum overlapping area after buffering)
    if flagcase==1 %if no match, ignore/remove the negative cluster. no dilation;
        Mnc_dil=Mnc.z;
    elseif flagcase==0 %default; there is dilation; if buff =0, dilated map is all false.
        Mnc_dil=(imdilate(Mnc.z,ones(round(buff*2/resr))));
    end
    CCp = bwconncomp(Mpc.z);
    overlap=zeros(CCp.NumObjects,1);
    clear Moverlap;
    for jp=1:CCp.NumObjects
        BW3p=false(size(Mpc.z));
        BW3p(CCp.PixelIdxList{jp})=Mpc.z(CCp.PixelIdxList{jp});
        BW3p_dil=(imdilate(BW3p,ones(round(buff*2/resr))));
       
        overlap(jp)=sum(sum(Mnc_dil&BW3p_dil));
        Moverlap{jp}=(imdilate(Mnc_dil&BW3p_dil,ones(round(buff*2/resr))));
    end

    [overlap_max,id1]=max(overlap);
    flagbad=0; %1 bad; 0 good
    if isempty(overlap_max)
        flagbad=1;
    elseif overlap_max<2 %%(less than 2 overlapping pixels)
        flagbad=1;
    end
    
    %collect the matched positive clusters (the largest one)
    BW3p=false(size(Mpc.z));
    if ~(CCp.NumObjects==0||flagbad==1) %There is a match .
        jp=id1;BW3p(CCp.PixelIdxList{jp})=Mpc.z(CCp.PixelIdxList{jp});
        BW3p=BW3p|Moverlap{jp}; %add the connecting pixels Moverlap{jp} to connect two clusters.
    else 
        fprintf('\n No positive cluster matched for the given negative cluster.\n')
        if flagcase==1 %if no match, ignore/remove the negative cluster.
            continue; %skip the rest
        end
    end
    Mpc1=Mpc;Mpc1.z=BW3p;
    
    %the common area
    Mnp=Mpc; Mnp.z=Mnc.z|Mpc1.z;
    
    %crop out void area to keep the matrix small.
    [X1,Y1]=meshgrid(Mnp.x,Mnp.y);
    xmin=min(X1(Mnp.z));xmax=max(X1(Mnp.z));ymin=min(Y1(Mnp.z));ymax=max(Y1(Mnp.z));
    rang2=[xmin-10  xmax+10 ymin-10 ymax+10]; %10 m buffer 
    Mxs=Mnp.x>=rang2(1)&Mnp.x<=rang2(2);Mys=Mnp.y>=rang2(3)&Mnp.y<=rang2(4);
    %cropped clusters
    Mnps.x=Mnp.x(Mxs);Mnps.y=Mnp.y(Mys);
    Mnps.z=Mnp.z(Mys,Mxs);   
    
    Mouts(ki)=Mnps;%small

    %make the output same size as Mout/Mpc0;
%     Mnpi=interp2(Mnp.x,Mnp.y,int8(Mnp.z),Mout.x,Mout.y','*nearest',0); %slow
     Mnpi=false(size(Mout.z));Mnpi(My,Mx)=Mnp.z; 

    Mout.z=Mout.z|logical(Mnpi);
        
end %ki

%dilate the clusters to make the matched pair connected, so that output is only one cluster.
% Mout.z=(imdilate(Mout.z,ones(round(buff*2/resr))));  %Too much dilation!

return
end
