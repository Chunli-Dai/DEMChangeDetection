function [idkp,flagcond]=outlier(demp,dempmt,epochorg,model)
%input: model='linear': using a linear trend model for outlier detection;
%       model='constant': using a constant model for outlier detection;
% Modification March  2021: Using Schaffrin's outlier detection

    %make sure the time is sorted.
    [epochorg,idsort]=sort(epochorg);
    demp=demp(idsort);dempmt=dempmt(idsort);

    flagplot=0;
    yr=365.25;
    t=epochorg;
    ni=length(demp);
	flagcond=0; %not used; % if condtion number of AM*P*AM is large, flagcond=1.
    
    Pd=ones(ni,1);Pd(dempmt==0)=0.01; %lesser weight at false match points
    idnn=find(demp~=-9999);idout=[];
    idkp= idnn(~ismember(idnn,idout)); % B(~ismember(B,A)) excluding A from B
    %remove the repeated epochs (dt<10min), and only use the last epoch of the repeated pair. 
    %Important for successful outlier detection (avoid added weight to outliers throught repeats)! 
    dt=diff(epochorg(idkp));Mout=(dt<10/60/24);Mout=[Mout(:);0];idkp=idkp(~Mout);
    
    idkpsv=idkp;

    if length(idkp) <=1; return;end

    %Outlier detection 
    % Change to no iteration, in case of big jump events
    %%%  Linear fitting model: y=a + b(t-tm) 
    
    algorithm=2; % 1, y=a+bt; 2, y=a; %Constant works better.
      %/home/dai.56/chunliwork/ice/Barnes1c/bp0/tp1_constant.fig vs tp1_trend.fig
    option=2; %default is 2; %1 linear or constant ; 2 constant +jump; The default is 2;
    if exist('model','var')
        if strcmp(model,'linear')
            algorithm=1;
            option=1;
        elseif strcmp(model,'constant')
            algorithm=2;
            option=2;
        end
    end
    

    epoch=epochorg/yr;

    %Option 1: use contant or liner trend to fit the time series
    %Option 2: first do difference; then to use a constant 0 to fit the
    %difference time series; (This is for a case of known large change)
    
    % detecting outliers
    if option ==1 % go for option 1

    yobs=demp(idkp);
    Pd_obs=Pd(idkp);
    tobs=epoch(idkp);
    [Mout]=outlierSchaffrin(tobs,yobs,Pd_obs,algorithm);
%     Mout=isoutlier(yobs); %fails when you have 10 measurements in 2010, and 1 measurement (with real change) in 2020.

    idout=idkp(Mout);

    idkp= idkpsv(~ismember(idkpsv,idout));
    
    else %option 2; refers to timeseries2bp2.m
        
    yobs=demp(idkp);
    td=t(idkpsv(1:end-1));
    yd=yobs(2:end)-yobs(1:(end-1));
    
    idkp=[1:length(yd)]';
%     [Mout]=outlierSchaffrin(td,yd,[],2); %% isoutlier works better for multiple outliers
    Mout=isoutlier(yd);
    idout=idkp(Mout); %idout's value is the index of idkp=[1:length(yd)]';
    idkp= idkp(~ismember(idkp,idout));

    %match the outlier in the difference time series to the original y time sereis
    if 0 %bad algorithm.
    idout=idkpsv(idout+1);%half of them are good ones, one of them might be the actual jump.
    else
    idout2=unique([idout(:);idout(:)+1]); %get the id of two points for all the detected outliers in the difference time series (the difference is large)
    idkp2=unique([idkp(:);idkp(:)+1]); %get the id of two points if the difference is small
    idout=idout2(~ismember(idout2,idkp2)); % exclude pairs with large difference, but keep pairs with small difference.
    %It could be the fisrt obs or the last obs of the time series that was bad.
    idd=[];
    for j=1:length(idout)
        if idout(j)==1 || idout(j)==length(idkpsv); idd=[idd;j];end
    end
    idout(idd)=[];
    
    idout=idkpsv(idout);
    end

    % give extreme low weight to the outliers.
    idkp= idkpsv(~ismember(idkpsv,idout));
     if flagplot==1
        figure
        set(gcf,'Color','white')
        set(gca,'FontSize', 18);
        set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
        hold all
        plot(t(idkp)/yr,demp(idkp),'b>:','MarkerSize',12,'linewidth',4)
        plot(t(idout)/yr,demp(idout),'ks','MarkerSize',12,'linewidth',4)
        title('Original time series and outliers')
        box on
        %     plot(t(idkp(idsort)),fit(idsort),'g*-','MarkerSize',12,'linewidth',4)
    end
        
    end

	close all
end
