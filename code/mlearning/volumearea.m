function [alpha,beta,Mflag]=volumearea(vol,surfacearea,volstd,flagmodel)
% generate the volume-area relationship (Jayboyedoff et al., 2020)
% vol, volume n by 1 vector, m^3;
% surfacearea, area n by 1 vector in m^2;
% Output: alpha, beta, coefficients of the model;
%           Mflag, 1 input data is within the model range, 0 input data is outside
%           of the model range.

    logvol=log(vol(:));
    logarea=log(surfacearea(:));
    Mflag=false(length(vol),1);
    alpha=nan;beta=nan;
    
    % logV=α+βlogA 
    
    %Initializing
    P=[];lenf=length(logvol);mp=2;
    stdi=1;flagstd=0; %default
    if exist('volstd', 'var')
        if ~isempty(volstd)
            flagstd=1;
        fprintf('\n Use given volume uncertainty. \n');
        ystd=log(volstd(:)); %m^3
        else
            flagstd=0;
        end
    end
    
    if flagstd==0
        fprintf('\n Use estimated sigma02hat for volume uncertainty. \n');
        ystd=[stdi*ones(size(logvol))];
    end

    if lenf <=2e3
        P=diag(ystd(:).^-2); 
    else %equal weight
        fprintf(['\n volumearea.m: input vector is too long, use equal weight, len=',num2str(lenf),'.\n'])
    end
    AM=zeros(lenf,mp);  %
    AM(:,1)=1.;
    AM(:,2)=logarea;

    yobs=logvol; xv=logarea;

    % Least squares; 
    % flagmodel=1; %1 use given model Jayboyedoff et al. (2020), 2, fit the model using input data.

    if flagmodel==2 && lenf> mp 
        fprintf(['\n Apply Least-squares Estimation. \n'])
        if lenf <=2e3
            cdA=cond(AM'*P*AM);
        else
            cdA=cond(AM'*AM);
        end
        if(cdA>1e5); display(['condtion number of AM*P*AM is large: ',num2str(cdA)]); return;end
        
        if lenf <=2e3
                var=inv(AM'*P*AM);
                est=var*AM'*P*yobs;
        else
                var=inv(AM'*AM);
                est=AM\yobs; % equal weight
        end
        etilde=yobs-AM*est;
    
        %reestimate the reference variance, posteriori;
        if lenf <=mp %lenf <=2 
            sigma02hat=NaN;%oflag=2;
        else
            if lenf <=2e3
                sigma02hat=etilde'*P*etilde/(lenf-mp); %mp=2
            else
                sigma02hat=etilde'*etilde/(lenf-mp); %mp=2
            end
        end
        stdi=sqrt(sigma02hat);
        var=var*sigma02hat;
    
        % reestimate std
        if exist('volstd', 'var')
%             fprintf('\n Use given volume uncertainty. \n');
        else        
            ystd=ystd*sqrt(sigma02hat);
        end
    
        for m=1:mp 
            eststd(m:m)=sqrt(var(m,m)); 
        end
        fprintf('\n volumearea.m alpha beta: %f +/- %f ,  %f +/- %f .\n',est(1), eststd(1),est(2),eststd(2))

%         xvfit=xv;
        xvfit=7:0.1:17; %log([2000 500e3*50])
        xvfit=xvfit(:);
        AMfit=zeros(length(xvfit),mp);  %
        AMfit(:,1)=1.;
        AMfit(:,2)=xvfit;
        
        fit=AMfit*est;
        fitstdall=AMfit*var*AMfit'; %Matrix can be too big.
        
        rate=est(2);ratestd=sqrt(var(2,2));
    
    else
        fprintf(['\n Use the values in Jayboyedoff et al., 2020. \n'])
%         The coefficient (β) for deep-seated bedrock landslides ranges from 1.25 to 1.5 (Jayboyedoff et al., 2020), 
%         and α ranges from -2.996 to -0.26. 
        est=[(-2.996-0.26)/2; (1.25 + 1.5)/2 ];
        eststd=[(2.996-0.26)/2; (-1.25 + 1.5)/2 ];

        % AM
        % 
%         xvfit=9:0.1:13;
        xvfit=7:0.1:17; %log([2000 500e3*50])
        xvfit=xvfit(:);
        AMfit=zeros(length(xvfit),mp);  %
        AMfit(:,1)=1.;
        AMfit(:,2)=xvfit;

        fit=AMfit*est;
        var=zeros(mp,mp);
        for m=1:mp 
            var(m,m)=eststd(m:m).^2; 
        end
        fitstdall=AMfit*var*AMfit';

        rate=est(2);ratestd=sqrt(var(2,2));
        fprintf('\n volumearea.m alpha beta (Jayboyedoff et al., 2020): %f +/- %f ,  %f +/- %f .\n',est(1), eststd(1),est(2),eststd(2))

    end
    alpha=est(1);beta=est(2);

    %% Get Mflag
    etilde=yobs-AM*est;
    if 0 %large memory use
        ymodelstd_all=AM*var*AM'; %Matrix can be too big.
        ymodelstd=sqrt(diag(ymodelstd_all));
    else %small memory use
        ymodelstd=zeros(lenf,1);
        for i=1:lenf
            ymodelstd(i)=sqrt(AM(i,:)*var*AM(i,:)');
        end
    end

    Mflag=abs(etilde)<=ymodelstd;

    %% plotting
    formatSpec = '%6.2f';

    fitstd=zeros(length(fitstdall(:,1)),1);
    for j=1:length(fitstdall(:,1))
        fitstd(j)=sqrt(fitstdall(j,j));
    end
    
    figure %(2)  
    set(gcf,'Color','white')
    set(gca,'FontSize', 12);
    set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
    hold all
    plot(xv(Mflag),yobs(Mflag),'b.','linewidth',2,'markersize',8) %observation used for fitting
    plot(xv(~Mflag),yobs(~Mflag),'r.','linewidth',2,'markersize',8) %observation used for fitting
    plot(xvfit,fit,'r-','linewidth',2,'markersize',14) %linear fit
    hold on
    if 0
        % xv needs to be ordered
        [xvs,ids]=sort(xv);
        shadedErrorBar(xvs,yobs(ids),ystd(ids),'b',1)
        [xvfits,ids]=sort(xvfit);
        shadedErrorBar(xvfits,fit(ids),fitstd(ids),'r',1)
    end

    text(max(xv)-(max(xv)-min(xv))*0.3,max(yobs)-(max(yobs)-min(yobs))*0.4 ,['Trend=',num2str(rate,formatSpec),'\pm',num2str(ratestd,formatSpec),''],'FontSize',12)
    legend('Observations','Linear Fit')
    xlabel(['Logarithm of Area (m^2)'],'FontSize',12)
    ylabel(['Logarithm of Volume (m^3)'],'FontSize',12)
    box on

    ofile=['volarea.fig'];
    saveas(gcf,ofile,'fig')


return
end




