function [oflag,trend,trest,rate,ratestd,eqm,eqs,eqe]=timeseries2(demp,dempmt,epochorg,timefix,eqepoch,flagplotsv,jq)
% v1: constant fitting of the time series
% v2: fitting of the difference of the time series. 
%   Output: 
%   oflag (right to left): bit 1, bit2, bit3, bit4, bit5, bit6
%	bit 1: 1, length(data) <=1; 0 otherwise;
%       bit 2: 1, length(data) ==2, 0 otherwise;
%       bit 3: 1, length(data) >=3, but no jump;
%       bit 4: 1, no measurement before or after event time;
%       bit 5: 1, there is a confident jump
%       bit 6: 1, the change is at the start or end of the time series.
%     Output cases flag (oflag): 1) no data, or length of data <=1; eqm=0;
%               2) only two data dt>0; 3) length of data >=3, no jump; eqm=0;; 
%               4)If no obs before and after event eqm=0;;5) there is a jump ; 
%               6) there is a jump, but jump is at the start or end.
%          trend, trest: the change of elevation and its uncertainty (meter).
%          rate,ratestd: the change rate of elevation and uncertianty (m/yr). 
%          eqm, eqs,eqe: the middle epoch, start epoch and end epoch of
%                        detected event (the event in year after AD 0).
%	  For trend estimation (ice melting): eqm, eqs, eqe is the average, start, end epoch of all used data.

%   Revision March 2019: Adapt for volcano eruption (a duration of time); exclude event jump from outlier detection.

%   Revision November 2019: revisit the outlier detection, to use the
%   original time series instead of its difference. The difference makes it
%   noiser, and one outlier is turned to two outliers.
%

      save outlierex1.mat demp dempmt epochorg timefix eqepoch flagplotsv jq
      if flagplotsv ==1
	ofile=['outliertp',num2str(jq),'.mat'];
        save(ofile,'demp','dempmt', 'epochorg','timefix','eqepoch','flagplotsv');%let jq=0
      end
      
      constant
      oflagbit=zeros(6,1);
      oflag=0;

    % test 
    if 0 %test run; initial parameters.
        lateq=60.172038; loneq= -141.172456; % Oct. 17 2015 , largest landslide
        [xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
        teq=datenum('20151017','yyyymmdd');
        texteq={'20151017Landslide'};
        eqepoch=teq;timefix=1; % 1 assume event time known; 0 assume time unknown
        eqepoch=0;timefix=0;
        flagplotsv=0;
    end
%     tic
    
    %make sure the time is sorted.
    demporg=demp;dempmtorg=dempmt;epochorgorg=epochorg;    
    [epochorg,idsort]=sort(epochorg);
    demp=demp(idsort);dempmt=dempmt(idsort);
    
    %remove the data inbetween eruptions ->1\ easier outlier detection;2\easier time series
    % demp,dempmt,epochorg
    if flagvolc==1
        idd=find(epochorg>eqepoch(1)&epochorg<eqepoch(2));
%         fprintf(['\n Files within the duration of eruption (to be deleted):','','\n'])
%         datestr([epochorg(idd)])
        epochorg(idd)=[];demp(idd)=[];dempmt(idd)=[];
    end
%     flagvolc=1; %1 volcano event, a duration of time;

    flagplot=0;
    yr=365.25;
    t=epochorg;
    maxcount=0;
    ni=length(demp);
    AMa=ones(ni,1);
    epoch=epochorg/yr;
    tm=mean(epoch);
    epochs=0; %datenum('2000/1/1')/yr;
    formatSpec = '%6.1f';
    neq=1;
    multi=3;
    tryall=0; %  1 try all epochs for jump detection; 0 faster, try detected range; 
    %Wrong pick based on sigma0hat (actual jump treated as outlier); use maximum jump as best event criteria, it works!
    option=1; % 1 faster, use the wider time interval for the event time, might lose the closest jump epoch;
               % 2 slower, but more refinement on the time interval and more
               % measurements included. 
    
    Pd=ones(ni,1);Pd(dempmt==0)=0.01; %lesser weight at false match points
    idnn=find(demp~=-9999);idout=[];idoutpre=[];
    idkp= idnn(~ismember(idnn,idout)); % B(~ismember(B,A)) excluding A from B
    idkpsv=idkp;

    if isempty(idkp)
        dt=0; trend=0;trest=0; eqm=0;eqs='0000/00/00';eqe='0000/00/00';oflag=1;rate=0;ratestd=0;
        oflagbit(1)=1;oflag=bit2value(oflagbit);
        return;
    elseif length(idkp) <=2 && length(idkp) >=1
    %write the trend, date, and flag=2
        eqs=datestr(epoch(idkp(1))*yr,26);
        eqe=datestr(epoch(idkp(end))*yr,26);
        dt=epoch(idkp(end))-epoch(idkp(1));
        eqm=(epoch(idkp(end))+epoch(idkp(1)))./2-epochs; % the event in year after AD 0 (2000/1/1/)
        if dt==0; 
            trend=0;trest=0;eqm=0;oflag=1;rate=0;ratestd=0;
            oflagbit(1)=1;oflag=bit2value(oflagbit);
	    return
        elseif timefix==0  %: in case of timefix=1 and no data befor/after event, there should be no output. 
            %oflag=2;
            oflagbit(2)=1;oflag=bit2value(oflagbit);
            if algorithmin==1 % 1 linear (ice melting); 2 constant (landslides) ; 3 constant + linear (Okmok volcano)
                %trend=0;trest=0;rate=(demp(idkp(end))-demp(idkp(1)))/dt;ratestd=NaN;
                trend=0;trest=0;rate=NaN;ratestd=NaN;
            elseif algorithmin==2
                trend=(demp(idkp(end))-demp(idkp(1)));
                trest=NaN; % no posteriori std.
                rate=0;ratestd=0;
            elseif algorithmin==3
                trend=(demp(idkp(end))-demp(idkp(1)));
                trest=NaN; % no posteriori std.
                rate=NaN;ratestd=NaN;
            end
	    return
        end
	
    end

    %Outlier detection 
    if algorithmin==2||algorithmin==3
        model='constant';
    elseif algorithmin==1
        model='linear';
    end
    [idkp,flagcond]=outlier(demp,dempmt,epochorg,model);
    idout=idkpsv(~ismember(idkpsv,idkp)); %not including demp==-9999;
    
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
    
    %% Detection of event time, assuming only one event.
    % Restore the jump of the event
    jumpflag=1;
    if timefix == 1 ||algorithmin==1 %timefix =1 event time known; 0 event time not known.
         if algorithmin==1 %ice melting
             neq=0;%oflag=3;
             oflagbit(3)=1;oflag=bit2value(oflagbit);
             % to be compatible
             no2=1;eqepoch2=eqepoch; 
    %        eqm=eqepoch/yr-epochs;eqs=datestr(eqepoch,26);eqe=eqs;
    %        eqm=eqepoch/yr-epochs;eqs=datestr(eqepoch,26);eqe=eqs;
             eqsi=min(epoch(idkp));eqei=max(epoch(idkp));%begin and end time of good measurements; 
             eqs=datestr(eqsi*yr,26);eqe=datestr(eqei*yr,26);
             dt=eqei-eqsi; % in years
%            eqm=(eqei+eqsi)./2-epochs; %in years
             eqm=mean(epoch(idkp))-epochs;; %in years %eqm will be replaced later.
             if ~isempty(eqm)
                 eqm2(1)=eqm;eqs2{1}=eqs;eqe2{1}=eqe;
             else
                 eqm2(1)=0;eqs2{1}='0000/00/00';eqe2{1}='0000/00/00';
             end
         else
         neq=length(eqepoch);
         eqm=eqepoch/yr-epochs;
%          algorithm=2; %constant
         eqs=datestr(eqepoch,26);eqe=eqs;
         jumpflag=1;
%        no2=1;eqepoch2(1)=eqepoch;eqm2(1)=eqm;eqs2{1}=eqs;eqe2{1}=eqe;
	
         if flagvolc==1 %one period of volcanic eruption. eqm2 eqs2 eqe2 are choices of output date.
             neq=1; teq=eqepoch(1);teq2=eqepoch(2); eqepoch2=eqepoch;
             no2=1;eqm2=eqm(2);eqs2{1}=datestr(teq,26);eqe2{1}=datestr(teq2,26);
         elseif neq==1; teq=eqepoch;teq2=eqepoch;
             no2=1;eqepoch2(1)=eqepoch;eqm2(1)=eqm;eqs2{1}=eqs;eqe2{1}=eqe;
         else
             no2=1; %number of attempts of searching.
             eqepoch2=eqepoch;
             fprintf(['\n Timeseries2.m warning: multiple events, not developed! \n'])
         end

          % use outliers when no enough obs before or after event.
         if sum(t(idkp)>teq2) && sum(t(idkp)<teq) 
             %do nothing
         elseif sum(t(idkp)>teq2)==0 && sum(t(idkp)<teq)% keep outliers after the event as obs.
             idkp=[idkp;idout(t(idout)>teq2)];
             idout=idout(t(idout)<=teq);
         elseif sum(t(idkp)<teq)==0 && sum(t(idkp)>teq2) % keep outliers before the event as obs.
             idkp=[idkp;idout(t(idout)<teq)];
             idout=idout(t(idout)>=teq2);
         else  % use all outliers
             idkp=[idkp;idout];idout=[];
         end    
         end %if algorithmin

    elseif timefix ==0 % need to detect the event epoch
        %July 2020, in search of jump time exclude bad matchtag points;
        if 1 %exclude bad matchtag points in search of jump time
        idkpsv1=idkp;
        idmt=find(dempmt==0);
        idkp=idkp(~ismember(idkp,idmt));
        end
        
        T6=demp(idkp); % idkp must contain 1 element.
        yd=T6(2:end)-T6(1:(end-1));

        id2=[];
        if tryall==1 % try all epochs for jump detection
            jumpflag=1;
            id2=sort([idkp;idout]); % try all valid data (not -9999) intervals.
            no2=(length(id2)-1);%number of attempts to find the jump.
            eqepoch2=zeros(no2,1);eqm2=zeros(no2,1);eqs2=cell(no2,1);eqe2=cell(no2,1);
            for io2=1:no2
                tpre=epoch(id2(io2));tpost=epoch(id2(io2+1));
                eqepoch2(io2)=(tpre+tpost)./2*yr;
                eqm2(io2)=eqepoch2(io2)/yr-epochs;
                eqs2{io2}=datestr(tpre*yr,26);
                eqe2{io2}=datestr(tpost*yr,26);
            end           
        elseif length(T6)<=1  % or try all
           jumpflag=0; 
        elseif length(T6)>=2 %lenf >=3 % % , flag=3.  % to do short length
           %detect jump using two constant fitting model
           % if jump exists, no trend; if no jump, use linear fitting model
           [ydmax,id]=max(abs(yd));
           idall=1:length(yd);
           idall= idall(~ismember(idall,id));
           ydstd=std(yd(idall)); %
           if (abs(yd(id)-mean(yd(idall)))>3*ydstd); jumpflag=1;else jumpflag=0;end
           if jumpflagc==1;jumpflag=1;end %always estimate the change.
           if jumpflag==1 
%                algorithm=2; %constant
               eqepoch=(epoch(idkp(id))+epoch(idkp(id+1)))./2*yr;
               eqm=eqepoch/yr-epochs;
               eqs=datestr(epoch(idkp(id))*yr,26);
               eqe=datestr(epoch(idkp(id+1))*yr,26);
               no2=1;eqepoch2(1)=eqepoch;eqm2(1)=eqm;eqs2{1}=eqs;eqe2{1}=eqe;
               
               % 1 faster, use the wider time interval for the event time, might lose the closest jump epoch;
               % 2 slower, but more refinement on the time interval and more
               % measurements included.
               if option ==2 % if there are outliers within the above range (from only good points), try each interval at outliers within the range.
                  tr1=epoch(idkp(id))*yr;tr2=epoch(idkp(id+1))*yr;
                  id2=sort([idkp(id);idkp(id+1);idout(t(idout)>tr1&t(idout)<tr2)]); % all within the range.
                  no2=(length(id2)-1);
                  eqepoch2=zeros(no2,1);eqm2=zeros(no2,1);eqs2=cell(no2,1);eqe2=cell(no2,1);
                  for io2=1:no2
                        tpre=epoch(id2(io2));tpost=epoch(id2(io2+1));
                        eqepoch2(io2)=(tpre+tpost)./2*yr;
                        eqm2(io2)=eqepoch2(io2)/yr-epochs;
                        eqs2{io2}=datestr(tpre*yr,26);
                        eqe2{io2}=datestr(tpost*yr,26);
                  end  
               end
%            else;algorithm=1; %linear
           end
        end
        
        if 0 %lenf >=4 && 0 % flag=4 % instablity when data points are less, trend
        clear AM
            AM=epoch(idkp(idsort(2:end)))-epoch(idkp(idsort(1:end-1)));
            est=AM\yd;
            etilde=yd-AM*est;
            sigma02hat=etilde'*etilde/(length(yd)-1);
            [emax,id]=max(abs(etilde));     
            if(abs(emax)>=3*sqrt(sigma02hat));jumpflag=1;else jumpflag=0;end
            if jumpflag==1
               eqs=datestr(epoch(idkp(idsort(id)))*yr,26);
               eqe=datestr(epoch(idkp(idsort(id+1)))*yr,26);
               eqepoch=(epoch(idkp(idsort(id)))+epoch(idkp(idsort(id+1))))./2*yr;
            end
            algorithm=1; %linear
        end
        
        idkp=idkpsv1;% restore old kept data which include bad matchtag points 
        
    end % timefix or not
    % End of Detection of the event time

    
          %% Jump fitting
 %%%  Linear fitting model: y=a + b(t-tm) 
%      +d11*(0 or 1);     
        if jumpflag ==0 % no detected jump
            dt=0; trend=0;trest=0; eqm=0;eqs='0000/00/00';eqe='0000/00/00';oflag=3;rate=0;ratestd=0;
            oflagbit(3)=1;oflag=bit2value(oflagbit);
            return
        end
        
%         mpnoj=1; % % number of parameters except the jump parameter.
        algorithm = algorithmin;
        if algorithmin==3||algorithmin==1
            mpnoj=2; %no jump
%             algorithm = 3; % 1 linear; 2 constant ; 3 constant + linear
        elseif algorithmin==2
            mpnoj=1;
        end
        
        % loop for trying different jump epochs.
        max2=0;io2m=1;trend=0;trest=0;rate=0;ratestd=0;
        for io2=1:no2
            eqepoch=eqepoch2(io2);%eqm=eqm2(io2);eqs=eqs2(io2);eqe=eqe2(io2);
            if flagvolc==1;eqepoch=eqepoch2(2);	 end   

            %Initializing
            P=[];
            % In case of obs are mistakenly taken as outliers, keep these
            % outliers but with lesser weight.
            epoch=epochorg([idkp;idout])/yr; 
            T6=demp([idkp;idout]);
            stdi=1; stdio=100; % %lesser weight for the detected outliers.
            T6std=[stdi*ones(size(idkp));stdio*ones(size(idout))];
	    
            if 1 %do not apply; in case of too many points given less weight, only 3 point with good weight, 
             %-> super small sigma0hat, excluding all points with less weight.
             % apply; if not, too many bad results.
            id=find(dempmt(idkp)==0); %if matchtag ==0, weight is less also.
            T6std(id)=stdio;
            end

            [epoch,idsort]=sort(epoch);T6=T6(idsort);T6std=T6std(idsort);
            
            if flagplot==1
                figure
                set(gcf,'Color','white')
                set(gca,'FontSize', 18);
                set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
                hold all
                plot(epoch,T6,'b>:','MarkerSize',12,'linewidth',4)
                errorbar(epoch,T6,T6std,'b')
                title(['Measurements and its uncertainty for jump estimation; try event epoch i=',num2str(io2)])
            end
            
%             T6std=stdi*Pd(idkp(idsort)).^-0.5;

            
            idout2=[];
        for iter=1:2 %1:2 % second iteration for outlier detection using the jump fitting.
 
            epoch(idout2)=[];T6(idout2)=[];T6std(idout2)=[]; lenf=length(T6);
            if iter==2;T6std=ones(size(T6));end % offer equal weights for all obs included.
            P=diag(T6std(:).^-2); 
            
            mp=mpnoj+neq;
            %construct ideq
            ideq=[];
            for i=1:neq
                id=find(epoch*yr >= eqepoch);
                if not(isempty(id)) && length(id)<lenf
                    ideq{i}=id;
                end    
            end
            %flag the types 
            if isempty(ideq) || lenf<mp %lenf<=mp %|| flagcond ==1 % If no obs before and after event.
%                fprintf(['\n No measurements before or after the event!\n'])
                if algorithmin==2 ||algorithmin==3
                    fprintf(['\n No measurements before or after the event!\n'])
                    dt=0; trend=0;trest=0; eqm=0;eqs='0000/00/00';eqe='0000/00/00';oflag=4;rate=0;ratestd=0;
        	    oflagbit(4)=1;oflag=bit2value(oflagbit);
                    continue;
                end
            elseif length(ideq{1})==1||length(ideq{1})==(lenf-1)
                %oflag=6;
                oflagbit(6)=1;oflag=bit2value(oflagbit);
            else
                %oflag=5;
            end

            AM=[];
            AM=zeros(lenf,mp);  %
            AM(:,1)=1.;
            if (algorithm == 1 ) %linear trend            
                AM(:,2)=epoch-tm; %unit year
            elseif (algorithm == 2) % constant
                %do nothing
            elseif (algorithm == 3) % constant + trend
                AM(:,2)=epoch-eqepoch/yr;
                Mid=false(lenf,1);Mid(ideq{1})=1;
                AM(~Mid,2)=0;
            end
            for j=1:neq
            AM(:,mpnoj+j)=0.; AM(ideq{j},mpnoj+j)=1.;
            end            
            
            yobs=T6;

            cdA=cond(AM'*P*AM);
            if(cdA>1e5); display(['condtion number of AM*P*AM is large: ',num2str(cdA),';jq:',num2str(jq)]); continue;end

            var=inv(AM'*P*AM);
            est=var*AM'*P*yobs;
%             est=AM\yobs; % equal weight
            
            etilde=yobs-AM*est;
            %reestimate the reference variance, posteriori;
            if lenf ==mp %lenf <=2 
                sigma02hat=NaN;%oflag=2;
                oflagbit(2)=1;oflag=bit2value(oflagbit);
            else
                sigma02hat=etilde'*P*etilde/(lenf-mp); %mp=2
            end
            stdi=sqrt(sigma02hat);
            idout2=find(abs(etilde)>=multi*stdi);
            
            var=var*sigma02hat;
            T6std=T6std*sqrt(sigma02hat);
            
            for m=1:mp 
%               eststd(m:m)=sqrt(var(m,m)); 
            end

            fit=AM*est;
            
            %1 linear (ice melting); 2 constant (landslides) ; 3 constant + linear (Okmok volcano)
            if algorithm == 2
                trend=est(mpnoj+1);trest=sqrt(var(mpnoj+1,mpnoj+1));
                rate=0;ratestd=0;
            elseif algorithm == 3
                trend=est(mpnoj+1);trest=sqrt(var(mpnoj+1,mpnoj+1));
                rate=est(mpnoj);ratestd=sqrt(var(mpnoj,mpnoj));
            elseif algorithm == 1
                rate=est(mpnoj);ratestd=sqrt(var(mpnoj,mpnoj));
                trend=0;trest=0;
            end

            jumpflag=1;
            if trest>=abs(trend)
                jumpflag=0; % considering rerun the fitting
            end
            if abs(trend)>3*trest & lenf>=4
                %oflag=5;
                oflagbit(5)=1;oflag=bit2value(oflagbit);
	    end
                
    % %%% plotting
        if flagplotsv ==1
	ofile=['outliertp',num2str(jq),'.mat'];
%        save(ofile,'demp','dempmt', 'epochorg')

        % % For plotting figures       
        epochp=epoch*yr; % change unit to day
        if eqepoch>=min(epochp) && eqepoch<=max(epochp)
            epochfit=sort([epochp;eqepoch(:)-0.1;eqepoch(:)+0.1;])/yr; %in year
        else
            epochfit=sort([epochp;])/yr; %in year
        end
        if flagvolc==1;epochfit=sort([epochp;eqepoch2(:);])/yr;end %in year
        ideq=[];
        for i=1:neq
            id=find(epochfit*yr >= eqepoch);
            if not(isempty(id)) && length(id)<length(epochfit)
                ideq{i}=id;
            end    
        end
        AM=[];
        AM=zeros(length(epochfit),mp);  %
        AM(:,1)=1.;
        if (algorithm == 1 ) %linear trend            
            AM(:,2)=epochfit-tm;
        elseif (algorithm == 2) % constant
            %do nothing
        elseif (algorithm == 3) % constant + trend
                AM(:,2)=epochfit-eqepoch/yr;
                Mid=false(lenf,1);Mid(ideq{1})=1;
                AM(~Mid,2)=0;
        end
        for j=1:neq
        AM(:,mpnoj+j)=0.; AM(ideq{j},mpnoj+j)=1.;
        end  
        fit=AM*est;
        fitstdall=AM*var*AM';
        fitstd=zeros(length(fitstdall(:,1)),1);
        for j=1:length(fitstdall(:,1))
            fitstd(j)=sqrt(fitstdall(j,j));
        end
        figure  
        set(gcf,'Color','white')
        set(gca,'FontSize', 12);
        set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
        hold all
        plot(epochp,T6(:),'bo','linewidth',2,'markersize',8) %observation used for fitting
        plot(epochfit*yr,fit,'r.-','linewidth',2,'markersize',14) %linear fit
        plot(epochorgorg,demporg(:),'ko','linewidth',2,'markersize',4)%original observation
       hold on
        if 1
%       plot(epochorg,demp(:),'ks','linewidth',2,'markersize',10)
        shadedErrorBar(epochp,T6(:),T6std(:),'b',1)
        shadedErrorBar(epochfit*yr,fit,fitstd,'r',1)
        end
        %             errorbar(epoch,T6(:),T6std(:),'b.')
        %             errorbar(epochfit,fit,fitstd,'r.')

        minmax=max(abs(T6));
        meanT=mean(T6);dev=max(T6)-min(T6);
        uplmt=ceil(max(T6)+dev*0.3);lowerlmt=floor(min(T6)-dev*0.2);
        %plot shade area for eruption duration
        if flagvolc==1
            ha = area([eqepoch2(1) eqepoch2(2)], [uplmt uplmt],'basevalue',lowerlmt);
            ha.FaceColor='k';
            ha.FaceAlpha=0.5;
        end
        %repeat plot to bring lines to forward
        plot(epochp,T6(:),'bo','linewidth',2,'markersize',8) %observation used for fitting
        plot(epochfit*yr,fit,'r.-','linewidth',2,'markersize',14) %linear fit
      
        for k1=1:neq
            jumpt=trend;jumpstdt=trest;
            plot([eqepoch(k1),eqepoch(k1)],[lowerlmt uplmt],'k','Linewidth',2)
%             text(eqepoch(k1)+15,max(T6),texteq(k1),'FontSize',10)
        text(max(epochp)-(max(epochp)-min(epochp))*0.3,max(T6)-(max(T6)-min(T6))*0.2 ,['Jump=',num2str(jumpt,formatSpec),'\pm',num2str(jumpstdt,formatSpec),' m'],'FontSize',12)
        end
        text(max(epochp)-(max(epochp)-min(epochp))*0.3,max(T6)-(max(T6)-min(T6))*0.4 ,['Trend=',num2str(rate,formatSpec),'\pm',num2str(ratestd,formatSpec),' m/yr'],'FontSize',12)
%         plot(epochorg,demp(:),'go','linewidth',2,'markersize',14) %original 
        plot(epochorg(dempmt==0),demp(dempmt==0),'m+','linewidth',2,'markersize',14)%bad match points
        %outlier from previous iteration
        M1=ismember(epochorg,epochp); M2=ismember(demp,T6);
        Mkp=M1&M2;idoutM=~Mkp;
        plot(epochorg(idoutM),demp(idoutM),'ks','MarkerSize',12,'linewidth',4)% outliers from previous iteration
        %new outlier from this iteration
        M1=ismember(epochorg,epochp(idout2)); M2=ismember(demp,T6(idout2));
        idoutM2=M1&M2;gray=[0.5 0.5 0.5];green=[0 1 0];
        plot(epochorg(idoutM2),demp(idoutM2),'s','color',green,'MarkerSize',12,'linewidth',4)%new outlier from this iteration
        legend('ArcticDEM Elevation','Linear Fit','Unfiltered ArcticDEM Elevation')
        ylabel(['Elevation (m)'],'FontSize',12)
        set(gca,'XTickLabel',[])
        set(gca,'XTick',[])
        %             eptick=datenum({'2014/01/01','2014/07/01','2015/01/01','2015/07/01','2016/01/01','2016/07/01'})
        %             set(gca,'XTick',eptick)
        %             set(gca,'XTickLabel',eptick,'FontSize',12,'XTickLabelMode','manual')
        %             datetick('x','mm/yy','keepticks')
        datetick('x','mm/yy')
        box on
        %             grid on
        if (uplmt> lowerlmt)
        axis([min(epochp)-30 max(epochp)+30*4 lowerlmt uplmt])
        end
        ofile=['tp',num2str(jq)];
        title(['try event epoch i=',num2str(io2),';iter=',num2str(iter),';',num2str([jq])])
        saveas(gcf,ofile,'fig')
        %             print('-dpdf','-r300',ofile)    
        end % plotting
        %END of jump fitting
        
        end %iter
%         if (sigma02hat<max2);io2m=io2;max2=sigma02hat;end
        if (abs(trend)>max2);io2m=io2;max2=abs(trend);end
       
        end % for io2=1:no2

        eqepoch=eqepoch2(io2m);eqm=eqm2(io2m);eqs=eqs2{io2m};eqe=eqe2{io2m};
%         fprintf('Time series analysis ...')
%         toc
    %hypothesis test the statistical significance of the detected jump. %see outlierSchaffrin.m

	%Ian asks for average time for ice melting outputs.
	if algorithmin==1
	   eqm=mean(epoch)-epochs; %time in year after AD 0.
	end

end
