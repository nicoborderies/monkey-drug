function [ data ] = processForce_choiceRE( data, gripdata )

% detrend
g.grip = gripdata.grip(1,:);
[ ~ , g1 ] = grip_correctbaseline_choiceRE( g );
g.grip = gripdata.grip(2,:);
[ ~ , g2 ] = grip_correctbaseline_choiceRE( g );
gripdata.grip(1,:) = g1(:);
gripdata.grip(2,:) = g2(:);


% smoothing
    maxColumn = 125;
    if length(gripdata.grip)<height(data)
        gripdata.grip(:,height(data)) = {[];[]};
    end
    data.leftRawForce = gripdata.grip(1,:)';
    data.rightRawForce = gripdata.grip(2,:)';

    for iT = 1:numel(data.leftRawForce)
        data.time{iT,1}  = nan(1,maxColumn);
        
        if ~isempty(gripdata.grip{1,iT})
            time = gripdata.time{1,iT};  
            data.time{iT,1}(1:numel(time)) = (time-time(1));
            try
                [ leftForce , leftForceSpeed ]   = smoothForceSignal(data.leftRawForce(iT) ,{time-time(1)} ) ;
                data.leftForce(iT,1) = { leftForce };
                data.leftForceSpeed(iT,1) = { leftForceSpeed };
                [ rightForce , rightForceSpeed ]   = smoothForceSignal(data.rightRawForce(iT) ,{time-time(1)} ) ;
                data.rightForce(iT,1) = {  rightForce  };
                data.rightForceSpeed(iT,1) = { rightForceSpeed  };
            end
        end
    end
    

% reformat  
    maxLength1 = max(cellfun(@numel,data.leftForce));
    maxLength2 = max(cellfun(@numel,data.rightForce));
    maxLength = max([ maxLength1 maxLength2 ]);
    maxColumn = 125;
    data.chosenForce = cell(height(data),1);
    data.yankPeak = nan(height(data),1);
    data.responseTime =  nan(height(data),1);
    
    
     for iT = 1:numel(data.leftRawForce)
        data.chosenForce{iT,1} = nan(1,maxColumn);
        data.chosenPerf{iT,1} = data.chosenForce{iT,1};
        data.chosenPerfSpeed{iT,1} = data.chosenForce{iT,1};
        data.chosenPerfRT{iT,1} = data.chosenForce{iT,1};

        
        if ~isempty(gripdata.grip{1,iT})
            switch data.sideChoice(iT)
                case 1
                    nt = numel(data.leftRawForce{iT,1});
                    data.chosenForce{iT,1}(1:nt) = data.leftRawForce{iT,1};
                    data.chosenPerf{iT,1}(1:nt) = data.chosenForce{iT,1}(1:nt)./data.calibLeft(iT);
                    data.chosenPerfSpeed{iT,1}(1:nt-1) = data.leftForceSpeed{iT,1}(1:nt-1)./data.calibLeft(iT);

                case 2
                   nt = numel(data.rightRawForce{iT,1});
                   data.chosenForce{iT,1}(1:nt) = data.rightRawForce{iT,1};
                   data.chosenPerf{iT,1}(1:nt) = data.chosenForce{iT,1}(1:nt)./data.calibRight(iT);
                   data.chosenPerfSpeed{iT,1}(1:nt-1) = data.rightForceSpeed{iT,1}(1:nt-1)./data.calibRight(iT);
            end
        
        
        % correct estimate of peak force / response time
        if data.errorType(iT)==0 
            responseTimeCriterion = 0.50;
            force = data.chosenForce{iT,1}(1:nt);
            data.force(iT) = max(force);
            data.perf(iT) = max(data.chosenPerf{iT,1}(1:nt));
            speed = data.chosenPerfSpeed{iT,1}(1:nt-1);
            [data.yankPeak(iT)] = max(speed);
            i_t = find(speed >= responseTimeCriterion*max(speed),1,'first');
            data.responseTime(iT,1)  = data.time{iT,1}(i_t);
            i_t = find(force == max(force),1,'first');
            data.time2peak(iT,1)  = data.time{iT,1}(i_t) - data.responseTime(iT,1) ;


            data.chosenPerfRT{iT,1} = nan(1,maxColumn);
            data.chosenPerfRT{iT,1}(1:(nt-i_t+1))  =  data.chosenPerf{iT,1}(i_t:nt);
            
%             % plot signal smoothing
%                 clf('reset');
%                 hold on;
%                 x=data.time{iT,1};
%                 y=data.chosenPerf{iT,1};
%                 z=data.chosenPerfSpeed{iT,1};
%                 rt = data.responseTime(iT,1);
%                 
%                 h0 = plot(x,y,'k');
%                 h1 = plot(x,z./max(z),'Color',[1 1 1]*0.5);
%                 max_s = data.perf(iT).*ones(1,numel(x));
%                 max_ds = data.yankPeak(iT).*ones(1,numel(x));
%                 rtline = [0 1.2];
%                 h2 = plot(x,max_s,'Color','r');
%                 h4 = plot([rt rt],rtline,'Color','b');
%                 
%                 xlabel('time');
%                 ylabel('signal');
%                 title('smoothing');
%                 legend([h0,h1,h2,h4],{'s','ds_dt','max(s)','rt'})
%                 hold off;
%                 pause
            
        end
        

        end

    end

end

