function displayRT(data)
%% define
   % data acquisition parameters
        forceSamplingFreq = 25; % Hz.
        rewardUnit2Volume = 0.9; % conversion rate between valve opening time (sec) & water volume (ml)

    % data selection
        selectionRepetition = (data.isRepeatedTrial==1) ;
        
        selectionMissed = ( data.errorType==1);
        selectionPremature = ( data.errorType==2);
        selectionSwitch = ( data.errorType==3);
        selectionIncorrect = ( data.errorType==4 | data.errorType==5 );

        selectionParticipate = ~selectionMissed ;
        sessionRT = nominal({'11_18_2015'});
        selectionResponseTime = ~selectionRepetition & ~selectionMissed & ~isnan(data.intertrialTime);
        selectionChosen = ( ~selectionRepetition & ~selectionMissed & ~selectionPremature & ~selectionSwitch  );
        selectionCorrect = ( selectionChosen  & ~selectionIncorrect );
        
     % data completion
        data.volumeRewardLeft = data.cardinalRewardLeft.*data.rewardUnit.*rewardUnit2Volume;
        data.volumeRewardRight= data.cardinalRewardRight.*data.rewardUnit.*rewardUnit2Volume;
        ind = find(ismember(data.Properties.VariableNames,'gain'));
        data.Properties.VariableNames{ind} = 'ordinalRewardOutcome';
        
    
        data.chosenOrdinalReward = nan(numel(data.trialNumber),1);
        data.chosenOrdinalReward(selectionChosen & data.sideChoice==1) = data.ordinalRewardLeft(selectionChosen & data.sideChoice==1);
        data.chosenOrdinalReward(selectionChosen & data.sideChoice==2) = data.ordinalRewardRight(selectionChosen & data.sideChoice==2);
        data.chosenOrdinalEffort = nan(numel(data.trialNumber),1);
        data.chosenOrdinalEffort(selectionChosen & data.sideChoice==1) = data.ordinalEffortLeft(selectionChosen & data.sideChoice==1);
        data.chosenOrdinalEffort(selectionChosen & data.sideChoice==2) = data.ordinalEffortRight(selectionChosen & data.sideChoice==2);

        data.chosenCardinalReward = nan(numel(data.trialNumber),1);
        data.chosenCardinalReward(selectionChosen & data.sideChoice==1) = data.cardinalRewardLeft(selectionChosen & data.sideChoice==1);
        data.chosenCardinalReward(selectionChosen & data.sideChoice==2) = data.cardinalRewardRight(selectionChosen & data.sideChoice==2);
        data.chosenCardinalEffort = nan(numel(data.trialNumber),1);
        data.chosenCardinalEffort(selectionChosen & data.sideChoice==1) = data.cardinalEffortLeft(selectionChosen & data.sideChoice==1);
        data.chosenCardinalEffort(selectionChosen & data.sideChoice==2) = data.cardinalEffortRight(selectionChosen & data.sideChoice==2);

        data.chosenVolumeReward  = data.chosenCardinalReward.*data.rewardUnit.*rewardUnit2Volume;
        data.volumeRewardOutcome = data.chosenVolumeReward.*(data.errorType==0); 

        
        data.deltaOrdinalReward = data.ordinalRewardRight - data.ordinalRewardLeft;
        data.deltaCardinalReward = data.cardinalRewardRight - data.cardinalRewardLeft;
        
        data.deltaOrdinalEffort = data.ordinalEffortRight - data.ordinalEffortLeft;
        data.deltaCardinalEffort = data.cardinalEffortRight - data.cardinalEffortLeft;
        
        [data.maxReward,iMaxReward] = max([data.ordinalRewardLeft , data.ordinalRewardRight ],[],2);
        [data.minReward,iMinReward] = min([data.ordinalRewardLeft , data.ordinalRewardRight ],[],2);
        [data.maxEffort,iMaxEffort] = max([data.ordinalEffortLeft , data.ordinalEffortRight ],[],2);
        [data.minEffort,iMinEffort] = min([data.ordinalEffortLeft , data.ordinalEffortRight ],[],2);

        % subjectList 
        subjectList = unique(data.subject);
        
        
        data.cumulativeReward = nan(height(data),1);
        data.cumulativeEffort = nan(height(data),1);
        for iSub = 1:numel(subjectList)
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            sessionList = unique(data.session(selectSubject));
            for iSess = 1:numel(sessionList)
                selectSession = (data.session==sessionList(iSess)) ;
                cumRew =   nancumsum(data.volumeRewardOutcome(selectSubject & selectSession & ~selectionRepetition)); 
                cumRew(2:end) = cumRew(1:end-1);cumRew(1)= 0;
                data.cumulativeReward(selectSubject & selectSession & ~selectionRepetition) = cumRew;
                cumPerf =  nancumsum(data.perf(selectSubject & selectSession & ~selectionRepetition)); 
                cumPerf(2:end) = cumPerf(1:end-1);cumPerf(1)= 0;
                data.cumulativeEffort(selectSubject & selectSession & ~selectionRepetition) = cumPerf;

            end
        end
        
    % subjectList 
    subjectList = unique(data.subject);
    
    
%% Response Times
 
 fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','CHOICE PATTERN');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.05 0.8 0.8];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % line plots
        rt = data.reactionTime(selectSubject & selectionResponseTime);
        offer = data.offerTime(selectSubject & selectionResponseTime);
        effort = data.chosenCardinalEffort(selectSubject & selectionResponseTime);        
        reward = data.chosenOrdinalReward(selectSubject & selectionResponseTime);        
        
        subplot(3,numel(subjectList),iSub);hold on
        h = scatter(offer,rt,'filled');
        h.MarkerFaceColor = 'k';
        plot(offer,offer,'k--');
        ax = gca; 
        xlabel('offer time (sec)'); 
        ylabel('response time (sec)');
        t = title(['monkey ' sub(1) ]); 
        
        subplot(3,numel(subjectList),iSub+numel(subjectList));hold on
        i = unique(reward);
        p = tapply(rt,{reward},@nanmean);
        e = tapply(rt,{reward},@sem);
        h = errorbar(i,p,e,'g-o','LineWidth',2);
        h.MarkerFaceColor = 'g';
        ax = gca; 
        ax.XLim = [0 4.5];
        ax.XTick = i(~isnan(i))';
        xlabel('chosen reward'); 
        ylabel('response time (sec)');
        
        subplot(3,numel(subjectList),iSub+2*numel(subjectList));hold on
        i = unique(effort);
        p = tapply(rt,{effort},@nanmean);
        e = tapply(rt,{effort},@sem);
        h = errorbar(i,p,e,'r-o','LineWidth',2);
        h.MarkerFaceColor = 'r';
        ax = gca; 
        ax.XLim = [0 1.1];
        ax.XTick = i(~isnan(i))';
        xlabel('chosen effort'); 
        ylabel('response time (sec)');
        
    end       
