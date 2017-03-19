function [fig] = displayChoiceRE_2(data,result)

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
        selectionResponseTime =  ~selectionMissed & ~isnan(data.intertrialTime);
        selectionChosen = ( ~selectionMissed & ~selectionPremature & ~selectionSwitch  );
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

%% Trial classification
    fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','PARTICIPATION PATTERN');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.15 0.8 0.4];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram pie chart
        subplot(1,numel(subjectList),iSub);hold on
        errors = data.errorType(selectSubject & ~selectionRepetition);
        empirDist = zeros(1,6);
        empirDist(unique(errors)+1) = tapply(errors,{errors},@numel);
        labels = {'chosen(correct)';'unchosen(missed)';'unchosen(premature)';'';'chosen(undershoot error)';'chosen (overshoot error)'};
        p = pie(empirDist,labels);
        ax = gca; 
        ax.XColor = 'w'; ax.YColor = 'w';
        ax.XLim = [-1.5 1.5 ];
        ax.YLim = [-1 1.5 ];
        xlabel('trial classification'); 
        t = title(['monkey ' sub(1) ]); t.Position(2) = 1.3;
    end
    
    
%% Trial Classes conditionnal frequency
    fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','PARTICIPATION PATTERN 2');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.15 0.8 0.4];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram pie chart
        subplot(1,numel(subjectList),iSub);hold on
        x = selectionRepetition(selectSubject);
        repeat = mean(x);
        e1 = sem(x);
        x = selectionParticipate(selectSubject & ~selectionRepetition);
        participation = mean(x);
        e2 = sem(x);
        x = selectionPremature(selectSubject & selectionParticipate & ~selectionRepetition);
        premature = mean(x);
        e3 = sem(x);
        x = selectionCorrect(selectSubject & selectionChosen & ~selectionRepetition);
        correctExecution = mean(x);
        e4 = sem(x);

        b = bar([1 2 3 4],[ repeat participation premature correctExecution ] );
        h = errbar([1 2 3 4],[ repeat participation premature correctExecution ],[ e1 e2 e3 e4 ],...
                 'Color',[0.5 0.5 0.5],'LineWidth',2);
        b.FaceColor = [ 0.3 0.75 0.9 ];
        b.EdgeColor = 'none';
%         h.LineStyle = '-';
%         h.Color = [0.5 0.5 0.5];
%         h.LineWidth = 2;
        
%         h.Color = [0.5 0.5 0.5];

        ax = gca; 
        ax.YLim = [0 1];
        ax.TickLength = [0 0]; 
        ax.XTick = [ 1 2 3 4];
        ax.XTickLabel = {'P(repetition)','P(participation)','P(prematrue | participation)','P(correct | chosen)'};
        ax.XTickLabelRotation = 45; 
        t = title(['monkey ' sub(1) ]); t.Position(2) = 1;
    end
    
    
%% Participation functions
    fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','PARTICIPATION PATTERN 2');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.05 0.8 0.8];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % image plot
        subplot(2,numel(subjectList),iSub);hold on
        participation = selectionParticipate(selectSubject & ~selectionRepetition);
        p = tapply(participation,{data.maxReward(selectSubject & ~selectionRepetition),...
                                    data.minReward(selectSubject & ~selectionRepetition)},@nanmean);
        cm = colormap('hot'); cm = flipud(cm);
        minCM = min(min(p))-0.3 ; maxCM = max(max(p))+0.3;
        heatmap(p,[1:4],[1:4],fliplr(p),'Colormap',cm,'NaNColor',[1 1 1],...
            'MaxColorValue',maxCM,'MinColorValue',minCM,'TextColor',[1 1 1]*0);
        colorbar;
        ax = gca; 
        xlabel('min(offered rewards)');
        ylabel('max(offered rewards)');
        t = title(['monkey ' sub(1) ]); t.Position(2) = 4.6;
        
        % image plot
        subplot(2,numel(subjectList),iSub+numel(subjectList));hold on
        participation = selectionParticipate(selectSubject & ~selectionRepetition);
        p = tapply(participation,{data.maxEffort(selectSubject & ~selectionRepetition),...
                                    data.minEffort(selectSubject & ~selectionRepetition)},@nanmean);
         cm = colormap('hot'); cm = flipud(cm);
        heatmap(p,[1:4],[1:4],fliplr(p),'Colormap',cm,'NaNColor',[1 1 1],...
            'MaxColorValue',1.3,'MinColorValue',-0.3,'TextColor',[1 1 1]*0);
        colorbar;
        ax = gca; 
        xlabel('min(offered efforts)');
        ylabel('max(offered efforts)');
        t = title(['monkey ' sub(1) ]); t.Position(2) = 4.6;
    end    
    
    
%% Participation functions
    fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','PARTICIPATION PATTERN 2');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.05 0.8 0.8];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % image plot
        participation = selectionParticipate(selectSubject & ~selectionRepetition);
        nt = data.trialNumber(selectSubject & ~selectionRepetition);
        cumRew  =   data.cumulativeReward(selectSubject & ~selectionRepetition); 
        cumPerf =  data.cumulativeEffort(selectSubject & ~selectionRepetition); 

        
        subplot(3,numel(subjectList),iSub);hold on
        i = tapply(nt,{nt},@nanmean,{'continuous'},10);
        p = tapply(participation,{nt},@nanmean,{'continuous'},10);
        e = tapply(participation,{nt},@sem,{'continuous'},10);
        [h,hp] = boundedline( i , p , e , 'alpha' );   
        set(h,'Color','k','LineWidth',2);set(hp,'FaceColor','k');
        h.LineStyle = '-'; 
        scatter( i , p , 'filled','k' );
        ax = gca; 
        xlabel('trial number');
        ylabel('participation (%)');
        t = title(['monkey ' sub(1) ]); 
        
        
        subplot(3,numel(subjectList),iSub+numel(subjectList));hold on
        i = tapply(cumRew,{cumRew},@nanmean,{'continuous'},10);
        p = tapply(participation,{cumRew},@nanmean,{'continuous'},10);
        e = tapply(participation,{cumRew},@sem,{'continuous'},10);
        [h,hp] = boundedline( i , p , e , 'alpha' );   
        set(h,'Color','g','LineWidth',2);set(hp,'FaceColor','g');
        h.LineStyle = '-'; 
        scatter( i , p , 'filled','g' );
        ax = gca; 
        xlabel('cumulative reward (ml)');
        ylabel('participation (%)');
        
        subplot(3,numel(subjectList),iSub+2*numel(subjectList));hold on
        i = tapply(cumPerf,{cumPerf},@nanmean,{'continuous'},10);
        p = tapply(participation,{cumPerf},@nanmean,{'continuous'},10);
        e = tapply(participation,{cumPerf},@sem,{'continuous'},10);
        [h,hp] = boundedline( i , p , e , 'alpha' );   
        set(h,'Color','r','LineWidth',2);set(hp,'FaceColor','r');
        h.LineStyle = '-'; 
        scatter( i , p , 'filled','r' );
        ax = gca; 
        xlabel('cumulative force (%fmax)');
        ylabel('participation (%)');
        
    end    

%% Force pattern
%  fig = figure; set(fig,'Color',[1,1,1]);
%     set(fig,'Name','PARTICIPATION PATTERN 2');
%     fig.Units = 'normalized';
%     fig.Position = [0.15 0.05 0.8 0.4];
%     for iSub = 1:numel(subjectList)
%         sub = subjectList{iSub};
%         selectSubject = ismember(data.subject, sub);
%         
%         % line plots
%         force = data.perf(selectSubject & selectionChosen);
%         dynamic = cell2mat(data.chosenPerfRT(selectSubject & selectionChosen,:));
%         choice = data.sideChoice(selectSubject & selectionChosen);
%         effort = data.chosenCardinalEffort(selectSubject & selectionChosen);
%         reward = data.chosenOrdinalReward(selectSubject & selectionChosen);
%         freq = 25;
%         
%         
%         subplot(1,numel(subjectList),iSub);hold on
%         i = unique(effort);
%         it=0;
%         for level = i'
%             it = it+1;
%             d = nanmean(dynamic(effort==level,:),1);
%             e = sem(dynamic(effort==level,:),1);
%             tmax = find(isnan(d),1,'first')-30;
%             [h(it),hp] = boundedline( [1:numel(d(1:tmax))].*(1/freq) , d(1:tmax) , e(1:tmax) , 'alpha' );   
%             col = [1-0.3*level 1-level 0.5*(1-level)];
%             set(h(it),'Color',col,'LineWidth',2);set(hp,'FaceColor',col);
%             h(it).LineStyle = '-'; 
%         end
%         legend([h(1) h(2) h(3) h(4)],'10%','40%','70%','100%');
%         ax = gca; 
%         xlabel('time (sec.)'); 
%         ylabel('force dynamic (% calibration)');
%         t = title(['monkey ' sub(1) ]);         
%      
%     end    

    
%% Force pattern
 fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','FORCE PATTERN 2');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.05 0.8 0.8];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % line plots
        force = data.perf(selectSubject & selectionChosen);
        correct =  selectionCorrect(selectSubject & selectionChosen);
        choice = data.sideChoice(selectSubject & selectionChosen);
        effort = data.chosenCardinalEffort(selectSubject & selectionChosen);        
        
        i = unique(effort);
        it=0;
        for level = i'
            it = it+1;
            subplot(4,numel(subjectList),iSub+(it-1)*numel(subjectList));hold on
            
            % correct
            f = force(effort==level & correct);
            l = repmat(0.05,numel(f),1);
            h(1) = histogram(f,'BinWidth',0.05);
            h(1).FaceColor = 'b' ;
            lim = max(histcounts(f,'BinWidth',0.05));
            
            % incorrect
            f = force(effort==level & ~correct);
            l = repmat(level,numel(f),1);
            h(2) = histogram(f,'BinWidth',0.05);
            h(2).FaceColor = [0.5 0.5 0.5] ;
            lim = max([histcounts(f,'BinWidth',0.05),lim]);
            
            % level separation
            plot([level level],[0 lim],'r','LineWidth',2);
            ax = gca; 
            ax.YLim = [0 lim+1];
            ax.XLim = [0 1.5];
            ax.TickLength = [0 0];
            ax.XTick = [];
            if it==1
                t = title(['monkey ' sub(1) ]); 
            elseif it==4
                ax.XTick = [2:2:14]*0.1;
                xlabel('force peak (% calibration)');
            end
            if iSub==1
                ylabel(['effort = '...
                num2str(round(level*100)) ' %']);
            end
        end
    end    
    legend([h(1) h(2)],'correct','incorrect');


%% Force distributions 

 fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','FORCE PATTERN 2');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.05 0.8 0.8];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % line plots
        force = data.perf(selectSubject & selectionChosen);
        choice = data.sideChoice(selectSubject & selectionChosen);
        effort = data.chosenCardinalEffort(selectSubject & selectionChosen);
        reward = data.chosenOrdinalReward(selectSubject & selectionChosen);
        
        
        subplot(2,numel(subjectList),iSub);hold on
        i = tapply(effort,{choice,effort},@nanmean); 
        p = tapply(force,{choice,effort},@nanmean);
        e = tapply(force,{choice,effort},@sem);
        h(1) = errorbar(i(1,:),p(1,:),e(1,:),'c-o','LineWidth',2);
        h(1).MarkerFaceColor = 'c';
        h(2) = errorbar(i(2,:),p(2,:),e(2,:),'m-o','LineWidth',2);
        h(2).MarkerFaceColor = 'm';
        plot(i(1,:),i(1,:),'k--');
        legend([h(1) h(2)],'left','right')
        ax = gca; 
        ax.XLim = [0 1.1];
        ax.XTick = i(1,:);
        xlabel('chosenEffort(% calibration)'); 
        ylabel('force (% calibration)');
        t = title(['monkey ' sub(1) ]); 
        
        subplot(2,numel(subjectList),iSub+numel(subjectList));hold on
        i = unique(reward);
        p = tapply(force,{reward},@nanmean);
        e = tapply(force,{reward},@sem);
        h = errorbar(i,p,e,'k-o','LineWidth',2);
        h.MarkerFaceColor = 'k';
        ax = gca; 
        ax.XLim = [0 4.5];
        ax.XTick = i';
        xlabel('chosenReward'); 
        ylabel('force (% calibration)');
        
        
    end                

 %% Choice pattern
 
 fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','CHOICE PATTERN');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.05 0.8 0.8];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % line plots
        force = data.perf(selectSubject & selectionChosen);
        choice = data.sideChoice(selectSubject & selectionChosen)-1;
        choiceR = (data.sideChoice(selectSubject & selectionChosen)...
                    == iMaxReward(selectSubject & selectionChosen));
        dR = abs(data.deltaCardinalReward(selectSubject & selectionChosen));
        DR = (data.deltaCardinalReward(selectSubject & selectionChosen));
        DR = round(DR);

        minR = data.minReward(selectSubject & selectionChosen);
        choiceE = (data.sideChoice(selectSubject & selectionChosen)...
                    == iMaxEffort(selectSubject & selectionChosen));
        dE = abs(data.deltaCardinalEffort(selectSubject & selectionChosen));
        DE = (data.deltaCardinalEffort(selectSubject & selectionChosen));
        DE = round(DE,1);
        minE = data.minEffort(selectSubject & selectionChosen);

        
        
        subplot(3,numel(subjectList),iSub);hold on
        ind = unique(dR);
        p = tapply(choiceR,{minR,dR},@nanmean);
        e = tapply(choiceR,{minR,dR},@sem);
        it=0;
        for level = ind(1:end-1)'
            it = it+1;
            col = [0.7-0.15*it 1-0.15*it 0.7-0.15*it];
            h(it) = errorbar(ind(2:end),p(it,2:end),e(it,2:end),'g-o','LineWidth',2);
            h(it).MarkerFaceColor = col ; h(it).Color = col;
        end
        legend([h(1) h(2) h(3)],'minR=1','minR=2','minR=3');

        ax = gca; 
        ax.XLim = [0 4];
        ax.XTick = ind';
        xlabel('delta reward'); 
        ylabel('max reward choice (%)');
        t = title(['monkey ' sub(1) ]); 

        
        subplot(3,numel(subjectList),iSub+numel(subjectList));hold on
        dE = round(dE,1);
        ind = unique(dE);
        p = tapply(choiceE,{minE,dE},@nanmean);
        e = tapply(choiceE,{minE,dE},@sem);
        it=0;
        for level = ind(1:end-1)'
            it = it+1;
            col = [1-0.15*it 0.7-0.15*it 0.7-0.15*it];
            h(it) = errorbar(ind(2:end),p(it,2:end),e(it,2:end),'r-o','LineWidth',2);
            h(it).MarkerFaceColor = col ; h(it).Color = col;
        end
        legend([h(1) h(2) h(3)],'minE=1','minE=2','minE=3');
        ax = gca; 
        ax.XLim = [0 1.1];
        ax.XTick = ind';
        xlabel('delta effort'); 
        ylabel('max effort choice (%)');
        
        
         % image plot
        subplot(3,numel(subjectList),iSub+2*numel(subjectList));hold on
        p = tapply(choice,{DR,DE},@nanmean);
        cm = colormap('hot'); cm = flipud(cm);
        minCM = min(min(p))-0.3 ; maxCM = max(max(p))+0.3;
        heatmap(p,unique(DE),unique(DR),fliplr(p),'Colormap',cm,'NaNColor',[1 1 1],...
            'MaxColorValue',maxCM,'MinColorValue',minCM,'TextColor',[1 1 1]*0);
        colorbar;
        ax = gca; 
        ax.XTick = [1:numel(unique(DE))];
        ax.XTickLabel = {'-0.9','-0.6','-0.3','0','+0.3','+0.6','+0.9'};
        ax.YTick = [1:numel(unique(DR))];
        ax.YTickLabel = {'-3','-2','-1','0','+1','+2','+3'};
        xlabel(texlabel('delta effort (right-left)'));
        ylabel(texlabel('delta reward (right-left)'));
        title('right choice (%)');


    end    
 
%% CHoice Pattern 2
 
 fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','CHOICE PATTERN');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.05 0.8 0.8];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % line plots
        nt = data.trialNumber(selectSubject & selectionChosen);
        choiceR = (data.sideChoice(selectSubject & selectionChosen)...
                    == iMaxReward(selectSubject & selectionChosen));
        choiceE = (data.sideChoice(selectSubject & selectionChosen)...
                    == iMaxEffort(selectSubject & selectionChosen));
        
        subplot(2,numel(subjectList),iSub);hold on
        i = tapply(nt,{nt},@nanmean,{'continuous'},10);
        p = tapply(choiceR,{nt},@nanmean,{'continuous'},10);
        e = tapply(choiceR,{nt},@sem,{'continuous'},10);
        [h,hp] = boundedline( i , p , e , 'alpha' );   
        col = 'g';
        set(h,'Color',col,'LineWidth',2);set(hp,'FaceColor',col);
        h.LineStyle = '-'; 
        scatter( i , p , 'filled',col );
        ax = gca; 
        xlabel('trial number'); 
        ylabel('max reward choice (%)');
        t = title(['monkey ' sub(1) ]); 

        
        subplot(2,numel(subjectList),iSub+numel(subjectList));hold on
        i = tapply(nt,{nt},@nanmean,{'continuous'},10);
        p = tapply(choiceE,{nt},@nanmean,{'continuous'},10);
        e = tapply(choiceE,{nt},@sem,{'continuous'},10);
        [h,hp] = boundedline( i , p , e , 'alpha' );   
        col = 'r';
        set(h,'Color',col,'LineWidth',2);set(hp,'FaceColor',col);
        h.LineStyle = '-'; 
        scatter( i , p , 'filled', col );
        ax = gca; 
        xlabel('trial number'); 
        ylabel('max effort choice (%)');
        
    end   

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



end