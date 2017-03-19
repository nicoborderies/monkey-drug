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
        selectionIncorrect = ( data.errorType==4);

        selectionParticipate = ~selectionMissed ;
        selectionChosen = ( ~selectionRepetition & ~selectionMissed & ~selectionPremature & ~selectionSwitch  );
        selectionCorrect = ( selectionChosen  & ~selectionIncorrect );
        
    % subjectList 
    subjectList = unique(data.subject);

%% participation figure
    fig = figure; set(fig,'Color',[1,1,1]);
    set(fig,'Name','PARTICIPATION PATTERN (all trials)');
    fig.Units = 'normalized';
    fig.Position = [0.15 0.15 0.8 0.5];
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        % histrogram pie chart
        subplot(1,numel(subjectList),iSub);hold on
        empirDist = zeros(1,5);
        empirDist(unique(data.errorType)+1) = tapply(data.errorType,{data.errorType},@numel);
        labels = {'chosen(correct)';'unchosen(missed)';'unchosen(premature)';'unchosen(side switch)';'chosen(execution error)'};
        pie(empirDist,labels);
        ax = gca; 
        xlabel('trial classification'); 
        title(['monkey ' sub(1) ]); 
        
    end

%%

            
       %
       subplot(2,2,[2]);hold on;
            selectionUnrepeated = (data.isRepeatedTrial==0 ) ;
            selectionUnmissed = ( data.errorType~=1);
            selection = ( selectionUnmissed & selectionUnrepeated ); 
            x = data.cardinalRewardLeft + data.cardinalRewardRight;
            x = x(selection);
            y = tools.tapply(data.errorType(selection)==2,{x},@nanmean);
            z = tools.tapply(data.errorType(selection)==2,{x},@sem);

            errorbar(unique(x),y,z,'ko','LineWidth',2);
            plot(unique(x),y,'k--');

            ax = gca; %set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('sum of offered rewards'); 
            ylabel('premature response(%)');
            

     subplot(2,2,[3]);hold on;
            x = data.cardinalRewardLeft + data.cardinalRewardRight; 
            y = tools.tapply(data.errorType~=1,{x},@nanmean);
            z = tools.tapply(data.errorType~=1,{x},@sem);

            errorbar(unique(x),y,z,'ko','LineWidth',2);
            plot(unique(x),y,'k--');

            ax = gca; %set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('sum of offered rewards'); 
            ylabel('participation(%)');
            
      subplot(2,2,[4]);hold on;
            selectionUnrepeated = (data.isRepeatedTrial==0 ) ;
            selectionUnmissed = ( data.errorType~=1);
            selection = ( selectionUnmissed & selectionUnrepeated ); 
            x = data.cardinalEffortLeft + data.cardinalEffortRight; 
            y = tools.tapply(data.errorType~=1,{x},@nanmean);
            z = tools.tapply(data.errorType~=1,{x},@sem);

            errorbar(unique(x),y,z,'ko','LineWidth',2);
            plot(unique(x),y,'k--');

            ax = gca; %set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('sum of offered efforts'); 
            ylabel('participation(%)');
       

            
    %% force figure
    fig = figure; set(fig,'Color',[1,1,1]);set(fig,'Name','FORCE PATTERN (chosen trials)')

        % lineplot
        subplot(3,2,1);hold on;
            X  = result.descriptive.quantile_perf_chosenEffort;
            meanY = result.descriptive.mean_perf_by_chosenEffort;
            errorY = result.descriptive.sem_perf_by_chosenEffort;
            predictedY = result.inferential.prediction_perf_by_chosenEffort;


            errorbar(X,meanY,errorY,'bo','LineWidth',2);
            plot(X,predictedY,'b--','LineWidth',2);
            plot(X,X,'k--');

            ax = gca; %set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('chosenEffort'); 
            ylabel('force (% calibration)');
            
            
        % lineplot
        subplot(3,2,2);hold on;
            X  = result.descriptive.quantile_perf_chosenReward;
            meanY = result.descriptive.mean_perf_by_chosenReward;
            errorY = result.descriptive.sem_perf_by_chosenReward;

            errorbar(X,meanY,errorY,'bo','LineWidth',2);

            ax = gca; %set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('chosenReward'); 
            ylabel('force  (% calibration)');
            
         % lineplot
        subplot(3,2,3);hold on;
            X  = result.descriptive.quantile_cumulativePerf_chosenEffort;
            meanY = result.descriptive.mean_cumulativePerf_by_chosenEffort;
            errorY = result.descriptive.sem_cumulativePerf_by_chosenEffort;

            errorbar(X,meanY,errorY,'bo','LineWidth',2);
            plot(X,X,'k--');

            ax = gca; %set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('chosenEffort'); 
            ylabel('force cumulée (% calibration)');
            
            
        % lineplot
        subplot(3,2,4);hold on;
            X  = result.descriptive.quantile_cumulativePerf_chosenReward;
            meanY = result.descriptive.mean_cumulativePerf_by_chosenReward;
            errorY = result.descriptive.sem_cumulativePerf_by_chosenReward;

            errorbar(X,meanY,errorY,'bo','LineWidth',2);

            ax = gca; %set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('chosenReward'); 
            ylabel('force cumulée(% calibration)');
            
            
        % barplot
        subplot(3,2,5);hold on;
            X  = result.descriptive.quantile_executionError_sideChoice;
            meanY = result.descriptive.mean_executionError_by_sideChoice ; 
            errorY = result.descriptive.sem_executionError_by_sideChoice ;

            bar(X,meanY,'b','LineWidth',2,'EdgeColor','k');
            errorbar(X,meanY,errorY,'k','LineWidth',2,'LineStyle','none');

            ax = gca;
            set(ax,'YLim',[0 1]);
            set(ax,'XTiCk',[1,2],'XTickLabel',{'left';'right'});

            xlabel('chosen side'); 
            ylabel('execution error (%)');
            
       % barplot
        subplot(3,2,6);hold on;
            X  = result.descriptive.quantile_perf_sideChoice;
            meanY = result.descriptive.mean_perf_by_sideChoice ; 
            errorY = result.descriptive.sem_perf_by_sideChoice ;

            bar(X,meanY,'b','LineWidth',2,'EdgeColor','k');
            errorbar(X,meanY,errorY,'k','LineWidth',2,'LineStyle','none');

            ax = gca;
%             set(ax,'YLim',[0 1]);
            set(ax,'XTiCk',[1,2],'XTickLabel',{'left';'right'});

            xlabel('chosen side'); 
            ylabel('force  (% calibration)');

 %% choice figure
    fig = figure; set(fig,'Color',[1,1,1]);set(fig,'Name','CHOICE PATTERN (correct trials)')

        % lineplot
        subplot(2,2,1);hold on;
            selectionUnrepeated = (data.isRepeatedTrial==0 ) ;
            selectionUnmissed = ( data.errorType~=1);
            selectionChosen = ( selectionUnmissed & data.errorType~=2  & data.errorType~=3  );
            selectionCorrect = ( selectionChosen  & selectionUnrepeated & data.errorType~=4  );
            selectionDR =  (data.deltaCardinalReward~=0);
            selection = (selectionChosen & selectionDR);
            
            design = data(selection,:);
            rewards = [design.cardinalRewardLeft, design.cardinalRewardRight];
           [xMax,iMaxReward] = max(rewards,[],2);
           [xMin,iMinReward] = min(rewards,[],2);
           x = xMax - xMin;
           
           y = (design.sideChoice==iMaxReward);
           yy = (design.predictedChoice==iMaxReward);
 
           X = tools.tapply(x,{x},@nanmean);
           Y = tools.tapply(y,{x},@nanmean);
           errorY = tools.tapply(y,{x},@sem);
%            predictedY = tools.tapply(yy,{x},@nanmean);
        
%             X  = result.descriptive.quantile_sideChoice_deltaOrdinalReward;
%             meanY = result.descriptive.mean_sideChoice_by_deltaOrdinalReward;
%             errorY = result.descriptive.sem_sideChoice_by_deltaOrdinalReward;
%             predictedY = result.inferential.prediction_sideChoice_by_deltaOrdinalReward;
%             errorbar(X,meanY,errorY,'ro','LineWidth',2);
%             plot(X,predictedY,'r--','LineWidth',2);
%             
            errorbar(X,Y,errorY,'go','LineWidth',2,'MarkerFaceColor','g');
            plot(X,Y,'g--','LineWidth',2);
            
%             plot(X,predictedY,'g--','LineWidth',2)

            ax = gca;
            set(ax,'YLim',[0 1]);%set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('max(reward) - min(reward)'); 
%             ylabel('right choice (%)');
            ylabel(' max reward choice (%)');

            
        % lineplot
        subplot(2,2,2);hold on;
%             X  = result.descriptive.quantile_sideChoice_deltaOrdinalEffort;
%             meanY = result.descriptive.mean_sideChoice_by_deltaOrdinalEffort;
%             errorY = result.descriptive.sem_sideChoice_by_deltaOrdinalEffort;
%             predictedY = result.inferential.prediction_sideChoice_by_deltaOrdinalEffort;
% 
%             errorbar(X,meanY,errorY,'ro','LineWidth',2);
%             plot(X,predictedY,'r--','LineWidth',2)
% 
%             ax = gca;
%             set(ax,'YLim',[0 1]);%set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});
% 
%             xlabel('max(reward) - min(reward)'); 
%             ylabel('right choice (%)');
            
            selectionUnrepeated = (data.isRepeatedTrial==0 ) ;
            selectionUnmissed = ( data.errorType~=1);
            selectionChosen = ( selectionUnmissed & data.errorType~=2  & data.errorType~=3  );
            selectionCorrect = ( selectionChosen  & selectionUnrepeated & data.errorType~=4  );
            selectionDE =  (data.deltaCardinalEffort~=0);
            selection = (selectionChosen & selectionDE);
            
            design = data(selection,:);
            xx = [design.cardinalEffortLeft, design.cardinalEffortRight];
           [xMax,iMax] = max(xx,[],2);
           [xMin,iMin] = min(xx,[],2);
           x = xMax - xMin;
           x = (round(x*10)/10);
           
           y = (design.sideChoice==iMax);
           yy = (design.predictedChoice==iMax);
 
           X = tools.tapply(x,{x},@nanmean);
           Y = tools.tapply(y,{x},@nanmean);
           errorY = tools.tapply(y,{x},@sem);
            errorbar(X,Y,errorY,'ro','LineWidth',2,'MarkerFaceColor','r');
            plot(X,Y,'r--','LineWidth',2);
           ax = gca;
            set(ax,'YLim',[0 1]);%set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('max(effort) - min(effort)'); 
            ylabel(' max effort choice (%)');
            
            
         % barplot
        subplot(2,2,3);hold on;
            X  = result.descriptive.quantile_sideChoice_previousSideChoice + 1;
            meanY = result.descriptive.mean_sideChoice_by_previousSideChoice; meanY(1) = 1 - meanY(1);
            errorY = result.descriptive.sem_sideChoice_by_previousSideChoice;

            bar(X,meanY,'r','LineWidth',2,'EdgeColor','k');
            errorbar(X,meanY,errorY,'k','LineWidth',2,'LineStyle','none');

            ax = gca;
            set(ax,'YLim',[0 1]);
            set(ax,'XTiCk',[1,2],'XTickLabel',{'left';'right'});

            xlabel('previous choice'); 
            ylabel('choice repetition (%)');
            



end