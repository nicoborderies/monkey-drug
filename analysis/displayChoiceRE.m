function [fig] = displayChoiceRE(result,data)

%% define
    sub = unique(data.subject);
    sub = sub{1};

    %% participation figure
    participation = figure; set(participation,'Color',[1,1,1]);set(participation,'Name',[sub ':PARTICIPATION PATTERN (all trials)'])

        % histrogram
        subplot(2,2,[1 2]);hold on;
            hist(data.errorType);

            ax = gca; 
            set(ax,'XLim',[-1 5],'XTiCk',[0,1,2,3,4],'XTickLabel',...
                {'chosen(correct)';'unchosen(missed)';'unchosen(premature)';'unchosen(side switch)';'chosen(execution error)'},...
                  'TickLength',[ 0 0 ]);

            xlabel('trial category'); 
            ylabel('number of trials');
            
       
    %% participation figure
    participation = figure; set(participation,'Color',[1,1,1]);set(participation,'Name','PARTICIPATION PATTERN (all trials)')
    
        % histrogram
        subplot(2,2,[1]);hold on;
            empirDist = zeros(1,5);
            empirDist(unique(data.errorType)+1) = tools.tapply(data.errorType,{data.errorType},@numel);
            labels = {'chosen(correct)';'unchosen(missed)';'unchosen(premature)';'unchosen(side switch)';'chosen(execution error)'};
            pie(empirDist,labels);
            ax = gca; set(ax,'XColor','none'); set(ax,'YColor','none');

       subplot(2,2,[2]);hold on;
            selectionUnrepeated = (data.isRepeatedTrial==0 ) ;
            selectionUnmissed = ( data.errorType~=1);
            selection = ( selectionUnmissed & selectionUnrepeated ); 
            x = data.cardinalRewardLeft + data.cardinalRewardRight;
            x = round(x*100)/100;
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
            x = round(x*100)/100;
            y = tools.tapply(data.errorType~=1,{x},@nanmean);
            z = tools.tapply(data.errorType~=1,{x},@sem);
            errorbar(unique(x),y,z,'ko','LineWidth',2);
            plot(unique(x),y,'k--');
            ax = gca; %set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});
            xlabel('sum of offered efforts'); 
            ylabel('participation(%)');
            
            
    %% force figure
    force = figure; set(force,'Color',[1,1,1]);set(force,'Name',[sub ':FORCE PATTERN (chosen trials)'])

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
    choice = figure; set(choice,'Color',[1,1,1]);set(choice,'Name',[sub ':CHOICE PATTERN (correct trials)'])

        % lineplot
        subplot(2,2,1);hold on;
            X  = result.descriptive.quantile_sideChoice_deltaOrdinalReward;
            meanY = result.descriptive.mean_sideChoice_by_deltaOrdinalReward;
            errorY = result.descriptive.sem_sideChoice_by_deltaOrdinalReward;
            predictedY = result.inferential.prediction_sideChoice_by_deltaOrdinalReward;


            errorbar(X,meanY,errorY,'ro','LineWidth',2);
            plot(X,predictedY,'r--','LineWidth',2)

            ax = gca;
            set(ax,'YLim',[0 1]);%set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('delta Reward (right-left)'); 
            ylabel('right choice (%)');
            
            
        % lineplot
        subplot(2,2,2);hold on;
            X  = result.descriptive.quantile_sideChoice_deltaOrdinalEffort;
            meanY = result.descriptive.mean_sideChoice_by_deltaOrdinalEffort;
            errorY = result.descriptive.sem_sideChoice_by_deltaOrdinalEffort;
            predictedY = result.inferential.prediction_sideChoice_by_deltaOrdinalEffort;

            errorbar(X,meanY,errorY,'ro','LineWidth',2);
            plot(X,predictedY,'r--','LineWidth',2)

            ax = gca;
            set(ax,'YLim',[0 1]);%set(ax,'XTiCk',[1,2,3,4],'XTickLabel',{'50';'65';'75';'85'});

            xlabel('delta Effort (right-left)'); 
            ylabel('right choice (%)');
            
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
            


fig.participation = participation;
fig.force = force;
fig.choice = choice;

end
