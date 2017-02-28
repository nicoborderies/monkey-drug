function [result,data] = processChoiceRE(data)
%% define
    % data acquisition parameters
        forceSamplingFreq = 25; % Hz.

    % data selection
        selectionUnrepeated = (data.isRepeatedTrial==0 ) ;
        selectionUnmissed = ( data.errorType~=1);
        selectionChosen = ( selectionUnmissed & data.errorType~=2  & data.errorType~=3  );
        selectionCorrect = ( selectionChosen  & selectionUnrepeated & data.errorType~=4  );

    % data completion
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

        data.deltaOrdinalReward = data.ordinalRewardRight - data.ordinalRewardLeft;
        data.deltaCardinalReward = data.cardinalRewardRight - data.cardinalRewardLeft;
        
        data.deltaOrdinalEffort = data.ordinalEffortRight - data.ordinalEffortLeft;
        data.deltaCardinalEffort = data.cardinalEffortRight - data.cardinalEffortLeft;

        
%% descriptive
    % define output
        result.summary.treatment = data.treatment(1);
        result.summary = struct2dataset(result.summary);
        taskClock = clock;
        taskTime = [ num2str(taskClock(3)) '_' num2str(taskClock(2)) '_' num2str(taskClock(1)) ];
        result.summary.Properties.ObsNames = taskTime;
    
    % univariate
        result.summary.nParticipation = sum(selectionUnmissed); % count omissions
        result.summary.omissionRate = sum(selectionUnmissed==0)/numel(selectionUnmissed);
        result.summary.prematureRate = sum(data.errorType(selectionUnmissed)==2 )/numel(find(selectionUnmissed));
        result.summary.sideErrorRate = sum(data.errorType(selectionUnmissed)==3 )/numel(find(selectionUnmissed));
        result.summary.executionErrorRate = sum(data.errorType(selectionChosen)==4 )/numel(find(selectionChosen));
        result.summary.correctRate = sum(selectionCorrect)/numel(selectionCorrect);


    % bivariate
        % errors observed
        selection = selectionChosen;
            y = (data.errorType(selection)==4);
                x = data.sideChoice(selection); 
                result.descriptive.quantile_executionError_sideChoice  = unique(x);
                result.descriptive.mean_executionError_by_sideChoice= tapply( y,{ x },@nanmean, {'discrete'});
                result.descriptive.sem_executionError_by_sideChoice = tapply( y ,{ x },@sem, {'discrete'});
    
                result.summary.participationRate_byEffort =  tools.tapply(selectionUnmissed,{data.ordinalEffortLeft+data.ordinalEffortLeft},@nanmean);
                
        % force observed
%         selection = selectionCorrect;
        selection = selectionChosen;
            y = data.perf(selection);
                x = data.chosenCardinalEffort(selection)*100; 
                result.descriptive.quantile_perf_chosenEffort = unique(x);
                result.descriptive.mean_perf_by_chosenEffort = tapply( y,{ x },@nanmean, {'discrete'});
                result.descriptive.sem_perf_by_chosenEffort = tapply( y ,{ x },@sem, {'discrete'});

                x = data.chosenCardinalReward(selection); 
                result.descriptive.quantile_perf_chosenReward = unique(x);
                result.descriptive.mean_perf_by_chosenReward = tapply( y,{ x },@nanmean, {'discrete'});
                result.descriptive.sem_perf_by_chosenReward = tapply( y ,{ x },@sem, {'discrete'});
                
                x = data.sideChoice(selection); 
                result.descriptive.quantile_perf_sideChoice = unique(x);
                result.descriptive.mean_perf_by_sideChoice= tapply( y,{ x },@nanmean, {'discrete'});
                result.descriptive.sem_perf_by_sideChoice = tapply( y ,{ x },@sem, {'discrete'});
                
             y = data.cumulativePerf(selection);
                x = data.chosenCardinalEffort(selection)*100*(data.effortUnit(1)*forceSamplingFreq); 
                result.descriptive.quantile_cumulativePerf_chosenEffort = unique(x);
                result.descriptive.mean_cumulativePerf_by_chosenEffort = tapply( y,{ x },@nanmean, {'discrete'});
                result.descriptive.sem_cumulativePerf_by_chosenEffort = tapply( y ,{ x },@sem, {'discrete'});

                x = data.chosenCardinalReward(selection); 
                result.descriptive.quantile_cumulativePerf_chosenReward = unique(x);
                result.descriptive.mean_cumulativePerf_by_chosenReward = tapply( y,{ x },@nanmean, {'discrete'});
                result.descriptive.sem_cumulativePerf_by_chosenReward = tapply( y ,{ x },@sem, {'discrete'});

        % choice observed
        selection = selectionChosen;
            y = data.sideChoice(selection)-1;
                x = data.deltaOrdinalReward(selection);
                result.descriptive.quantile_sideChoice_deltaOrdinalReward = unique(x);
                result.descriptive.mean_sideChoice_by_deltaOrdinalReward = tapply( y,{ x },@nanmean, {'discrete'});
                result.descriptive.sem_sideChoice_by_deltaOrdinalReward = tapply( y ,{ x },@sem, {'discrete'});

                x = data.deltaOrdinalEffort(selection); 
                result.descriptive.quantile_sideChoice_deltaOrdinalEffort = unique(x);
                result.descriptive.mean_sideChoice_by_deltaOrdinalEffort = tapply( y,{ x },@nanmean, {'discrete'});
                result.descriptive.sem_sideChoice_by_deltaOrdinalEffort = tapply( y ,{ x },@sem, {'discrete'});
                
            y = y(2:end);
                x = data.sideChoice(selection)-1; x = x(1:end-1);
                result.descriptive.quantile_sideChoice_previousSideChoice = unique(x);
                result.descriptive.mean_sideChoice_by_previousSideChoice = tapply( y,{ x },@nanmean, {'discrete'});
                result.descriptive.sem_sideChoice_by_previousSideChoice = tapply( y ,{ x },@sem, {'discrete'});
                
                
            result.summary.choiceRepetitionRate = nanmean([ 1 - result.descriptive.mean_sideChoice_by_previousSideChoice(1) , result.descriptive.mean_sideChoice_by_previousSideChoice(2) ] );
        


%% inferential
  % force observed
        %         selection = selectionCorrect;
        selection = selectionChosen;
            % fit perf
                y = data.perf(selection);
                x1 = data.chosenCardinalEffort(selection)*100; x2 = data.chosenCardinalReward(selection); x3 = data.sideChoice(selection); 
                    result.inferential.glm_perf = GeneralizedLinearModel.fit([x1,x2,x3],y,...
                        'force ~ chosenEffort + chosenReward + chosenSide + chosenEffort:chosenSide',...
                        'Distribution','normal','CategoricalVars',[3],'VarNames',{'chosenEffort';'chosenReward';'chosenSide';'force'});
            % predictions
                y = result.inferential.glm_perf.Fitted.Response;
                data.predictedPerf = y;
                result.inferential.prediction_perf_by_chosenEffort = tapply( y,{ x1 },@nanmean, {'discrete'});
                
            % fit cumulative perf
                y = data.cumulativePerf(selection);
                x1 = data.chosenCardinalEffort(selection)*100*(data.effortUnit(1)*forceSamplingFreq); x2 = data.chosenCardinalReward(selection); x3 = data.sideChoice(selection); 
                    result.inferential.glm_cumulativePerf = GeneralizedLinearModel.fit([x1,x2,x3],y,...
                        'cumulativeForce ~ chosenEffort + chosenReward + chosenSide + chosenEffort:chosenSide',...
                    'Distribution','normal','CategoricalVars',[3],'VarNames',{'chosenEffort';'chosenReward';'chosenSide';'cumulativeForce'});

    % choice observed
        selection = selectionChosen;
              % fit choice
                 y = data.sideChoice(selection)-1;
                 x1 = data.deltaOrdinalReward(selection); x2 = data.deltaOrdinalEffort(selection); x3 = y(1:end-1); x3 = [NaN; x3];
                    result.inferential.glm_choice = GeneralizedLinearModel.fit([x1,x2,x3],y,...
                        'choice ~ deltaReward + deltaEffort + previousChoice',...
                        'Distribution','binomial','link','logit','CategoricalVars',[3],'VarNames',{'deltaReward';'deltaEffort';'previousChoice';'choice'});
                    
                % predictions
                    y = result.inferential.glm_choice.Fitted.Response;
                    data.predictedChoice= y;
                    result.inferential.prediction_sideChoice_by_deltaOrdinalReward = tapply( y,{ x1 },@nanmean, {'discrete'});
                    result.inferential.prediction_sideChoice_by_deltaOrdinalEffort = tapply( y,{ x2 },@nanmean, {'discrete'});
                    
                result.summary.maxReward_choiceRate = nanmean([ nanmean(y(x1>0)) , nanmean(1-y(x1<0)) ]);
                result.summary.minEffort_choiceRate = nanmean([ nanmean(y(x2<0)) , nanmean(1-y(x2>0)) ]);
                result.summary.choicePredictedByPreviousChoiceRate = mean([ mean(y(x3==1)) , mean(1-y(x3==0)) ]);


%% report

    clc;
    fprintf(' %% ----- résumé comportemental ----- %% \n');
    fprintf('\n');
    disp(result.summary);
    fprintf('force regression\n');
    disp(result.inferential.glm_perf.Coefficients);
    fprintf('choice regression\n');
    disp(result.inferential.glm_choice.Coefficients);


   
end