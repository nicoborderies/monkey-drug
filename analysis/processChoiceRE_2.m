function [result,data] = processChoiceRE_2(data)
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
        sessionRT = nominal({'11_18_2015'});
        selectionResponseTime = ~selectionRepetition & ~selectionMissed & ~isnan(data.intertrialTime) & ismember(data.session,sessionRT);
      
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
        
        
%% compute result output
    % define output
        result = struct;
        clc; fprintf('\n');

        
    %% design
        result.design.subject = data.subject(1);
        result.design.experimenter = data.experimenter(1);
        result.design.treatment = data.treatment(1);
        result.design.calibLeft = data.calibLeft(1);
        result.design.calibRight = data.calibRight(1);
        result.design.rewardUnit = data.rewardUnit(1);
        result.design.effortUnit = data.effortUnit(1);
        result.design.nTrial = max(data.trialNumber);
        
        taskClock = clock;
        taskTime = [ num2str(taskClock(3)) '/' num2str(taskClock(2)) '/' num2str(taskClock(1)) ' ' ...
                     num2str(taskClock(4)) 'h' num2str(taskClock(5)) 'min'  ];
%         result.design.Properties.ObsNames = taskTime;
        result.design.taskTime = taskTime;
        
        % display
        fprintf(' Experimental design : \n');
        fprintf('_______________________\n\n');

        disp(result.design);
        result.design = struct2dataset(result.design);


        fprintf('\n  Results : \n');
          fprintf('____________\n\n');

        

    %% Participation Pattern
        result.inferential.participation.nParticipationWithRepetition = sum(selectionParticipate); % count omissions
        result.inferential.participation.nParticipation = sum(selectionParticipate(~selectionRepetition)); % count unrepeated omissions
        result.inferential.participation.nChoice = sum(selectionChosen(~selectionRepetition));  % count choice
        result.inferential.participation.totalReward = [ num2str(round(nansum(data.volumeRewardOutcome))) 'ml' ];
        
        result.inferential.participation.omissionRate = mean(selectionMissed(~selectionRepetition));
        result.inferential.participation.prematureRate =  mean(selectionPremature(~selectionRepetition));
        result.inferential.participation.sideErrorRate = mean(selectionSwitch(~selectionRepetition));
        result.inferential.participation.executionErrorRate = mean(selectionIncorrect(~selectionRepetition & selectionChosen));
        
        % participation analysis
            % select
             selection = ~selectionRepetition;
             y = data.errorType(selection)~=1;         
             x1 = data.maxReward(selection);
             x2 = data.minReward(selection);
             x3 = data.maxEffort(selection);
             x4 = data.minEffort(selection);
             x5 = data.trialNumber(selection);
             x6 = y(1:end-1); x6 = [NaN; x6]; % previousParticipation
             x7 = (data.errorType(selection)==0); x7(end) = []; x7 = [NaN; x7]; % previousCorrect
             x8 = data.ordinalRewardOutcome(selection) ; x8(end) = []; x8 = [NaN; x8]; % previousReward
             x9 = data.cumulativeReward(selection);
             x10 = data.cumulativeEffort(selection);

              % fit
              f = figure; hold on;
              f.Units = 'normalized';
              f.Position = [0.05 0.05 0.6 0.9];
              for iSub = 1:numel(subjectList)
                    sub = subjectList{iSub};
                    selectSubject = ismember(data.subject, sub);

                     predictor = [zscore(x1),zscore(x2),zscore(x3),zscore(x4),zscore(x5),x6,x7,nanzscore(x8)];
                     predictor = [zscore(x1),zscore(x2),zscore(x3),zscore(x4),zscore(x5),x6,x7,nanzscore(x8),nanzscore(x9),nanzscore(x10)];

                     predictor = predictor(selectSubject(selection),:);
                     predicted = y(selectSubject(selection));
                     
%                      formula = 'participation ~ maxR + minR + maxE + minE + nTrial*(maxR + minR + maxE + minE) + previousParticipation + previousCorrect + previousReward';
                     formula = 'participation ~ maxR + minR + maxE + minE + nTrial + cumulativeReward + cumulativeEffort  + previousCorrect';

                     result.model.participation.glm = fitglm( predictor, predicted,...
                        formula,...
                        'Distribution','binomial','link','logit','CategoricalVars',[6,7],...
                        'VarNames',{'maxR';'minR';'maxE';'minE';'nTrial';'previousParticipation';'previousCorrect';'previousReward';'cumulativeReward';'cumulativeEffort';'participation'});
                    glm  = result.model.participation.glm;
                    
                    subplot(numel(subjectList),1,iSub);hold on
                    ax = gca;
                    option.pValue=0;
                    f = displayGLMcoeff( glm , f , ax, option );
                    ax.FontSize = 14;
                    ax.Color = 'none';
                    if iSub~= numel(subjectList)
                        ax.XTick = [];
                        ax.XLabel = [];
                    end
                    ylabel(ax,'% of variance');
                    t = title(ax,['monkey ' sub(1) ]); 
              end
        
              
%         % premature analysis
%             % select
%              selection = ~selectionRepetition & ~selectionMissed;
%              y = data.errorType(selection)==2;         
%              x1 = data.maxReward(selection);
%              x2 = data.minReward(selection);
%              x3 = data.maxEffort(selection);
%              x4 = data.minEffort(selection);
%              x5 = data.trialNumber(selection);
%              x6 = y(1:end-1); x6 = [NaN; x6]; % previousPremature
%              x7 = (data.errorType(selection)==0); x7(end) = []; x7 = [NaN; x7]; % previousCorrect
%              x8 = data.ordinalRewardOutcome(selection) ; x8(end) = []; x8 = [NaN; x8]; % previousReward
% 
%             
%              % fit
%               for iSub = 1:numel(subjectList)
%                     sub = subjectList{iSub};
%                     selectSubject = ismember(data.subject, sub);
% 
%                     predictor = [zscore(x1),zscore(x2),zscore(x3),zscore(x4),zscore(x5),x6,x7,nanzscore(x8)];
%                     predictor = predictor(selectSubject(selection),:);
%                     predicted = y(selectSubject(selection));
%                     
%                      result.model.premature.glm = fitglm( predictor, predicted,...
%                         'premature ~ maxR + minR + maxE + minE + nTrial*(maxR + minR + maxE + minE) + previousPremature + previousCorrect + previousReward',...
%                         'Distribution','binomial','link','logit','CategoricalVars',[6,7],...
%                         'VarNames',{'maxR';'minR';'maxE';'minE';'nTrial';'previousPremature';'previousCorrect';'previousReward';'premature'});
%                     glm  = result.model.premature.glm;
%                     
%                     f = displayGLMcoeff( glm );
%                     ax = get(f,'child');
%                     ylabel(ax,'regression coefficients (%variance)');
%                     t = title(ax,['monkey ' sub(1) ]); 
%               end
         
        
        % display
        fprintf(' Participation pattern : \n\n');
        disp(result.inferential.participation);

   %%  Choice Pattern
        % select
         selection = selectionChosen & ~selectionRepetition ;
         y = data.sideChoice(selection)-1;
         x1 = data.deltaOrdinalReward(selection);
         x2 = data.deltaOrdinalEffort(selection);
         x3 = y(1:end-1); x3 = [NaN; x3]; % previousChoice
         x4 = (data.errorType(selection)==0); x4(end) = []; x4 = [NaN; x4]; % previousCorrect
         x5 = data.trialNumber(selection);
         x6 = data.maxReward(selection);
         x7 = data.minReward(selection);
         x8 = data.maxEffort(selection);
         x9 = data.minEffort(selection);
         
         % describe
         result.inferential.choice.maxRewardChoiceRate = mean([ mean(y(x1>0)) , mean(1-y(x1<0)) ]);
         result.inferential.choice.minEffortChoiceRate = mean([ mean(y(x2<0)) , mean(1-y(x2>0)) ]);
         result.inferential.choice.rightChoiceRate = mean(y);
         result.inferential.choice.repetitionChoiceRate = mean([ mean(y(x3==1)) , mean(1-y(x3==0)) ]);
         result.inferential.choice.correctStayRate = mean([ mean(y(x3==1 & x4==1 )) , mean(1-y(x3==0 & x4==1)) ]);
         result.inferential.choice.incorrectStayRate = mean([ mean(y(x3==1 & x4==0 )) , mean(1-y(x3==0 & x4==0)) ]);

         % fit
         f = figure; hold on;
              f.Units = 'normalized';
              f.Position = [0.05 0.05 0.6 0.9];
         data.predictedSideChoice = nan(size(data,1),1);
        for iSub = 1:numel(subjectList)
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);

             predictor = [zscore(x1),zscore(x2),(x3),(x4),zscore(x5)];
             predictor = [zscore(x1),zscore(x2),(x3),(x4),zscore(x5),zscore(x6),zscore(x7),zscore(x8),zscore(x9)];

             predictor = predictor(selectSubject(selection),:);
             predicted = y(selectSubject(selection));
             
%              formula = 'choice ~ dR + dE + dR:dE + dR:nTrial + dE:nTrial + previousChoice + previousChoice:previousCorrect';
             formula = 'choice ~ dR + dE + dR:dE + dR:minR + dE:minE + dR:nTrial + dE:nTrial + previousChoice*previousCorrect';

             
             result.model.choice.glm = fitglm( predictor, predicted,...
                formula,...
                'Distribution','binomial','link','logit','CategoricalVars',[3,4],...
                'VarNames',{'dR';'dE';'previousChoice';'previousCorrect';'nTrial';'maxR';'minR';'maxE';'minE';'choice'});
            glm  = result.model.choice.glm;
            data.predictedSideChoice(selection & selectSubject) = glm.Fitted.Response;
            result.inferential.choice.predictedChoice = balanced_accuracy(glm.Fitted.Response,y);
            
            subplot(numel(subjectList),1,iSub);hold on
                    ax = gca;
                    option.pValue=0;
                    f = displayGLMcoeff( glm , f , ax, option );
                    ax.FontSize = 14;
                    ax.Color = 'none';
                    if iSub~= numel(subjectList)
                        ax.XTick = [];
                        ax.XLabel = [];
                    end
                    ylabel(ax,'% of variance');
                    t = title(ax,['monkey ' sub(1) ]); 
        end
         
        % display
        fprintf(' Choice pattern : \n\n');
        disp(result.inferential.choice);
        
    %%  Force Pattern
       % select
         selection = selectionChosen & ~selectionRepetition;
         y = data.perf(selection);
         x1 = data.chosenOrdinalEffort(selection);
         x2 = data.chosenOrdinalReward(selection);
         x3 = data.sideChoice(selection); x3 = ((x3-1)*2)-1; % x3 = { 1 / -1}
         x4 = data.trialNumber(selection);
         x5 = y(1:end-1); x5 = [NaN; x5]; % previousForce
         x6 = (data.errorType(selection)==0); x6(end) = []; x6 = [NaN; x6]; x6 = ((x6)*2)-1; % previousCorrect
         x7 = x3(1:end-1); x7 = [NaN; x7]; % previousSide
         x8 = data.cumulativeReward(selection);
         x9 = data.cumulativeEffort(selection);
        
            
            
          % fit
          f = figure; hold on;
              f.Units = 'normalized';
              f.Position = [0.05 0.05 0.6 0.9];
              data.predictedForce = nan(size(data,1),1);
              for iSub = 1:numel(subjectList)
                    sub = subjectList{iSub};
                    selectSubject = ismember(data.subject, sub);

                    predictor = [zscore(x1),zscore(x2),zscore(x3),zscore(x4),nanzscore(x5),nanzscore(x6),nanzscore(x7),nanzscore(x8),nanzscore(x9)];
                    predictor = predictor(selectSubject(selection),:);
                    predicted = y(selectSubject(selection));
                    
%                     formula = 'force ~ chosenEffort + chosenReward + chosenSide + chosenEffort:chosenSide + nTrial*(previousForce) + previousSide:chosenSide + previousCorrect';
                    formula = 'force ~ chosenEffort + chosenReward + chosenReward:chosenEffort + chosenSide  + nTrial + cumulativeEffort + cumulativeReward ';

                    
                    result.model.force.glm = fitglm( predictor, predicted,...
                        formula,...
                        'Distribution','normal','VarNames',{'chosenEffort';'chosenReward';'chosenSide';'nTrial';'previousForce';'previousCorrect';'previousSide';'cumulativeReward';'cumulativeEffort';'force'});
                    glm  = result.model.force.glm;
                    data.predictedForce(selection & selectSubject) = glm.Fitted.Response;
            
                    subplot(numel(subjectList),1,iSub);hold on
                    ax = gca;
                    option.pValue=0;
                    f = displayGLMcoeff( glm , f , ax, option );
                    ax.FontSize = 14;
                    ax.Color = 'none';
                    if iSub~= numel(subjectList)
                        ax.XTick = [];
                        ax.XLabel = [];
                    end
                    ylabel(ax,'% of variance');
                    t = title(ax,['monkey ' sub(1) ]); 
              end
      
        % display
        fprintf(' Force pattern : \n\n');
%         disp(result.inferential.force);

%% RT analysis
            % select
             selection =  selectionResponseTime;
             y = data.reactionTime(selection);         
             x1 = data.maxReward(selection);
             x2 = data.minReward(selection);
             x3 = data.maxEffort(selection);
             x4 = data.minEffort(selection);
             x5 = data.trialNumber(selection);
             x6 = data.offerTime(selection);
             x7 = data.intertrialTime(selection);

              % fit
              f = figure; hold on;
              f.Units = 'normalized';
              f.Position = [0.05 0.05 0.6 0.9];
              for iSub = 1:numel(subjectList)
                    sub = subjectList{iSub};
                    selectSubject = ismember(data.subject, sub);

                     predictor = [zscore(x1),zscore(x2),zscore(x3),zscore(x4),zscore(x5),nanzscore(x6),nanzscore(x7)];

                     predictor = predictor(selectSubject(selection),:);
                     predicted = y(selectSubject(selection));
                     
%                      formula = 'participation ~ maxR + minR + maxE + minE + nTrial*(maxR + minR + maxE + minE) + previousParticipation + previousCorrect + previousReward';
                     formula = 'responseTime ~ maxR + minR + maxE + minE + nTrial + offerTime + intertrialTime';

                     result.model.responseTime.glm = fitglm( predictor, predicted,...
                        formula,...
                        'Distribution','gamma','link','log',...
                        'VarNames',{'maxR';'minR';'maxE';'minE';'nTrial';'offerTime';'intertrialTime';'responseTime'});
                    glm  = result.model.responseTime.glm;
                    
                    subplot(numel(subjectList),1,iSub);hold on
                    ax = gca;
                    option.pValue=0;
                    f = displayGLMcoeff( glm , f , ax, option );
                    ax.FontSize = 14;
                    ax.Color = 'none';
                    if iSub~= numel(subjectList)
                        ax.XTick = [];
                        ax.XLabel = [];
                    end
                    ylabel(ax,'% of variance');
                    t = title(ax,['monkey ' sub(1) ]); 
              end
   
end