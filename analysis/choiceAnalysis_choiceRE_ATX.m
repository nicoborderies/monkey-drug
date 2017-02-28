%% choiceAnalysis_choiceRE_ATX
% this script execute the analysis of choice of the choiceRE task
% for monkeys & visualize the results
%
% Nicolas Borderies 05/2016
%
  
post = struct;
data.choiceValue = nan(height(data),1);
data.choicePredicted = nan(height(data),1);

for iSub = 1:numel(subjectList)
         % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen ;
        
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            choice = data.sideChoice(selection)-1;   
            x1 = data.ordinalRewardLeft(selection);
            x2 = data.ordinalRewardRight(selection);
            x3 = data.ordinalEffortLeft(selection);
            x4 = data.ordinalEffortRight(selection);
            nt = data.trialNumber(selection);
            nnt = data.normalizedTrialNumber(selection);
            ee =  data.cumulativeEffort(selection);
            rr =  data.cumulativeReward(selection);
        
        % fit
             predictor = array2table([rangescore(x1),rangescore(x2),rangescore(x3),rangescore(x4),rangescore(nt),...
                                double(treatment) ,rangescore(rr),rangescore(ee), choice]);
             predictor.Properties.VariableNames = {'R1';'R2';'E1';'E2';'ntrial';'treatment';
                                                    'Rtot';'Etot';
                                                    'choice'};
%              formula = [' choice ~ -1 + R1 + R2 + E1 + E2 + ntrial:(R1 + R2 + E1 + E2) '];
%              formula = [' choice ~ -1 + E1 + E2 + ntrial:(E1 + E2) + treatment:( E1 + E2 ) '];

             formula = [' choice ~ -1 + R1 + R2 + E1 + E2 + Rtot:(R1 + R2)  + Etot:(E1 + E2) ',...
                        ' + treatment:( R1 + R2 + E1 + E2 + Rtot:(R1 + R2)  + Etot:(E1 + E2))'];
             predictor = predictor(selectSubject(selection),:);
%         % classical mll 
             glm = fitglm( predictor, formula,'Distribution','binomial','link','logit','CategoricalVars',[6]);
             y2  = glm.Fitted.Response;
             displayGLMcoeff(glm);
             data.choiceValue(selection) = glm.Fitted.LinearPredictor ;
             data.choicePredicted(selection) = y2 ;
%              
%         % vba inversion
%              paramNames= {'kR','kR2','kE','kE2','kRt','kEt',...
%                     'kR_atx','kRt_atx','kE_atx','kEt_atx','k_atx'};
%              [stat] = VBA_choice_ATX(predictor,[],paramNames);
%             
%             c = [1 1 1 1 1 1  1 1 1 1 1 ;
%                 
%                  1 1 1 1 1 1  1 0 0 0 0 ;
%                  1 1 1 1 1 1  0 1 0 0 0 ;
%                  1 1 1 1 1 1  0 0 1 0 0 ;
%                  1 1 1 1 1 1  0 0 0 1 0 ;
%                  1 1 1 1 1 1  0 0 0 0 1 ;
%                  
%                  1 1 1 1 1 1  0 0 0 0 0 ];
%                  
%             [stat] = VBA_choice_ATX(predictor,[],paramNames,c,stat);
%             result.choice.model.(sub).stat = stat;
            
         
end

%         % averaging across subjects
%         post = cell(2,1);
%         post{1} = result.choice.model.aliosha.stat.posterior;
%         post{2} = result.choice.model.bob.stat.posterior;
%         F = [ result.choice.model.aliosha.stat.logE ,...
%              result.choice.model.bob.stat.logE];
% %          
%         F = max(result.choice.model.aliosha.stat.contrast.logE);
%         beta = result.choice.model.aliosha.stat.beta;
%         post{1} = result.choice.model.aliosha.stat.posterior;
%         post{1}.muPhi = beta{1,:}';
%         post{1}.SigmaPhi =  repmat(beta{2,:},numel(post{1}.muPhi),1)'.*eye(numel(post{1}.muPhi));
% 
%         beta = result.choice.model.bob.stat.beta;
%         post{2} = result.choice.model.bob.stat.posterior;
%         post{2}.muPhi = beta{1,:}';
%         post{2}.SigmaPhi =  repmat(beta{2,:},numel(post{2}.muPhi),1)'.*eye(numel(post{2}.muPhi));
%          F = [F ,max(result.choice.model.bob.stat.contrast.logE)];
% 
%          [meanPost] = VBA_BMA(post,F);
%         beta = array2table([meanPost.muPhi';diag(meanPost.SigmaPhi)' ; ones(1,numel(meanPost.muPhi)) ],...
%                             'VariableNames',paramNames,...
%                             'RowNames',{'mu','std','p'});         
%                         
%            indpos = [1 2 4 5 6];
%            indneg = [3];
%           for ind=1:numel(beta{1,:})
%             if ismember(ind,indpos)
%                 beta{1,ind} = safepos(beta{1,ind});
%             elseif ismember(ind,indneg)
%                 beta{1,ind} = -safepos(-beta{1,ind});
%             end
%         end                
%         
%                         
%          f = displayLogit( beta ) ;
        
        
%%
fig = figure; set(fig,'Name','choice20');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen  ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        choice = data.sideChoice(selection)-1;
        accuracy = (abs(choice - data.choicePredicted(selection))<0.5);
        stateTime = data.stateDuration(selection);
  
        nbin = 5;
        x = tapply(stateTime,{stateTime,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        mu =  tapply(accuracy,{stateTime,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        err =   tapply(accuracy,{stateTime,treatment},@sem,{'continuous','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('state duration'); 
        ylabel('accuracy of choice'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    %%
fig = figure; set(fig,'Name','choice21');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen  ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        choice = data.sideChoice(selection)-1;
        accuracy = (abs(choice - data.choicePredicted(selection))<0.5);
        rt = data.responseTime(selection);

        stateTime = data.stateDuration(selection);
        % previousParticipation
        tmax=1;
        x5 = nan(numel(selectionParticipate),tmax);
        for it=1:tmax
            x5(1+it:end,it) = selectionParticipate(1:end-it);
            x5(data.trialNumber<=it,it) = NaN;
        end
        pt1 = x5(selection);
        
        nbin = 2;
        x = tapply(pt1,{pt1,treatment},@nanmean,{'discrete','discrete'},[nbin 2]);
        mu =  tapply(rt,{pt1,treatment},@nanmean,{'discrete','discrete'},[nbin 2]);
        err =   tapply(rt,{pt1,treatment},@sem,{'discrete','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('participation(t-1)'); 
        ylabel('accuracy of choice'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
