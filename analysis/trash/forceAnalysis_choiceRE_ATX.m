%% forceAnalysis_choiceRE_ATX
% this script execute the analysis of choice of the choiceRE task
% for monkeys & visualize the results
%
% Nicolas Borderies 05/2016
%
  
post = struct;
data.forceResidual = nan(height(data),1);
data.forcePredicted = nan(height(data),1);


for iSub = 1:numel(subjectList)
         % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen ;
        
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            
            choice = data.sideChoice(selection)-1;   
            force = data.perf(selection);   
            
            sv = data.participationValue(selection);
            r1 = data.ordinalRewardLeft(selection);
            r2 = data.ordinalRewardRight(selection);
            rchoice = data.chosenCardinalReward(selection);
            e1 = data.ordinalEffortLeft(selection);
            e2 = data.ordinalEffortRight(selection);
            echoice = data.chosenCardinalEffort(selection);
            dr = round(r2-r1);
            de = round(e2-e1);
            sumr = round(r2+r1);
            sume = round(e2+e1);
            conflict = dr.*de;
%             conflict = ((dr.*de)>0);
            time = data.offerTime(selection);
            nt = data.trialNumber(selection);
            rt = data.responseTime(selection)+0.1;        

        % fit
             predictor = array2table([zscore(rchoice),zscore(echoice),zscore(sv),zscore(nt),zscore(choice),double(treatment) , force ]);
             predictor.Properties.VariableNames = {'rchoice';'echoice';'offerValue';'trialNumber';'side';'treatment';
                                                    'force'};
%              formula = [' force ~ 1 + side + rchoice*echoice + trialNumber ',...
%                         ' + treatment:( 1 + side + rchoice*echoice + trialNumber)'];
             formula = [' force ~ 1 + side + rchoice*echoice  ',...
                        ' + treatment:( 1 + side + rchoice*echoice )'];
                    %              formula = [' force ~ 1 + side + rchoice*echoice + trialNumber '];
%              formula = [' force ~ 1 + side + echoice*offerValue ',...
%                         ' + treatment:( 1 + side + echoice*offerValue)'];
             predictor = predictor(selectSubject(selection),:);
        % classical mll 
             glm = fitglm( predictor, formula,'Distribution','gamma','link','identity','CategoricalVars',[6]);
             data.forcePredicted(selection) =  glm.Fitted.Response;
             displayGLMcoeff(glm)
             data.forceResidual(selection) = glm.Residuals.Raw;
             
             figure; hold on;
             for i =1:2
                x = data.forceResidual(selection & data.treatment==trtList(i));
                histogram(x,'BinWidth',0.02,'Normalization','Probability','FaceColor',col{i},'EdgeColor',col{i},'FaceAlpha',0.5);
                xlabel('force residuals');
             end
             
         % fit residuals
             error = abs(data.forceResidual(selection));
             predictor = array2table([zscore(rchoice),zscore(echoice),zscore(sv),double(treatment) , error ]);
             predictor.Properties.VariableNames = {'rchoice';'echoice';'offerValue';'treatment';
                                                    'error'};
             formula = [' error ~ 1 + treatment '];
             predictor = predictor(selectSubject(selection),:);
            % classical mll 
             glm = fitglm( predictor, formula,'Distribution','gamma','link','identity','CategoricalVars',[4]);
             displayGLMcoeff(glm)
         
%         % vba inversion
%              paramNames= {'kR','kR2','kE','kE2','kRt','kEt',...
%                     'kR_atx','kRt_atx','kE_atx','kEt_atx','k_atx'};
%              [stat] = VBA_rt_ATX(predictor,[],paramNames);
%             
%             c = [1 1 1 1 1 1  1 0 0 0 0 ;
%                  1 1 1 1 1 1  0 1 0 0 0 ;
%                  1 1 1 1 1 1  0 0 1 0 0 ;
%                  1 1 1 1 1 1  0 0 0 1 0 ;
%                  1 1 1 1 1 1  0 0 0 0 1 ;
%                  1 1 1 1 1 1  0 0 0 0 0 ];
%                  
%             [stat] = VBA_rt_ATX(predictor,[],paramNames,c,stat);
%             result.choice.model.(sub).stat = stat;
            
         
end

%         % averaging across subjects
%         post = cell(2,1);
%         post{1} = result.choice.model.aliosha.stat.posterior;
%         post{2} = result.choice.model.bob.stat.posterior;
%         F = [ result.choice.model.aliosha.stat.logE ,...
%              result.choice.model.bob.stat.logE];
%         
%         [meanPost] = VBA_BMA(post,F);
%         beta = array2table([meanPost.muPhi';diag(meanPost.SigmaPhi)' ; ones(1,numel(meanPost.muPhi)) ],...
%                             'VariableNames',paramNames,...
%                             'RowNames',{'mu','std','p'});         
%          f = displayLogit( beta ) ;
        
        
%%
fig = figure; set(fig,'Name','force_errors');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.forcePredicted(selection); 
%         force = data.perf(selection);   
        effort = data.chosenCardinalEffort(selection);   

        error = abs(data.forceResidual(selection));   

        
  
        nbin = 4;
        x = tapply(force,{effort,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        mu =  tapply(error,{effort,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        err =   tapply(error,{effort,treatment},@sem,{'continuous','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('predicted force'); 
        ylabel('residual force'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
   
%%
fig = figure; set(fig,'Name','force_24');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;

        session = data.session(selection);
        treatment = data.treatment(selection);
%         force = data.forcePredicted(selection); 
        force = data.perf(selection);   
        effort = data.chosenCardinalEffort(selection);   
        nt = data.normalizedTrialNumber(selection);
        error = abs(data.forceResidual(selection));   

        
  
        nbin = [4 10 2];
        x = tapply(nt,{effort,nt,treatment},@nanmean,{'discrete','continuous','discrete'},nbin);
        mu =  tapply(force,{effort,nt,treatment},@nanmean,{'discrete','continuous','discrete'},nbin);
        err =   tapply(force,{effort,nt,treatment},@sem,{'discrete','continuous','discrete'},nbin);
        
        for i = 1:2
            for e =1:4
                 [~,~,h(i)] = errorscat(x(e,:,i) ,mu(e,:,i), err(e,:,i),col{i});
                 h(i).LineWidth = e;
            end
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('session progression (total trial)'); 
        ylabel('force (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
fig = figure; set(fig,'Name','force_25');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
                
%         selection = selectSubject & selectionChosen & data.chosenOrdinalEffort==4 ;
        selection = selectSubject & selectionChosen ;
%         selection = selectSubject  ;

        session = data.session(selection);
        treatment = data.treatment(selection);
%         force = data.forcePredicted(selection); 
        force = data.perf(selection);
        nt = data.trialNumber(selection);
        participation =  selectionChosen(selection);
        outcome = data.ordinalRewardOutcome(selection);
        effort = data.chosenCardinalEffort(selection);   
        nt = data.trialNumber(selection);
        error = abs(data.forceResidual(selection));  
        err2 = force - effort;
        sd = data.stateDuration(selection);
        
        % previous force
            ft = [NaN ; force(1:end-1)];
            ft(nt==1) = NaN ;
        % previous part
            pt = [NaN ; participation(1:end-1)];
            pt(nt==1) = NaN ;
        % previous outcome
            ot = [NaN ; outcome(1:end-1)];
            ot(nt==1) = NaN ;
        % previous err
            errt = [NaN ; force(1:end-1) - effort(1:end-1)];
            errt(nt==1) = NaN ;
            
 

        
  
%         nbin = [2 2];
%         x = tapply(pt,{pt,treatment},@nanmean,{'discrete','discrete'},nbin);
%         mu =  tapply(error,{pt,treatment},@nanmean,{'discrete','discrete'},nbin);
%         err =   tapply(error,{pt,treatment},@sem,{'discrete','discrete'},nbin);
        
        nbin = [6 2];
        x = tapply(sd,{sd,treatment},@nanmean,{'continuous','discrete'},nbin);
        mu =  tapply(error,{sd,treatment},@nanmean,{'continuous','discrete'},nbin);
        err =   tapply(error,{sd,treatment},@sem,{'continuous','discrete'},nbin);
        
        for i = 1:2
             [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});


        
        ax = gca;
        xlabel('force (t-1)'); 
        ylabel('force error (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
        
    end
    