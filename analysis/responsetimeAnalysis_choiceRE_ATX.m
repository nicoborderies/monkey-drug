%% responsetimeAnalysis_choiceRE_ATX
% this script execute the analysis of choice of the choiceRE task
% for monkeys & visualize the results
%
% Nicolas Borderies 05/2016
%
  
post = struct;
data.responseTimePredicted = nan(height(data),1);
data.responseTimeResidual = nan(height(data),1);

for iSub = 1:numel(subjectList)
         % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
%             selection = selectSubject & selectionChosen ;
            selection =   selectionChosen ;

        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            choice = data.sideChoice(selection)-1;   
            sv = data.participationValue(selection);
            r1 = data.ordinalRewardLeft(selection);
            r2 = data.ordinalRewardRight(selection);
            e1 = data.ordinalEffortLeft(selection);
            e2 = data.ordinalEffortRight(selection);
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
             predictor = array2table([rangescore(sumr),rangescore(sume),rangescore(sv),rangescore(conflict),...
                                      rangescore(time),double(treatment) ,rangescore(nt), rt]);
%              predictor = array2table([zscore(sumr),zscore(sume),zscore(sv),zscore(conflict),zscore(time),double(treatment) , rt]);
             predictor.Properties.VariableNames = {'sumR';'sumE';'offerValue';'offerConflict';'offerTime';
                                    'treatment';'trialNumber';
                                                    'responseTime'};
%              formula = [' responseTime ~ 1 + sumR + sumE + offerTime ',...
%                         ' + treatment:( 1 + sumR + sumE + offerTime )'];
%              formula = [' responseTime ~ 1 + offerValue + offerTime + trialNumber ',...
%                         ' + treatment:( 1 + offerValue + offerTime + trialNumber )'];
             formula = [' responseTime ~ 1 + offerValue  + trialNumber ',...
                    ' + treatment:( 1 + offerValue  + trialNumber )'];
        % classical mll 
             glm = fitglm( predictor, formula,'Distribution','gamma','link','log','CategoricalVars',[6]);
             y2  = glm.Fitted.Response;
             displayGLMcoeff(glm)
             data.responseTimePredicted(selection) =  glm.Fitted.Response;
             displayGLMcoeff(glm)
             data.responseTimeResidual(selection) = glm.Residuals.Raw;
             figure; hold on;
             for i =1:2
                x = data.responseTimeResidual(selection & data.treatment==trtList(i));
                histogram(x,'BinWidth',0.02,'Normalization','Probability','FaceColor',col{i},'EdgeColor',col{i},'FaceAlpha',0.5);
                xlabel('rt residuals');
             end
         % fit residuals
             error = abs(data.responseTimeResidual(selection));
             predictor = array2table([rangescore(sumr),rangescore(sume),exp(-rangescore(sv)),rangescore(conflict),...
                                      rangescore(time),double(treatment) ,rangescore(nt) , error ]);
             predictor.Properties.VariableNames = {'sumR';'sumE';'offerValue';'offerConflict';'offerTime';
                                                   'treatment';'trialNumber';
                                                   'error'};
             formula = [' error ~ 1 + offerValue + treatment:(1+offerValue) '];
            % classical mll 
             glm = fitglm( predictor, formula,'Distribution','gamma','link','identity','CategoricalVars',[6]);
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
        
        