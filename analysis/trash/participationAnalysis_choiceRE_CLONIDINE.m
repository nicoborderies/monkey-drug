%% participationAnalysis_choiceRE_CLONIDINE
% this script execute the analysis of participation of the choiceRE task
% for monkeys & visualize the results
%
% Nicolas Borderies 05/2016
%
%% Mean Effect
    fig = figure; set(fig,'Name','participation');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject ;
       
        % variables
            errors = data.errorType(selection);
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
        
        % compute stats
            mu =  tapply(participation,{session,treatment},@nanmean);
            [h,p] = ttest2(exnan(mu(:,1)),exnan(mu(:,2)));
            err =   sem(mu,1);
            mu  =   nanmean(mu,1);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [ b(i),~ ] = barplot( 1- 0.25 + 0.5*(i-1) ,mu(i), err(i) , col{i} );
              b(i).BarWidth = 0.3;
            end
            s = sigstar({[0.75,1.25]},p);
            
            
            
        % legending
            legend([b(1) b(2)],{'placebo','clonidine'});
            ax = gca; 
            ax.YLim = [0.5 1];
            ax.TickLength = [0 0];
            ax.XTick = [];

            ylabel('participation (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
    

%% 1. Effect of option value

    fig1 = figure; set(fig1,'Name','participation2');
    fig2 = figure; set(fig2,'Name','participation3');

    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject ;
       
        % variables
            errors = data.errorType(selection);
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
            sumE = round(data.ordinalEffortLeft(selection) + data.ordinalEffortRight(selection));
            sumR = round(data.ordinalRewardLeft(selection) + data.ordinalRewardRight(selection));

            
        % compute stats
            x = tapply(sumE,{sumE,treatment},@nanmean);
            y =  tapply(participation,{sumE,treatment},@nanmean);
            z =   tapply(participation,{sumE,treatment},@sem);

        % display
            figure(fig1);
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [~,~,h(i)] = errorscat(x(:,i) ,y(:,i), z(:,i),col{i});
            end
            
        % legending
            legend([h(1) h(2)],{'placebo','clonidine'});
            ax = gca; 
            ax.YLim = [0 1];
            ax.XLim = [1  9];
            xlabel(texlabel('E_left + E_right'));
            ylabel('participation (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]);
            
            
        % compute stats
            x = tapply(sumR,{sumR,treatment},@nanmean);
            y =  tapply(participation,{sumR,treatment},@nanmean);
            z =   tapply(participation,{sumR,treatment},@sem);

        % display
            figure(fig2);
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [~,~,h(i)] = errorscat(x(:,i) ,y(:,i), z(:,i),col{i});
            end
            
        % legending
            legend([h(1) h(2)],{'placebo','clonidine'});
            ax = gca; 
            ax.YLim = [0 1];
            ax.XLim = [1  9];
            xlabel(texlabel('R_left + R_right'));
            ylabel('participation (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
%%
fig1 = figure; set(fig1,'Name','participation4');

    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject  ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
            participation = selectionParticipate(selection);

            time = data.intertrialTime(selection) ;
%             time = data.offerTime(selection) ;
            sd = data.stateDuration(selection);
            predicted = data.participationPredicted(selection);
            % previousParticipation
            tmax=1;
            x5 = nan(numel(participation),tmax);
            for it=1:tmax
                x5(1+it:end,it) = participation(1:end-it);
                x5(nt<=it,it) = NaN;
            end
            stay = (participation==x5);
            
        % compute stats
            nbin = [10 2];
            x = tapply(sd,{sd,treatment},@nanmean,{'continuous','discrete'},nbin);
            y =  tapply(participation,{sd,treatment},@nanmean,{'continuous','discrete'},nbin);
            y2 =  tapply(predicted,{sd,treatment},@nanmean,{'continuous','discrete'},nbin);
            z =   tapply(participation,{sd,treatment},@sem,{'continuous','discrete'},nbin);
            
%             y =  tapply(stay,{sd,treatment},@nanmean,{'continuous','discrete'},nbin);
%             z =   tapply(stay,{sd,treatment},@sem,{'continuous','discrete'},nbin);
        % display
            figure(fig1);
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [~,~,h(i)] = errorscat(x(:,i) ,y(:,i), z(:,i),col{i}); 
              h(i).LineStyle = 'none';
              h(i) = plot(x(:,i) ,y2(:,i),'Color',col{i});
            end
            
        % legending
            legend([h(1) h(2)],{'placebo','atx'});
            ax = gca; 
            ax.YLim = [0 1];
%             xlabel(texlabel('inter-trial time (s.)'));
            xlabel(texlabel('state duration (number of trial)'));
            ylabel('participation (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]);
            
            
            
    end

    
%
%
%% 2. Effect of dynamics
% bin size for sliding statistics
    nbin = 100;
    
    
%% 2.1. participation  as a function of trial number


%% 2.2. participation offset as a function of trial number
% 2.3. participation reward sensitivity as a function of trial number
% 2.4. participation effort sensitivity as a function of trial number
% 2.5. participation option value sensitivity as a function of trial number
%

f1 = figure; set(f1,'Name','participation_2_2');
f2 = figure; set(f2,'Name','participation_2_3');
f3 = figure; set(f3,'Name','participation_2_4');
f4 = figure; set(f4,'Name','participation_2_5');

for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject ;
        
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
            sumR = (data.cardinalRewardLeft(selection) + data.cardinalRewardRight(selection));
            sumR = round(sumR,3);
            sumE = (data.cardinalEffortLeft(selection) + data.cardinalEffortRight(selection));
            sumE = round(sumE,3);
            nt = data.trialNumber(selection);
            
        % init
            X = nan(numel(sessionList),max(nt));
            Y = nan(numel(sessionList),max(nt));
            Y1 = nan(numel(sessionList),max(nt));
            Y2 = nan(numel(sessionList),max(nt));
            Y3 = nan(numel(sessionList),max(nt));
            Y4 = nan(numel(sessionList),max(nt));

        for is = 1:numel(sessionList)
           % select
               x = nt(session==sessionList(is));
               x2 = sumR(session==sessionList(is));
               x3 = sumE(session==sessionList(is));
               y = participation(session==sessionList(is));
           % sliding fit
                y1 = zeros(1,max(x));
                y2 = zeros(1,max(x));
                y3 = zeros(1,max(x));
                y4 = zeros(1,max(x));
               for j = [1:max(x)]
                     ind = [j-nbin/2-1:j+nbin/2]; ind = ind(ind>0); ind = ind(ind<=max(x));

                    [beta,~,stat] = glmfit([rangescore(x2(ind)),rangescore(x3(ind))],y(ind),'binomial','link','logit');
                    y1(j) = stat.t(1);
                    y2(j) = stat.t(2);
                    y3(j) = stat.t(3);
                    y4(j) = stat.t(2)-stat.t(3);
                    if isnan(y1(j)); y1(j)=0;end
                    if isnan(y2(j)); y2(j)=0;end
                    if isnan(y3(j)); y3(j)=0;end
                    if isnan(y4(j)); y4(j)=0;end

               end
           
           c = col{ ismember(is,trtSessions{2}) + 1 };
       % display individual sessions
%            scatter(x,y1,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%            plot(x,y1,'LineWidth',2,'Color',c);
%            pause;
%         normalization of session length
            ind = round(x*max(nt)/max(x));
           X(is,ind) = x./max(x);
        % without normalization of session length
%            ind = round(x);
%            X(is,ind) = x;
           [ y1 ] = interpslide( y1,'foreward' );
           [ y2 ] = interpslide( y2,'foreward' );
           [ y3 ] = interpslide( y3,'foreward' );
           [ y4 ] = interpslide( y4,'foreward' );
           Y1(is,ind) = y1;
           Y2(is,ind) = y2;
           Y3(is,ind) = y3;
           Y4(is,ind) = y4;
           

           
           
        end
        
        % plots
        figure(f1);hold on
                subplot(1,numel(subjectList),iSub);hold on
                for i =1:2
                    ind = trtSessions{i};

                    % stats
                        x = nanmean(X(ind,:),1);
                        err = sem(Y1(ind,:),1);
                        mu = nanmean(Y1(ind,:),1);

                    [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 

                    cl = col{i};
                    set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
                    h(i).LineStyle = '-'; 

                end
            % legending
                legend([h(1) h(2)],{'placebo','atx'});
                ax = gca;
                xlabel('session progression (total trial)');
                ylabel('offset to option value for participation (t-stat)'); 
                yy =ylim;
                t = title(['monkey ' sub(1) ]);    
                
        figure(f2);hold on
                subplot(1,numel(subjectList),iSub);hold on
                for i =1:2
                    ind = trtSessions{i};

                    % stats
                        x = nanmean(X(ind,:),1);
                        err = sem(Y2(ind,:),1);
                        mu = nanmean(Y2(ind,:),1);

                    [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 

                    cl = col{i};
                    set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
                    h(i).LineStyle = '-'; 

                end
            % legending
                legend([h(1) h(2)],{'placebo','atx'});
                ax = gca;
                xlabel('session progression (total trial)');
                ylabel('sensitivity to reward for participation (t-stat)'); 
                yy =ylim;
                t = title(['monkey ' sub(1) ]);     
                
         figure(f3);hold on
                subplot(1,numel(subjectList),iSub);hold on
                for i =1:2
                    ind = trtSessions{i};

                    % stats
                        x = nanmean(X(ind,:),1);
                        err = sem(Y3(ind,:),1);
                        mu = nanmean(Y3(ind,:),1);

                    [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 

                    cl = col{i};
                    set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
                    h(i).LineStyle = '-'; 

                end
            % legending
                legend([h(1) h(2)],{'placebo','atx'});
                ax = gca;
                xlabel('session progression (total trial)');
                ylabel('sensitivity to effort for participation (t-stat)'); 
                yy =ylim;
                t = title(['monkey ' sub(1) ]);     
                
         figure(f4);hold on
                subplot(1,numel(subjectList),iSub);hold on
                for i =1:2
                    ind = trtSessions{i};

                    % stats
                        x = nanmean(X(ind,:),1);
                        err = sem(Y4(ind,:),1);
                        mu = nanmean(Y4(ind,:),1);

                    [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 

                    cl = col{i};
                    set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
                    h(i).LineStyle = '-'; 

                end
            % legending
                legend([h(1) h(2)],{'placebo','atx'});
                ax = gca;
                xlabel('session progression (total trial)');
                ylabel('sensitivity to option value for participation (t-stat)'); 
                yy =ylim;
                t = title(['monkey ' sub(1) ]);     
                
end

%% 3. model selection
data.participationValue = nan(height(data),1);
data.participationPredicted = nan(height(data),1);

for iSub = 1:numel(subjectList)
         % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject ;
        
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
            x1 = data.ordinalRewardLeft(selection);
            x2 = data.ordinalRewardRight(selection);
            x3 = data.ordinalEffortLeft(selection);
            x4 = data.ordinalEffortRight(selection);
            nt = data.trialNumber(selection);
            nnt = data.normalizedTrialNumber(selection);
            ee =  data.cumulativeEffort(selection);
            rr =  data.cumulativeReward(selection);

            % previousParticipation
                tmax=10;
                x5 = nan(numel(participation),tmax);
                for it=1:10
                    x5(1+it:end,it) = participation(1:end-it);
                    x5(nt<=it,it) = NaN;
                end
        
        % fit
             predictor = array2table([rangescore(x1),rangescore(x2),rangescore(x3),rangescore(x4),rangescore(nnt),double(treatment) ,(x5*2 -1),rangescore(rr),rangescore(ee), participation]);
             predictor.Properties.VariableNames = {'R1';'R2';'E1';'E2';'ntrial';'treatment';
                                                   'part_t_1';'part_t_2';'part_t_3';'part_t_4';'part_t_5';'part_t_6';'part_t_7';'part_t_8';'part_t_9';'part_t_10';
                                                   'Rtot';'Etot';
                                                   'participation'};
%                  formula = [' participation ~ 1 + R1 + R2 + E1 + E2 ',...
%                             ' + part_t_1 + part_t_2 + part_t_3 + part_t_4 + part_t_5 + part_t_6 + part_t_7 + part_t_8 + part_t_9 + part_t_10',...
%                             ' + part_t_1:( part_t_2 + part_t_3 + part_t_4 + part_t_5 + part_t_6 + part_t_7 + part_t_8 + part_t_9 + part_t_10)',...
%                             ' + treatment:( 1 + R1 + R2 + E1 + E2 ',...
%                             ' + part_t_1 + part_t_2 + part_t_3 + part_t_4 + part_t_5 + part_t_6 + part_t_7 + part_t_8 + part_t_9 + part_t_10',...
%                             ' + part_t_1:( part_t_2 + part_t_3 + part_t_4 + part_t_5 + part_t_6 + part_t_7 + part_t_8 + part_t_9 + part_t_10))'];
                formula = [' participation ~ 1 + R1 + R2 + E1 + E2 + ntrial ',...
                           ' + treatment:(1 + R1 + R2 + E1 + E2 + ntrial) '];
%                 formula = [' participation ~ 1 + R1 + R2 + E1 + E2 + Rtot:(R1 + R2) + Etot:(E1 + E2) ',...
%                            ' + treatment:( R1 + R2 + E1 + E2 + Rtot:(R1 + R2) + Etot:(E1 + E2)) '];
%                 formula = [' participation ~ 1 + R1 + R2 + E1 + E2 + Rtot:(R1 + R2) + Etot:(E1 + E2) '];
                       predictor = predictor(selectSubject(selection),:);
        % classical mll fit
             glm = fitglm( predictor, formula,'Distribution','binomial','link','logit','CategoricalVars',[6]);
             data.participationValue(selection) = glm.Fitted.LinearPredictor ;
             y2  = glm.Fitted.Response;
             data.participationPredicted(selection) = y2 ;
              displayGLMcoeff(glm);


             
%         % vba fit
%             paramNames = {'kr1','kr2','ke1','ke2',...
%                           'p0','pt_1','kt','kf',...
%                           'kr1_atx','kr2_atx','ke1_atx','ke2_atx',...
%                           'p0_atx','pt_1_atx','kt_atx','kf_atx',...
%                           'krt','ket','krt_atx','ket_atx'};
%             [stat] = VBA_participation_ATX(predictor,[],paramNames);
%         
%              c = [1 1 1 1  1 1 1 1  1 0 1 0  1 1 1 1  1 1 1 1;
%                  
%                   1 1 1 1  1 1 1 1  1 0 1 0  0 0 0 0  1 1 1 1;
%                   1 1 1 1  1 1 1 1  1 0 1 0  1 0 0 0  1 1 0 0;
%                   1 1 1 1  1 1 1 1  0 0 0 0  0 0 0 1  1 1 1 1;
%                   1 1 1 1  1 1 1 1  0 0 0 0  1 1 1 0  1 1 0 0;
%                   
%                   1 1 1 1  1 1 1 1  1 0 0 0  0 0 0 0  1 1 0 0;
%                   1 1 1 1  1 1 1 1  0 0 1 0  0 0 0 0  1 1 0 0;
%                   1 1 1 1  1 1 1 1  0 0 0 0  1 0 0 0  1 1 0 0;
%                   1 1 1 1  1 1 1 1  0 0 0 0  0 1 0 0  1 1 0 0;
%                   1 1 1 1  1 1 1 1  0 0 0 0  0 0 1 0  1 1 0 0;
%                   1 1 1 1  1 1 1 1  0 0 0 0  0 0 0 1  1 1 0 0;
%                   1 1 1 1  1 1 1 1  0 0 0 0  0 0 0 0  1 1 1 0;
%                   1 1 1 1  1 1 1 1  0 0 0 0  0 0 0 0  1 1 0 1;
%                   
%                   1 1 1 1  1 1 1 1  0 0 0 0  0 0 0 0  1 1 0 0];
%               
%             [stat] = VBA_participation_ATX(predictor,[],paramNames,c,stat);
%             result.participation.model.(sub).stat = stat;
end
%         % averaging across subjects
%         post = cell(2,1);
%         post{1} = result.participation.model.aliosha.stat.posterior;
%         post{2} = result.participation.model.bob.stat.posterior;
%         F = [ result.participation.model.aliosha.stat.logE ,...
%              result.participation.model.bob.stat.logE];
%         
%         [meanPost] = VBA_BMA(post,F);
%         beta = array2table([meanPost.muPhi';diag(meanPost.SigmaPhi)' ; ones(1,numel(meanPost.muPhi)) ],...
%                             'VariableNames',paramNames,...
%                             'RowNames',{'mu','std','p'});         
%                         
%         indpos = [1 2 5 6 7   9 10 11 12   13 16 17 18 19 20];
%         indneg = [3 4 8];
%         for ind=1:numel(beta{1,:})
%             if ismember(ind,indpos)
%                 beta{1,ind} = safepos(beta{1,ind});
%             elseif ismember(ind,indneg)
%                 beta{1,ind} = -safepos(-beta{1,ind});
%             end
%         end                
%         
%          f = displayLogit( beta ) ;
%          
%          % model selection
%           F = [ result.participation.model.aliosha.stat.contrast.logE  ,...
%                 result.participation.model.bob.stat.contrast.logE ];
%          [ms_post,out] = VBA_groupBMC(F);
% %          modelList = {'full model',...
%                        ''};
%             
         

%% 3.history effects

for iSub = 1:numel(subjectList)
         % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject ;
        
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
            rewardOutcome = data.ordinalRewardOutcome(selection);
            sucess = selectionCorrect(selection);
            force =  data.perf(selection);
            residual = participation - data.participationPredicted(selection);
            
            % previousParticipation
                tmax=10;
                part_t = nan(numel(participation),tmax);
                reward_t = nan(numel(participation),tmax);
                success_t = nan(numel(participation),tmax);
                force_t = nan(numel(participation),tmax);
                for it=1:tmax
                    part_t(1+it:end,it) = participation(1:end-it);
                    part_t(nt<=it,it) = NaN;
                    reward_t(1+it:end,it) = rewardOutcome(1:end-it);
                    reward_t(nt<=it,it) = NaN;
                    success_t(1+it:end,it) = sucess(1:end-it);
                    success_t(nt<=it,it) = NaN;
                    force_t(1+it:end,it) = force(1:end-it);
                    force_t(nt<=it,it) = NaN;
                end
                
            % p
            y = nan(1+tmax*3,2);
            z = nan(1+tmax*3,2);

        % fit
%          predictor = array2table(nanzscore([part_t, reward_t , success_t , force_t , residual]));
%          predictor = array2table(nanzscore([part_t, reward_t  , force_t , residual]));
         predictor = array2table(rangescore([part_t, reward_t  , force_t , residual],1));

%              predictor.Properties.VariableNames = {...; 'participation'};
        formula = 'linear' ;
            isIntercept = 0;
            y = nan(isIntercept+tmax*3,2);
            z = nan(isIntercept+tmax*3,2);
                     
        for it = 1:2
           predictor = array2table(nanzscore([part_t, reward_t  , force_t , residual]));
           predictor = predictor(treatment==trtList(it),:);

            
            % classical mll fit
             glm = fitglm( predictor, formula,'Distribution','normal','link','identity','Intercept',logical(isIntercept));
             y2  = glm.Fitted.Response;
%               displayGLMcoeff(glm);
            y(:,it) =  glm.Coefficients.Estimate;
            z(:,it) =  glm.Coefficients.SE;
        end

        
        % plots
         f = figure; set(f,'Name',['participation_history_' sub ]);
         
            for i =1:2

                % stats
                if  isIntercept
                subplot(1,4,1);hold on
                    ind = 1;
                   [ b(i),~ ] = barplot( 1- 0.25 + 0.5*(i-1) , y(ind,i), z(ind,i) , col{i} );
                     b(i).BarWidth = 0.3;
                    xlabel('residual variance');
                    ylabel('influence on participation '); 
                    ax=gca; ax.YLim = [-0.5 0.5];
                end
                
                % stats
                subplot(1,4,2);hold on
                    ind = isIntercept + [1:10];
                    x = [1:numel(ind)];
                    plot(x,zeros(1,numel(ind)),'k--');
                    [h(i),hp] = boundedline( x , y(ind,i) , z(ind,i) , 'alpha','transparency', 0.5 ); 
                    cl = col{i};
                    set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
                    h(i).LineStyle = '-'; 
                    xlabel('participation(t-x)');
                    ylabel('influence on participation '); 
                    ax=gca; ax.YLim = [-0.5 0.5];

                    

                % stats
                subplot(1,4,3);hold on
                    ind = isIntercept + [11:20];
                    x = [1:numel(ind)];
                    [h(i),hp] = boundedline( x , y(ind,i) , z(ind,i) , 'alpha','transparency', 0.5 ); 
                    cl = col{i};
                    set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
                    h(i).LineStyle = '-'; 
                    xlabel('reward(t-x)');
                    ax=gca; ax.YLim = [-0.5 0.5];
                    

                % stats
                subplot(1,4,4);hold on
                    ind = isIntercept + [21:30];
                    x = [1:numel(ind)];
                    [h(i),hp] = boundedline( x , y(ind,i) , z(ind,i) , 'alpha','transparency', 0.5 ); 
                    cl = col{i};
                    set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
                    h(i).LineStyle = '-'; 
                    xlabel('force(t-x)');
                    ax=gca; ax.YLim = [-0.5 0.5];
                    
            end
            
            % legending
                legend([h(1) h(2)],{'placebo','atx'});
                ax = gca;
                yy =ylim;
                t = title(['monkey ' sub(1) ]); 
end



%% 3. Effect of internal value
%
%%
fig = figure; set(fig,'Name','participation20');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject  ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        % previousParticipation
        tmax=1;
        x5 = nan(numel(participation),tmax);
        for it=1:tmax
            x5(1+it:end,it) = participation(1:end-it);
            x5(nt<=it,it) = NaN;
        end
        stateTime = data.stateDuration(selection);
        stay = (participation==x5);
  
        nbin = 5;
        x = tapply(stateTime,{stateTime,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        mu =  tapply(stay,{stateTime,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        err =   tapply(stay,{stateTime,treatment},@sem,{'continuous','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('state duration'); 
        ylabel('p[participation(t) = participation(t-1)]'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end

fig = figure; set(fig,'Name','participation21');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject  ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        % previousParticipation
        tmax=1;
        x5 = nan(numel(participation),tmax);
        for it=1:tmax
            x5(1+it:end,it) = participation(1:end-it);
            x5(nt<=it,it) = NaN;
        end
        stateTime = data.stateDuration(selection);
  
        nbin = 5;
        x = tapply(x5,{x5,treatment},@nanmean,{'discrete','discrete'},[nbin 2]);
        mu =  tapply(participation,{x5,treatment},@nanmean,{'discrete','discrete'},[nbin 2]);
        err =   tapply(participation,{x5,treatment},@sem,{'discrete','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('participation(t-1)'); 
        ylabel('participation(t)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
%%
    fig = figure; set(fig,'Name','participation22');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject  ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        nt = data.trialNumber(selection);
        % previousParticipation
        tmax=1;
        x5 = nan(numel(participation),tmax);
        for it=1:tmax
            x5(1+it:end,it) = participation(1:end-it);
            x5(nt<=it,it) = NaN;
        end
        stateTime = data.stateDuration(selection);
        stay = (participation==x5);
  
        nbin = 4;
        x = tapply(stateTime,{stateTime,x5,treatment},@nanmean,{'discrete','discrete','discrete'},[nbin 2 2]);
        mu =  tapply(participation,{stateTime,x5,treatment},@nanmean,{'discrete','discrete','discrete'},[nbin 2 2]);
        err =   tapply(participation,{stateTime,x5,treatment},@sem,{'discrete','discrete','discrete'},[nbin 2 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,1,i) ,mu(:,1,i), err(:,1,i),col{i}); h(i).LineStyle = '--';
           [~,~,h(i)] = errorscat(x(:,2,i) ,mu(:,2,i), err(:,2,i),col{i}); h(i).LineStyle = '-';

        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('state duration'); 
        ylabel('p[participation(t)]'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
    fig = figure; set(fig,'Name','participation23');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject  ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionChosen(selection);
        nt = data.trialNumber(selection);
        rt = data.responseTime(selection);

        % previousParticipation
        tmax=1;
        x5 = nan(numel(participation),tmax);
        for it=1:tmax
            x5(1+it:end,it) = participation(1:end-it);
            x5(nt<=it,it) = NaN;
        end
        stateTime = data.stateDuration(selection);
        stay = (participation==x5);
  
        nbin = 4;
        x = tapply(stateTime,{stateTime,x5,treatment},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
        mu =  tapply(rt,{stateTime,x5,treatment},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
        err =   tapply(rt,{stateTime,x5,treatment},@sem,{'continuous','discrete','discrete'},[nbin 2 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,1,i) ,mu(:,1,i), err(:,1,i),col{i}); %h(i).LineStyle = '--';
           [~,~,h(i)] = errorscat(x(:,2,i) ,mu(:,2,i), err(:,2,i),col{i}); h(i).LineStyle = '-';

        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('state duration'); 
        ylabel('response time (sec)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
%
%% 4. model selection
%
%
%