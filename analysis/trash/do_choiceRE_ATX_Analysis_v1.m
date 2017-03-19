%% do_choiceRE_ATX_Analysis_v1

%     clear all;
%     close all;
%     clear h ;
%     clear b;



%% local variables
    % data acquisition parameters
        forceSamplingFreq = 25; % Hz.
        rewardUnit2Volume = 0.9; % conversion rate between valve opening time (sec) & water volume (ml)
        drugName = 'atx';
        
        col = {[1 1 1]*0.5,[1 0.5 0]};
        
        subjectList = unique(data.subject);
        sessionList = unique(data.session);
        trtList = unique(data.treatment);
        trtSessions = {[1:5,15:18,24:29],[6:14,19:23]};
        sideList = {'left','right'};
        
        
%% renaming
    data.treatment = nominal(data.treatment);
    data.treatment(data.treatment=='free')='placebo';
    data.treatment(data.treatment=='atomoxetine')='atx';
    
    data.treatment(data.treatment=='placebo')='1_placebo';
    data.treatment(data.treatment=='atx')='2_atx';
    data.treatment = removecats(data.treatment);
    
    data.treatment(ismember(data.session,sessionList(trtSessions{1}))) = '1_placebo';
    data.treatment(ismember(data.session,sessionList(trtSessions{2}))) = '2_atx';

    ind = find(ismember(data.Properties.VariableNames,'gain'));
    data.Properties.VariableNames{ind} = 'ordinalRewardOutcome';

%% data selection
        selectionAll =  ones(height(data),1) ;
        selectionRepetition = (data.isRepeatedTrial==1) ;
        
        selectionMissed = ( data.errorType==1);
        selectionPremature = ( data.errorType==2);
        selectionSwitch = ( data.errorType==3);
        selectionIncorrect = ( data.errorType==4);

        selectionParticipate = ~selectionMissed ;
        selectionResponseTime =  ~selectionMissed & ~isnan(data.intertrialTime) ;
      
        selectionChosen = ( ~selectionMissed & ~selectionPremature & ~selectionSwitch  );
        selectionCorrect = ( selectionChosen  & ~selectionIncorrect );

%% data completion
        data.perf(isnan(data.perf)) = 0;

        data.volumeRewardLeft = data.cardinalRewardLeft.*data.rewardUnit.*rewardUnit2Volume;
        data.volumeRewardRight= data.cardinalRewardRight.*data.rewardUnit.*rewardUnit2Volume;

        
    
        data.chosenOrdinalReward = nan(numel(data.trialNumber),1);
        data.chosenOrdinalReward(selectionChosen & data.sideChoice==1) = data.ordinalRewardLeft(selectionChosen & data.sideChoice==1);
        data.chosenOrdinalReward(selectionChosen & data.sideChoice==2) = data.ordinalRewardRight(selectionChosen & data.sideChoice==2);
        data.chosenOrdinalEffort = nan(numel(data.trialNumber),1);
        data.chosenOrdinalEffort(selectionChosen & data.sideChoice==1) = data.ordinalEffortLeft(selectionChosen & data.sideChoice==1);
        data.chosenOrdinalEffort(selectionChosen & data.sideChoice==2) = data.ordinalEffortRight(selectionChosen & data.sideChoice==2);

        data.chosenCardinalReward = zeros(numel(data.trialNumber),1);
        data.chosenCardinalReward(selectionChosen & data.sideChoice==1) = data.cardinalRewardLeft(selectionChosen & data.sideChoice==1);
        data.chosenCardinalReward(selectionChosen & data.sideChoice==2) = data.cardinalRewardRight(selectionChosen & data.sideChoice==2);
        data.chosenCardinalEffort = zeros(numel(data.trialNumber),1);
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

        data.cumulativeReward = zeros(height(data),1);
        data.cumulativeEffort = zeros(height(data),1);
        data.cumulativeTime = zeros(height(data),1);
        data.normalizedTrialNumber = zeros(height(data),1);
        data.stateDuration = ones(height(data),1);


        for iSub = 1:numel(subjectList)
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            sessList = unique(data.session(selectSubject));
            for iSess = 1:numel(sessList)
                selectSession = (data.session==sessList(iSess)) ;
                select = selectSubject & selectSession;
                
                cumRew =   nancumsum(data.volumeRewardOutcome(select)); 
%                 cumRew(2:end) = cumRew(1:end-1);cumRew(1)= 0;
                data.cumulativeReward(select) = cumRew;
                
                cumPerf =  nancumsum(data.perf(select)); 
%                 cumPerf(2:end) = cumPerf(1:end-1);cumPerf(1)= 0;
                data.cumulativeEffort(select) = cumPerf;
                
                cumT1 =   nancumsum(data.intertrialTime(select)); 
                cumT2 =   nancumsum(data.offerTime(select)); 
                cumT3 =   nancumsum(data.responseTime(select)); 
                cumT4 =   nancumsum(data.volumeRewardOutcome(select)./rewardUnit2Volume); 

%                 cumT(2:end) = cumT(1:end-1);cumT(1)= 0;
                data.cumulativeTime(select) = cumT1 + cumT2 + cumT3 + cumT4 ;
                
                
                data.normalizedTrialNumber(select)  = data.trialNumber(select)./max(data.trialNumber(select)) ;
                
                stateDuration = ones(numel(data.trialNumber(select)),1);
                for it = 2:max(data.trialNumber(select))
                    if selectionChosen(select & data.trialNumber==it)==selectionChosen(select & data.trialNumber==it-1)
                        stateDuration(it) = stateDuration(it-1)+1;
                    else
                        stateDuration(it) = 1;
                    end
                end
                data.stateDuration(select) = stateDuration;

            end
        end
        
        data.rewardRate = data.cumulativeReward./data.cumulativeTime;
        data.effortRate = data.cumulativeEffort./data.cumulativeTime;

    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        selection = selectSubject & selectionChosen;
        calib = [unique(data.calibLeft(selection)),unique(data.calibRight(selection))];
        data.perf(selection) = data.force(selection)./(calib(data.sideChoice(selection))');
    end

%% data access
dataName = 'choiceRE_bob_aliosha_ATX_PLACEBO_19_05_2016';
% % load
% load(dataName,'data');
% % save
save(dataName,'data');
     

%% 1) Trial classification
    fig = figure; set(fig,'Name','n_trials');
    
    
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        errors = data.errorType(selectSubject & ~selectionRepetition);
        session = data.session(selectSubject & ~selectionRepetition);
        treatment = data.treatment(selectSubject & ~selectionRepetition);
        
        nErr = tapply(errors,{errors,treatment,session},@numel);
        nErr(nErr==0) = NaN;
        meanErr = nanmean(nErr,3);
        semErr = sem(nErr,3);
        
        for i = 1:2
          b(i) = bar( [1:4] - 0.25 + 0.5*(i-1) , [meanErr(1,i) meanErr(4,i) meanErr(2,i) meanErr(3,i)]) ;
          h = errbar([1:4] - 0.25 + 0.5*(i-1) ,...
              [meanErr(1,i) meanErr(4,i) meanErr(2,i) meanErr(3,i)],...
              [semErr(1,i) semErr(4,i) semErr(2,i) semErr(3,i)],...
                 'Color',[0.5 0.5 0.5],'LineWidth',2);
             b(i).FaceColor = col{i};
             b(i).EdgeColor = 'none';
             b(i).BarWidth = 0.3;
        end

        labels = {'chosen(correct)';'chosen(incorrect)';'unchosen(omitted)';'unchosen(premature)'};
        legend([b(1) b(2)],{'placebo','atx'});
        ax = gca; 
        ax.XTick = [1:4];
        ax.XTickLabels = labels;
        ax.XTickLabelRotation = 45;
        xlabel('trial classification'); 
        ylabel('number of trials'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
        
%% 2) Participation
     
fig = figure; set(fig,'Name','participation2');
    
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        errors = data.errorType(selectSubject & ~selectionRepetition);
        session = data.session(selectSubject & ~selectionRepetition);
        treatment = data.treatment(selectSubject & ~selectionRepetition);
        participation = selectionParticipate(selectSubject & ~selectionRepetition);
        
        x =  tapply(session,{session},@unique);
        mu =  tapply(participation,{session},@nanmean);
        err =   tapply(participation,{session},@sem);
        
        for i = 1:2
            ind = find(ismember(x,trtSessions{i}));
            h(i) = errorbar( ind ,mu(ind), err(ind) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = [1:numel(sessionList)];
        ax.XTickLabels = cellstr(sessionList');
        ax.XTickLabelRotation = 45;
        xlabel('session date');
        ylabel('participation rate'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    
    %%
    fig = figure; set(fig,'Name','participation_5');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        selection = selectSubject ;
        
        sessList = unique(data.session(selectSubject));
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        sumR = (data.cardinalRewardLeft(selection) + data.cardinalRewardRight(selection));
        sumR = round(sumR,3);
        nt = data.trialNumber(selection);
        
        nbin = 40;
        X = nan(numel(sessList),max(nt));
        Y = nan(numel(sessList),max(nt));
        S = nan(numel(sessList),max(nt));

        for is = 1:numel(sessList)
           selectSession = session==sessList(is);
           x = nt(selectSession);
%            x = x./max(x);
%            y = participation(session==sessList(is));
           % sliding average
%                % backward averaging
%                y = filter([repmat((1/nbin),1,nbin)],1,y);
%                for j = [1:(nbin-1)]
%                   y(j) = y(j)*(nbin)/j; 
%                end
               
               % backward-foreward averaging
               y = zeros(1,max(x));
               y0 = participation(selectSession);
               for j = [1:max(x)]
                     ind = [j-nbin/2-1:j+nbin/2]; ind = ind(ind>0); ind = ind(ind<=max(x));
                    y(j) = nanmean(y0(ind));
                    if isnan(y(j)); y(j)=0;end

               end
           
           c = col{ ismember(is,trtSessions{2}) + 1 };
%            scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%            plot(x,y,'LineWidth',2,'Color',c);
%            pause;
        % normalization of session length
%             ind = round(x*max(nt)/max(x));
%            X(is,ind) = x./max(x);
        % without normalization of session length
           ind = round(x);
           X(is,ind) = x;
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
           Y(is,isnan(Y(is,:))) = Y(is,find(isnan(Y(is,:)))-1 );
           end
        end
        
        for i =1:2            
            ind = find(ismember(sessList,sessionList(trtSessions{i})));
            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) ,...
                'alpha','transparency', 0.5 ); 
            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
%             h(i).LineStyle = '-'; 
            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.YLim = [0 1];
        xlabel('session progression (total trial)');
        ylabel('participation (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
    fig = figure; set(fig,'Name','rt_dynamic');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        selection = selectSubject & selectionChosen;
        
        sessList = unique(data.session(selectSubject));
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        sumR = (data.cardinalRewardLeft(selection) + data.cardinalRewardRight(selection));
        sumR = round(sumR,3);
        nt = data.trialNumber(selection);
        rt = data.responseTime(selection);
        
        nbin = 40;
        X = nan(numel(sessList),max(nt));
        Y = nan(numel(sessList),max(nt));
        S = nan(numel(sessList),max(nt));

        for is = 1:numel(sessList)
           selectSession = session==sessList(is);
           x = nt(selectSession);
%            x = x./max(x);
%            y = participation(session==sessList(is));
           % sliding average
%                % backward averaging
%                y = filter([repmat((1/nbin),1,nbin)],1,y);
%                for j = [1:(nbin-1)]
%                   y(j) = y(j)*(nbin)/j; 
%                end
               
               % backward-foreward averaging
               y = zeros(1,numel(x));
               y0 = rt(selectSession);
               for j = [1:numel(x)]
                     ind = [j-nbin/2-1:j+nbin/2]; ind = ind(ind>0); ind = ind(ind<=numel(x));
                    y(j) = nanmean(y0(ind));
                    if isnan(y(j)); y(j)=0;end

               end
           
           c = col{ ismember(is,trtSessions{2}) + 1 };
%            scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%            plot(x,y,'LineWidth',2,'Color',c);
%            pause;
        % normalization of session length
%             ind = round(x*max(nt)/max(x));
%            X(is,ind) = x./max(x);
        % without normalization of session length
           ind = round(x);
           X(is,ind) = x;
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
               missing = isnan(Y(is,:));
               neighbor = find(missing)-1;
               neighbor(neighbor==0) = 2;
               Y(is,missing) = Y(is,neighbor);
           end
        end
        
        for i =1:2            
            ind = find(ismember(sessList,sessionList(trtSessions{i})));
            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) ,...
                'alpha','transparency', 0.5 ); 
            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
%             h(i).LineStyle = '-'; 
            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('session progression (total trial)');
        ylabel('response time (s.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    
%%
    fig = figure; set(fig,'Name','participation_19');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
         selection = selectSubject ;
         
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        x1 = data.ordinalRewardLeft(selection);
        x2 = data.ordinalRewardRight(selection);
        x3 = data.ordinalEffortLeft(selection);
        x4 = data.ordinalEffortRight(selection);
        nt = data.trialNumber(selection);
        nnt = data.normalizedTrialNumber(selection);
        x5 = participation(1:end-1);  % previousParticipation
        x5 = [NaN ; x5];        x5(nt==1) = NaN;
        
        % fit
         predictor = table(rangescore(x1),rangescore(x2),rangescore(x3),rangescore(x4),rangescore(nnt),treatment ,(x5*2)-1, participation);
         predictor.Properties.VariableNames = {'R1';'R2';'E1';'E2';'ntrial';'treatment';'part_t_1';'participation'};
         formula = 'participation ~ 1 + part_t_1 + R1 + R2 + E1 + E2 + ntrial:(1 + part_t_1) + treatment:( 1 + part_t_1 + R1 + R2 + E1 + E2 + ntrial:(1 + part_t_1) )';
         predictor = predictor(selectSubject(selection),:);
         glm = fitglm( predictor, formula,'Distribution','binomial','link','logit','CategoricalVars',[6]);
         y2  = glm.Fitted.Response;


        nbin = 10;
        x = tapply(nnt,{nnt,x5,treatment},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
        mu = tapply(participation,{nnt,x5,treatment},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
        err = tapply(participation,{nnt,x5,treatment},@sem,{'continuous','discrete','discrete'},[nbin 2 2]);
        mu2 = tapply(y2,{nnt,x5,treatment},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);

        subplot(1,numel(subjectList),+iSub);hold on
        for i =1:2
            
%             subplot(2,numel(subjectList),(i-1)*2+iSub);hold on

           [~,~,h(1)] = errorscat(x(:,1,i) ,mu(:,1,i), err(:,1,i),col{i});
           h(1).YData = mu2(:,1,i);
           [~,~,h(2)] = errorscat(x(:,2,i) ,mu(:,2,i), err(:,2,i),col{i});
           h(2).YData = mu2(:,2,i);

            h(1).LineStyle='--';
            
            
            legend([h(1) h(2)],{'part(t-1)=0 ','part(t-1)=1'});
            
            ax = gca;
            ax.YLim = [0 1];
            xlabel('session progression (total trial)');
            ylabel('participation (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
           
        end
    end
    
      %%
    fig = figure; set(fig,'Name','participation_15');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
%          selection = selectSubject & ~selectionRepetition;
        selection = selectSubject ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        sumR = (data.cardinalRewardLeft(selection) + data.cardinalRewardRight(selection));
        sumR = round(sumR,3);
        nt = data.trialNumber(selection);
        
        nbin = 40;
        X = nan(numel(sessionList),max(nt));
        Y = nan(numel(sessionList),max(nt));
        S = nan(numel(sessionList),max(nt));

        for is = 1:numel(sessionList)
           x = nt(session==sessionList(is));
           x2 = sumR(session==sessionList(is));
           y2 = participation(session==sessionList(is));
           
%            x = x./max(x);
%            y = participation(session==sessionList(is));
%            % sliding average
%            y = filter([repmat((1/nbin),1,nbin)],1,y);
%            for j = [1:(nbin-1)]
%               y(j) = y(j)*(nbin)/j; 
%            end
           % sliding fit
           y = zeros(1,max(x));
           for j = [1:max(x)]
%                  ind = [j-nbin+1:j]; ind = ind(ind>0);
%                  ind = [j:j+nbin-1]; ind = ind(ind<=max(x));
                 ind = [j-nbin/2-1:j+nbin/2]; ind = ind(ind>0); ind = ind(ind<=max(x));

                [beta,~,stat] = glmfit(rangescore(x2(ind)),y2(ind),'binomial','link','logit');
%                 y(j) = beta(2)*(numel(ind)/nbin);
                y(j) = stat.t(2);
                if isnan(stat.t(2)); y(j)=0;end
%                 y(j) = beta(2);

           end
           
           c = col{ ismember(is,trtSessions{2}) + 1 };
%            scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%            plot(x,y,'LineWidth',2,'Color',c);
%            pause;
        % normalization of session length
%             ind = round(x*max(nt)/max(x));
%            X(is,ind) = x./max(x);
        % without normalization of session length
           ind = round(x);
           X(is,ind) = x;
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
           Y(is,isnan(Y(is,:))) = Y(is,find(isnan(Y(is,:)))-1 );
           end
        end
        
        for i =1:2
            ind = trtSessions{i};

            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
%             h(i) = plot(x,mu,'LineWidth',2,'Color',col{i});
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 
%             [h(i),hp] = boundedline( [1:sum(~isnan(x))]./sum(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 );   

            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
            h(i).LineStyle = '-'; 
            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.YLim = [0 1];
        xlabel('session progression (% total trial)');
        ylabel('sensitivity of participation to reward (t-stat)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
    fig = figure; set(fig,'Name','participation_16');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
%          selection = selectSubject & ~selectionRepetition;
        selection = selectSubject ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        sumE = (data.cardinalEffortLeft(selection) + data.cardinalEffortRight(selection));
        sumE = round(sumE,3);
        nt = data.trialNumber(selection);
        
        nbin = 50;
        X = nan(numel(sessionList),max(nt));
        Y = nan(numel(sessionList),max(nt));
        S = nan(numel(sessionList),max(nt));

        for is = 1:numel(sessionList)
           x = nt(session==sessionList(is));
           x2 = sumE(session==sessionList(is));
           y2 = participation(session==sessionList(is));
           
%            x = x./max(x);
%            y = participation(session==sessionList(is));
%            % sliding average
%            y = filter([repmat((1/nbin),1,nbin)],1,y);
%            for j = [1:(nbin-1)]
%               y(j) = y(j)*(nbin)/j; 
%            end
           % sliding fit
           y = zeros(1,max(x));
           for j = [1:max(x)]
%                  ind = [j-nbin+1:j]; ind = ind(ind>0);
%                  ind = [j:j+nbin-1]; ind = ind(ind<=max(x));
                 ind = [j-nbin/2-1:j+nbin/2]; ind = ind(ind>0); ind = ind(ind<=max(x));

                [beta,~,stat] = glmfit(zscore(x2(ind)),y2(ind),'binomial','link','logit');
%                 y(j) = beta(2)*(numel(ind)/nbin);
                y(j) = stat.t(2);
                if isnan(stat.t(2)); y(j)=0;end
%                 y(j) = beta(2);

           end
           
           c = col{ ismember(is,trtSessions{2}) + 1 };
%            scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%            plot(x,y,'LineWidth',2,'Color',c);
%            pause;
            ind = round(x*max(nt)/max(x));
           X(is,ind) = x./max(x);
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
           Y(is,isnan(Y(is,:))) = Y(is,find(isnan(Y(is,:)))+1 );
           end
        end
        
        for i =1:2
            ind = trtSessions{i};

            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
%             h(i) = plot(x,mu,'LineWidth',2,'Color',col{i});
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 
%             [h(i),hp] = boundedline( [1:sum(~isnan(x))]./sum(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 );   

            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
            h(i).LineStyle = '-'; 
            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.YLim = [0 1];
        xlabel('session progression (% total trial)');
        ylabel('sensitivity of participation to effort (t-stat)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    %%
    fig = figure; set(fig,'Name','participation_17');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
%          selection = selectSubject & ~selectionRepetition;
%         selection = selectSubject ;
        selection = selectSubject & selectionChosen ;

        session = data.session(selection);
        treatment = data.treatment(selection);
        dE = data.deltaCardinalEffort(selection);
        choice = data.sideChoice(selection)-1;

        nt = data.trialNumber(selection);
        
        nbin = 50;
        X = nan(numel(sessionList),max(nt));
        Y = nan(numel(sessionList),max(nt));
        S = nan(numel(sessionList),max(nt));

        for is = 1:numel(sessionList)
           x = nt(session==sessionList(is));
           x2 = dE(session==sessionList(is));
           y2 = choice(session==sessionList(is));
           
%            x = x./max(x);
%            y = participation(session==sessionList(is));
%            % sliding average
%            y = filter([repmat((1/nbin),1,nbin)],1,y);
%            for j = [1:(nbin-1)]
%               y(j) = y(j)*(nbin)/j; 
%            end
           % sliding fit
           y = zeros(1,numel(x));
           for j = [1:numel(x)]
%                  ind = [j-nbin+1:j]; ind = ind(ind>0);
%                  ind = [j:j+nbin-1]; ind = ind(ind<=max(x));
                 ind = [j-nbin/2-1:j+nbin/2]; ind = ind(ind>0); ind = ind(ind<=numel(x));

                [beta,~,stat] = glmfit(zscore(x2(ind)),y2(ind),'binomial','link','logit');
%                 y(j) = beta(2)*(numel(ind)/nbin);
                y(j) = stat.t(2);
                if isnan(stat.t(2)); y(j)=0;end
%                 y(j) = beta(2);

           end
           
           c = col{ ismember(is,trtSessions{2}) + 1 };
%            % TO COMMENT
%                scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%                plot(x,y,'LineWidth',2,'Color',c);
%                pause;
            ind = round(x*max(nt)/max(x));
           X(is,ind) = x./max(x);
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
           Y(is,isnan(Y(is,:))) = Y(is,find(isnan(Y(is,:)))+1 );
           end
        end
        
        for i =1:2
            ind = trtSessions{i};

            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
%             h(i) = plot(x,mu,'LineWidth',2,'Color',col{i});
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 
%             [h(i),hp] = boundedline( [1:sum(~isnan(x))]./sum(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 );   

            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
            h(i).LineStyle = '-'; 
            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
%         ax.YLim = [0 1];
        xlabel('session progression (% total trial)');
        ylabel('sensitivity of choice to effort (t-stat)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
      %%
    fig = figure; set(fig,'Name','participation_18');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
%          selection = selectSubject & ~selectionRepetition;
%         selection = selectSubject ;
        selection = selectSubject & selectionChosen ;

        session = data.session(selection);
        treatment = data.treatment(selection);
        dR = data.deltaCardinalReward(selection);
        choice = data.sideChoice(selection)-1;

        nt = data.trialNumber(selection);
        
        nbin = 50;
        X = nan(numel(sessionList),max(nt));
        Y = nan(numel(sessionList),max(nt));
        S = nan(numel(sessionList),max(nt));

        for is = 1:numel(sessionList)
           x = nt(session==sessionList(is));
           x2 = dR(session==sessionList(is));
           y2 = choice(session==sessionList(is));
           
%            x = x./max(x);
%            y = participation(session==sessionList(is));
%            % sliding average
%            y = filter([repmat((1/nbin),1,nbin)],1,y);
%            for j = [1:(nbin-1)]
%               y(j) = y(j)*(nbin)/j; 
%            end
           % sliding fit
           y = zeros(1,numel(x));
           for j = [1:numel(x)]
%                  ind = [j-nbin+1:j]; ind = ind(ind>0);
%                  ind = [j:j+nbin-1]; ind = ind(ind<=max(x));
                 ind = [j-nbin/2-1:j+nbin/2]; ind = ind(ind>0); ind = ind(ind<=numel(x));

                [beta,~,stat] = glmfit(zscore(x2(ind)),y2(ind),'binomial','link','logit');
%                 y(j) = beta(2)*(numel(ind)/nbin);
                y(j) = stat.t(2);
                if isnan(stat.t(2)); y(j)=0;end
%                 y(j) = beta(2);

           end
           
           c = col{ ismember(is,trtSessions{2}) + 1 };
%            % TO COMMENT
%                scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%                plot(x,y,'LineWidth',2,'Color',c);
%                pause;
            ind = round(x*max(nt)/max(x));
           X(is,ind) = x./max(x);
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
           Y(is,isnan(Y(is,:))) = Y(is,find(isnan(Y(is,:)))+1 );
           end
        end
        
        for i =1:2
            ind = trtSessions{i};

            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
%             h(i) = plot(x,mu,'LineWidth',2,'Color',col{i});
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 
%             [h(i),hp] = boundedline( [1:sum(~isnan(x))]./sum(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 );   

            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
            h(i).LineStyle = '-'; 
            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
%         ax.YLim = [0 1];
        xlabel('session progression (% total trial)');
        ylabel('sensitivity of choice to reward (t-stat)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    %%
    fig = figure; set(fig,'Name','participation_11');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
%          selection = selectSubject & ~selectionRepetition;
        selection = selectSubject ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        reward = data.volumeRewardOutcome(selection);
        effort = data.perf(selection);
        nt = data.trialNumber(selection);
        
        nbin = 40;
        X = nan(numel(sessionList),max(nt));
        Y = nan(numel(sessionList),max(nt));
        S = nan(numel(sessionList),max(nt));

        for is = 1:numel(sessionList)
           x = nt(session==sessionList(is));
%            x = x./max(x);
           y = reward(session==sessionList(is));
           % sliding average
           y = filter([repmat((1/nbin),1,nbin)],1,y);
           for j = [1:(nbin-1)]
              y(j) = y(j)*(nbin)/j; 
           end
%            c = col{ ismember(is,trtSessions{2}) + 1 };
%            scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%            plot(x,y,'LineWidth',2,'Color',c);
%            pause;
            ind = round(x*max(nt)/max(x));
           X(is,ind) = x./max(x);
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
           Y(is,isnan(Y(is,:))) = Y(is,find(isnan(Y(is,:)))+1 );
           end
        end
        
        for i =1:2
            ind = trtSessions{i};

            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
%             h(i) = plot(x,mu,'LineWidth',2,'Color',col{i});
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 
%             [h(i),hp] = boundedline( [1:sum(~isnan(x))]./sum(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 );   

            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
            h(i).LineStyle = '-'; 
            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
%         ax.YLim = [0 1];
        xlabel('session progression (% total trial)');
        ylabel('reward (ml)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    %%
        fig = figure; set(fig,'Name','participation_12');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        selection = selectSubject & selectionChosen;

        sessList = unique(data.session(selection));
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        effort = data.perf(selection);
        nt = data.trialNumber(selection);
        
        nbin = 40;
        X = nan(numel(sessList),max(nt));
        Y = nan(numel(sessList),max(nt));
        S = nan(numel(sessList),max(nt));

        for is = 1:numel(sessList)
           selectSession = session==sessList(is);
           x = nt(selectSession);
%            x = x./max(x);% backward-foreward averaging
           y = zeros(1,numel(x));
           y0 = effort(selectSession);
           for j = [1:numel(x)]
                 ind = [j-nbin/2-1:j+nbin/2]; ind = ind(ind>0); ind = ind(ind<=numel(x));
                y(j) = nanmean(y0(ind));
                if isnan(y(j)); y(j)=0;end

           end
%            c = col{ ismember(is,trtSessions{2}) + 1 };
%            scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%            plot(x,y,'LineWidth',2,'Color',c);
%            pause;
%             ind = round(x*max(nt)/max(x));
           ind = round(x);clc
           X(is,ind) = x;
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
           Y(is,isnan(Y(is,:))) = Y(is,find(isnan(Y(is,:)))-1 );
           end
        end
        
        subplot(1,numel(subjectList),iSub);hold on
        for i =1:2
            ind = find(ismember(sessList,sessionList(trtSessions{i})));

            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 
            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
%         ax.YLim = [0 1];
        xlabel('session progression (% total trial)');
        ylabel('force (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
      %%
        fig = figure; set(fig,'Name','participation_14');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
%          selection = selectSubject & selectionChosen;
        selection = selectSubject ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        reward = data.volumeRewardOutcome(selection);
        effort = data.cumulativeEffort(selection);
        nt = data.trialNumber(selection);
        
        nbin = 1;
        X = nan(numel(sessionList),max(nt));
        Y = nan(numel(sessionList),max(nt));
        S = nan(numel(sessionList),max(nt));

        for is = 1:numel(sessionList)
           x = nt(session==sessionList(is));
%            x = x./max(x);
           y = effort(session==sessionList(is));
           % sliding average
           y = filter([repmat((1/nbin),1,nbin)],1,y);
           for j = [1:(nbin-1)]
              y(j) = y(j)*(nbin)/j; 
           end
%            c = col{ ismember(is,trtSessions{2}) + 1 };
%            scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%            plot(x,y,'LineWidth',2,'Color',c);
%            pause;
            ind = round(x*max(nt)/max(x));
           X(is,ind) = x./max(x);
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
           Y(is,isnan(Y(is,:))) = Y(is,find(isnan(Y(is,:)))+1 );
           end
        end
        
        for i =1:2
            ind = trtSessions{i};

            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
%             h(i) = plot(x,mu,'LineWidth',2,'Color',col{i});
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 
%             [h(i),hp] = boundedline( [1:sum(~isnan(x))]./sum(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 );   

            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
            h(i).LineStyle = '-'; 
            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
%         ax.YLim = [0 1];
        xlabel('session progression (% total trial)');
        ylabel('cumulative force (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
 %%
fig = figure; set(fig,'Name','participation_13');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
         selection = selectSubject & selectionChosen;
%         selection = selectSubject ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        correct = selectionCorrect(selection);
        reward = data.volumeRewardOutcome(selection);
        effort = data.perf(selection);
        nt = data.trialNumber(selection);
        
        nbin = 40;
        X = nan(numel(sessionList),max(nt));
        Y = nan(numel(sessionList),max(nt));
        S = nan(numel(sessionList),max(nt));

        for is = 1:numel(sessionList)
           x = nt(session==sessionList(is));
%            x = x./max(x);
           y = correct(session==sessionList(is));
           % sliding average
           y = filter([repmat((1/nbin),1,nbin)],1,y);
           for j = [1:(nbin-1)]
              y(j) = y(j)*(nbin)/j; 
           end
%            c = col{ ismember(is,trtSessions{2}) + 1 };
%            scatter(x,y,'LineWidth',2,'MarkerFaceColor',c,'MarkerEdgeColor',c);
%            plot(x,y,'LineWidth',2,'Color',c);
%            pause;
            ind = round(x*max(nt)/max(x));
           X(is,ind) = x./max(x);
           Y(is,ind) = y;
           while sum(isnan(Y(is,:)))>0
           Y(is,isnan(Y(is,:))) = Y(is,find(isnan(Y(is,:)))+1 );
           end
        end
        
        for i =1:2
            ind = trtSessions{i};

            x = nanmean(X(ind,:),1);
            err = sem(Y(ind,:),1);
            mu = nanmean(Y(ind,:),1);
            
%             h(i) = plot(x,mu,'LineWidth',2,'Color',col{i});
            [h(i),hp] = boundedline( x(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 ); 
%             [h(i),hp] = boundedline( [1:sum(~isnan(x))]./sum(~isnan(x)) , mu(~isnan(x)) , err(~isnan(x)) , 'alpha','transparency', 0.5 );   

            cl = col{i};
            set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
            h(i).LineStyle = '-'; 
            
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
%         ax.YLim = [0 1];
        xlabel('session progression (% total trial)');
        ylabel('correct force (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    %%
fig = figure; set(fig,'Name','participation_19');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
%          selection = selectSubject & selectionChosen;
        selection = selectSubject ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionChosen(selection);
        nt = data.trialNumber(selection);
        
        nbin = 40;
        X1 = nan(numel(sessionList),1);
        X2 = nan(numel(sessionList),1);

        for is = 1:numel(sessionList)
           x = nt(session==sessionList(is));
           y = participation(session==sessionList(is));
           % sliding classification
           work = 0; iw=1;
           leisure = 0; il=1;
           for j = [1:numel(x)]
               if j>1 
                   if y(j) > y(j-1)
                       iw=iw+1;
                   elseif y(j) < y(j-1)
                       il=il+1;
                   end
               end
               work = work + y(j);
               leisure = leisure + (1-y(j));
               X1(is,iw) = work;
               X2(is,il) = leisure;
           end
        end
        X1(X1==0)=NaN;
        X2(X2==0)=NaN;
        
        for i =1:2
            ind = trtSessions{i};

            x = [1- 0.25 + 0.5*(i-1)];
            mu = nanmean(nanmean(X1(ind,:),2),1);
            err = sem(nanmean(X1(ind,:),2),1);
            b(i) = bar( x , mu) ;
            h = errbar( x ,mu, err ,...
                 'Color',[1 1 1].*0,'LineWidth',2);
             b(i).FaceColor = col{i};
             b(i).EdgeColor = 'none';
             b(i).BarWidth = 0.3;
             
            x = [3- 0.25 + 0.5*(i-1)];
            mu = nanmean(nanmean(X2(ind,:),2),1);
            err = sem(nanmean(X2(ind,:),2),1);
            b(i) = bar( x , mu) ;
            h = errbar( x ,mu, err ,...
                 'Color',[1 1 1].*0,'LineWidth',2);
             b(i).FaceColor = col{i};
             b(i).EdgeColor = 'none';
             b(i).BarWidth = 0.3;
        end

        ax = gca;
        ax.XTick = [ 1 3 ];
        ax.XTickLabel = {'participation','pause'};

        ylabel(' period length (trial number)'); 
        yy =ylim;
        legend([b(1) b(2)],{'placebo','atx'});
        t = title(['monkey ' sub(1) ]); 
    end

    
     %%
    fig = figure; set(fig,'Name','participation_6');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        selection = selectSubject ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        nt = data.trialNumber(selection);
        rr = data.cumulativeReward(selection);
        ee = data.cumulativeEffort(selection);
        
        nbin = 10;
        x = tapply(nt,{nt,session},@nanmean,{'continuous','discrete'},[nbin nbin]);
        y =  tapply(rr,{nt,treatment,session},@nanmean,{'continuous','discrete','discrete'},[nbin 2 nbin]);
        err =   tapply(rr,{nt,treatment,session},@sem,{'continuous','discrete','discrete'},[nbin 2 nbin]);
        x = nanmean(x,2);
        mu = nanmean(y,3);
        err = sem(y,3);
        
        for i = 1:2
            ind = [1:numel(x)];
            h(i) = errorbar( x ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
%             for j = 1:size(y,3)
%                 h(i) = plot( x , y(ind,i,j),'Color',col{i},'LineWidth',2 );
%             end
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('trial number');
        ylabel('cumulative reward'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
%%
 fig = figure; set(fig,'Name','cumReward');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
            nt = data.trialNumber(selection);
            rr = data.cumulativeReward(selection);
            ee = data.cumulativeEffort(selection);

        
        % compute stats
            nbin = 10;
            xx = tapply(nt,{nt,treatment,session},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
            yy =  tapply(rr,{nt,treatment,session},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
            x = nanmean(xx,3);
            y = nanmean(yy,3);
            z = sem(yy,3);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
               cl = col{i};
               [h(i),hp] = boundedline( x(:,i) , y(:,i) , z(:,i) ,...
                                    'alpha','transparency', 0.5 ); 
                set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
              
            end
            
            
            
        % legending
            legend([h(1) h(2)],{'placebo','atx'});
            ax = gca; 
%             ax.TickLength = [0 0];
%             ax.XTick = [];
            xlabel('session progression (total trial)');
            ylabel('cumulated reward (ml)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
    
%%
     fig = figure; set(fig,'Name','cumEffort');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
            nt = data.trialNumber(selection);
            rr = data.cumulativeReward(selection);
            ee = data.cumulativeEffort(selection);

        
        % compute stats
            nbin = 10;
            xx = tapply(nt,{nt,treatment,session},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
            yy =  tapply(ee,{nt,treatment,session},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
            x = nanmean(xx,3);
            y = nanmean(yy,3);
            z = sem(yy,3);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
               cl = col{i};
               [h(i),hp] = boundedline( x(:,i) , y(:,i) , z(:,i) ,...
                                    'alpha','transparency', 0.5 ); 
                set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
              
            end
            
            
            
        % legending
            legend([h(1) h(2)],{'placebo','atx'});
            ax = gca; 
%             ax.TickLength = [0 0];
%             ax.XTick = [];
            xlabel('session progression (total trial)');
            ylabel('cumulated force (%fmax)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end

    
%%
fig = figure; set(fig,'Name','participation10');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
            nt = data.trialNumber(selection);
            rr = data.cumulativeReward(selection);
            ee = data.cumulativeEffort(selection);
            tt = data.cumulativeTime(selection)./60;
            fmax = data.perf(selection);
            
        % compute stats
        varMat = [ rr , ee , fmax , nt, tt] ;
        yList = {'total reward (ml)','total force (%fmax)','max force (%fmax)',...
            'total number of trials','total time (min)'};
            for ivar=1:5
                var = varMat(:,ivar);
                yy =  tapply(var,{session,treatment},@max);
                z = sem(yy,1);
                y = nanmean(yy,1);
                [h,p] = ttest2(exnan(yy(:,1)),exnan(yy(:,2)));

            % display
                subplot(5,numel(subjectList),(ivar-1)*2+iSub);hold on
                for i = 1:2
                   x = [1- 0.25 + 0.5*(i-1)];
                  [ b(i),~ ] = barplot( x ,y(i), z(i) , col{i} );
                  b(i).BarWidth = 0.3;
                end
                if p<=0.05; s = sigstar({[0.75,1.25]},p); end

            % legending
                if ivar==1;
                    legend([b(1) b(2)],{'placebo','atx'});
                     t = title(['monkey ' sub(1) ]); 
                end
                ylabel(yList{ivar});
                ax = gca; 
                ax.TickLength = [0 0];
                ax.XTick = [];
             end
            
    end
    
    %%
      fig = figure; set(fig,'Name','participation_7');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
%          selection = selectSubject & ~selectionRepetition;
        selection = selectSubject ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        nt = data.trialNumber(selection);
        rr = data.cumulativeReward(selection);
        ee = data.cumulativeEffort(selection);
        
        nbin = 10;
        x = tapply(nt,{nt,session},@nanmean,{'continuous','discrete'},[nbin nbin]);
        y =  tapply(ee,{nt,treatment,session},@nanmean,{'continuous','discrete','discrete'},[nbin 2 nbin]);
        err =   tapply(ee,{nt,treatment,session},@sem,{'continuous','discrete','discrete'},[nbin 2 nbin]);
        x = nanmean(x,2);
        mu = nanmean(y,3);
        err = sem(y,3);
        
        for i = 1:2
            ind = [1:numel(x)];
            h(i) = errorbar( x ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
%             for j = 1:size(y,3)
%                 h(i) = plot( x , y(ind,i,j),'Color',col{i},'LineWidth',2 );
%             end
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('trial number');
        ylabel('cumulative effort'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
fig = figure; set(fig,'Name','participation_8');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject  ;
        session = data.session(selection);
        novelty = (~selectionRepetition(selection));

        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
  
        mu =  tapply(participation,{novelty,treatment},@nanmean);
        err =   tapply(participation,{novelty,treatment},@sem);
        
       for i = 1:2
            ind = [1:2];
%             h(i) = scatter( x(ind,i) ,mu(ind,i),10,...
%                      'MarkerFaceColor',col{i});
             h(i) = plot( ind ,mu(ind,i),...
                     'Color',col{i},'LineWidth',2);
             errbar( ind ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
       end
%         plot(x(ind,i),x(ind,i),'k--');
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = [1 2];
        ax.XLim = [0 3];
        ax.XTickLabel = {'repetition','novel'};
        ylabel('participation (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end

    
%% compare predictors of participation    
fig = figure; set(fig,'Name','participation_predictors1');
    
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        selection = selectSubject & ~selectionRepetition;
        selection = selectSubject;

        errors = data.errorType(selection);
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        
        
        sumR = (data.cardinalRewardLeft(selection) + data.cardinalRewardRight(selection));
        sumR = round(sumR,3);
        maxR = max([data.cardinalRewardLeft(selection) , data.cardinalRewardRight(selection)],[],2);
        minR = min([data.cardinalRewardLeft(selection) , data.cardinalRewardRight(selection)],[],2);
        delatR = abs(data.cardinalRewardRight(selection) - data.cardinalRewardLeft(selection));
        delatR = round(delatR,3);
        var = {'sumR','maxR','minR','delatR'};
       
        
        for ip = 1:4
            subplot(numel(subjectList),4,ip+4*(iSub-1));hold on

            
           eval([ 'predictor =  ' var{ip} ';' ]);
           x =  tapply(predictor,{predictor,treatment},@nanmean);
           mu =  tapply(participation,{predictor,treatment},@nanmean);
           err =   tapply(participation,{predictor,treatment},@sem);
            
            for i = 1:2
                ind = [1:numel(unique(predictor))];
                h(i) = errorbar( x(ind,i) ,mu(ind,i), err(ind,i) ,...
                         'Color',col{i},'LineWidth',2);
            end
%             legend([h(1) h(2)],{'placebo','atx'});
            ax = gca; 
            ax.YLim = [0.5 1];
            ylabel('participation rate'); 
            xlabel(var{ip}); 
            yy =ylim;
            if ip ==1; t = title(['monkey ' sub(1) ]); end
        end
        
    end
    
fig = figure; set(fig,'Name','participation_predictors2');
    
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        selection = selectSubject & ~selectionRepetition;
        selection = selectSubject;

        errors = data.errorType(selection);
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        
        
        sumE = (data.cardinalEffortLeft(selection) + data.cardinalEffortRight(selection));
        sumE = round(sumE,3);
        maxE = max([data.cardinalEffortLeft(selection) , data.cardinalEffortRight(selection)],[],2);
        minE = min([data.cardinalEffortLeft(selection) , data.cardinalEffortRight(selection)],[],2);
        delatE = abs(data.cardinalEffortRight(selection) - data.cardinalEffortLeft(selection));
        delatE = round(delatE,3);
        var = {'sumE','maxE','minE','delatE'};
       
        
        for ip = 1:4
            subplot(numel(subjectList),4,ip+4*(iSub-1));hold on

            
           eval([ 'predictor =  ' var{ip} ';' ]);
           x =  tapply(predictor,{predictor,treatment},@nanmean);
           mu =  tapply(participation,{predictor,treatment},@nanmean);
           err =   tapply(participation,{predictor,treatment},@sem);
            
            for i = 1:2
                ind = [1:numel(unique(predictor))];
                h(i) = errorbar( x(ind,i) ,mu(ind,i), err(ind,i) ,...
                         'Color',col{i},'LineWidth',2);
            end
%             legend([h(1) h(2)],{'placebo','atx'});
            ax = gca; 
            ax.YLim = [0.5 1];
            ylabel('participation rate'); 
            xlabel(var{ip}); 
            yy =ylim;
            if ip ==1; t = title(['monkey ' sub(1) ]); end
        end
        
    end
 
%%
fig = figure; set(fig,'Name','participation_predictors3');
    
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        selection = selectSubject & ~selectionRepetition;

        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        errors = data.errorType(selection);
        session = data.session(selection);
        treatment = data.treatment(selection);
        participation = selectionParticipate(selection);
        
        sumR = (data.cardinalRewardLeft(selection) + data.cardinalRewardRight(selection));
        sumR = round(sumR,3);
        sumE = (data.cardinalEffortLeft(selection) + data.cardinalEffortRight(selection));
        sumE = round(sumE,3);
        delatR = (data.cardinalRewardRight(selection) - data.cardinalRewardLeft(selection));
        delatE = (data.cardinalEffortRight(selection) - data.cardinalEffortLeft(selection));
        delatR = round(delatR,3);
        delatE = round(delatE,3);
        congruenceRE = delatR.*(-delatE);
        congruenceRE = sumR.*(sumE);

        
        x =  tapply(congruenceRE,{congruenceRE,treatment},@nanmean,{'continuous','discrete'},[6 2]);
        mu =  tapply(participation,{congruenceRE,treatment},@nanmean,{'continuous','discrete'},[6 2]);
        err =   tapply(participation,{congruenceRE,treatment},@sem,{'continuous','discrete'},[6 2]);
        
        
        for i = 1:2
          b(i) = bar( x(:,i) + 0.25*(i-1) , mu(:,i)) ;
          h = errbar( x(:,i) + 0.25*(i-1) ,mu(:,i), err(:,i) ,...
                 'Color',[0.5 0.5 0.5],'LineWidth',2);
             b(i).FaceColor = col{i};
             b(i).EdgeColor = 'none';
             b(i).BarWidth = 0.5;
        end

        legend([b(1) b(2)],{'placebo','atx'});
        ax = gca; 
        ax.YLim = [0.5 1];
        ylabel('participation rate'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    
%% 3) Correct execution
    
fig = figure; set(fig,'Name','correct_force');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen ;
       
        % variables
            errors = data.errorType(selection);
            session = data.session(selection);
            treatment = data.treatment(selection);
            correct = selectionCorrect(selection);
        
        % compute stats
            mu =  tapply(correct,{session,treatment},@nanmean);
            [h,p] = ttest2(exnan(mu(:,1)),exnan(mu(:,2)));
            err =   sem(mu,1);
            mu  =   nanmean(mu,1);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [ b(i),~ ] = barplot( 1- 0.25 + 0.5*(i-1) ,mu(i), err(i) , col{i} );
              b(i).BarWidth = 0.3;
            end
            if p<=0.05 ; s = sigstar({[0.75,1.25]},p); end
            
            
            
        % legending
            legend([b(1) b(2)],{'placebo','atx'});
            ax = gca; 
            ax.YLim = [0.5 1];
            ax.TickLength = [0 0];
            ax.XTick = [];

            ylabel('correct force (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
    

%%
fig = figure; set(fig,'Name','correct execution2');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        correct = selectionCorrect(selectSubject & selectionChosen &  ~selectionRepetition);
  
        mu =  tapply(correct,{session},@nanmean);
        err =   tapply(correct,{session},@sem);
        
        for i = 1:2
            ind = trtSessions{i};
            h(i) = errorbar( ind ,mu(ind), err(ind) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = [1:numel(sessionList)];
        ax.XTickLabels = cellstr(sessionList');
        ax.XTickLabelRotation = 45;
        xlabel('session date');
        ylabel('correct force execution (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%% 4) Force execution
fig = figure; set(fig,'Name','force execution');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;

        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        
        mu =  tapply(force,{treatment},@nanmean);
        err =   tapply(force,{treatment},@sem);
        
        
        for i = 1:2
          b(i) = bar( 1- 0.25 + 0.5*(i-1) , mu(i)) ;
          h = errbar( 1- 0.25 + 0.5*(i-1) ,mu(i), err(i) ,...
                 'Color',[0.5 0.5 0.5],'LineWidth',2);
             b(i).FaceColor = col{i};
             b(i).EdgeColor = 'none';
             b(i).BarWidth = 0.3;
             
        end

        legend([b(1) b(2)],{'placebo','atx'});
        ax = gca; 
        ylabel('force peak (% fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
    fig = figure; set(fig,'Name','force');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen ;
       
        % variables
            errors = data.errorType(selection);
            session = data.session(selection);
            treatment = data.treatment(selection);
            force = data.perf(selection);
            effort = data.chosenCardinalEffort(selection);

        % compute stats
            mu =  tapply(force,{session,treatment,effort},@nanmean);
            mu = nanmean(mu,3);
            mu =  tapply(force,{session,treatment},@nanmean);
            [h,p] = ttest2(exnan(mu(:,1)),exnan(mu(:,2)));
            err =   sem(mu,1);
            mu  =   nanmean(mu,1);
            
        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [ b(i),~ ] = barplot( 1- 0.25 + 0.5*(i-1) ,mu(i), err(i) , col{i} );
              b(i).BarWidth = 0.3;
            end
            if p<=0.05 ; s = sigstar({[0.75,1.25]},p); end
            
            
            
        % legending
            legend([b(1) b(2)],{'placebo','atx'});
            ax = gca; 
%             ax.YLim = [0.5 1];
            ax.TickLength = [0 0];
            ax.XTick = [];

            ylabel('force (%fmax)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
    
    
%%
fig = figure; set(fig,'Name','force execution2');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        force = data.perf(selectSubject & selectionChosen &  ~selectionRepetition);
  
        mu =  tapply(force,{session},@nanmean);
        err =   tapply(force,{session},@sem);
        
        for i = 1:2
            ind = trtSessions{i};
            h(i) = errorbar( ind ,mu(ind), err(ind) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = [1:numel(sessionList)];
        ax.XTickLabels = cellstr(sessionList');
        ax.XTickLabelRotation = 45;
        xlabel('session date');
        ylabel('force peak (% fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
fig = figure; set(fig,'Name','force');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        force = data.force(selectSubject & selectionChosen &  ~selectionRepetition);
        choice = data.sideChoice(selectSubject & selectionChosen & ~selectionRepetition)-1;

        
        x = [1,2];
        mu =  tapply(force,{choice,treatment},@nanmean);
        err =   tapply(force,{choice,treatment},@sem);
        
        for i = 1:2
            ind = [1:2];
            h(i) = errorbar( x ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = x;
        ax.XTickLabels = {'left','right'};
        xlabel(texlabel('side'));
        ylabel('force peak (au)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end    
    
    
    %%
     fig = figure; set(fig,'Name','force');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            participation = selectionParticipate(selection);
            nt = data.trialNumber(selection);
            force = data.perf(selection);

        
        % compute stats
            nbin = 15;
            xx = tapply(nt,{nt,treatment,session},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
            yy =  tapply(force,{nt,treatment,session},@nanmean,{'continuous','discrete','discrete'},[nbin 2 2]);
            x = nanmean(xx,3);
            y = nanmean(yy,3);
            z = sem(yy,3);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
               cl = col{i};
               [h(i),hp] = boundedline( x(:,i) , y(:,i) , z(:,i) ,...
                                    'alpha','transparency', 0.5 ); 
                set(h(i),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
              
            end
            
            
            
        % legending
            legend([h(1) h(2)],{'placebo','atx'});
            ax = gca; 
%             ax.TickLength = [0 0];
%             ax.XTick = [];
            xlabel('session progression (total trial)');
            ylabel('force (%fmax)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
    %%
    
fig = figure; set(fig,'Name','force3');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        selection = (selectSubject & selectionChosen) ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        novelty = selectionRepetition(selection);

        
        x = [1,2];
        mu =  tapply(force,{novelty,treatment},@nanmean);
        err =   tapply(force,{novelty,treatment},@sem);
        
        for i = 1:2
            ind = [1:2];
            h(i) = errorbar( x ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = x;
        ax.XTickLabels = {'novel trial','repeated trial'};
        xlabel(texlabel('side'));
        ylabel('force peak (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end    
    
fig = figure; set(fig,'Name','force4');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        selection = (selectSubject & selectionChosen) ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        novelty = selectionRepetition(selection);

        
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(force,{effort,novelty,treatment},@nanmean);
        err =   tapply(force,{effort,novelty,treatment},@sem);
        
        for i = 1:2
            ind = [1:4];
            h(i) = errorbar( x(ind) ,mu(ind,i,1), err(ind,i,1) ,...
                     'Color',col{1},'LineWidth',2);
             if i==2;h(i).LineStyle='--';end
        end
        legend([h(1) h(2)],{'novel','repeated'});
        for i = 1:2
            ind = [1:4];
            h(i) = errorbar( x(ind) ,mu(ind,i,2), err(ind,i,2) ,...
                     'Color',col{2},'LineWidth',2);
             if i==2;h(i).LineStyle='--';end
        end

        ax = gca;
        ax.XTick = x;
%         ax.XTickLabels = {'novel trial','repeated trial'};
        xlabel(texlabel('effort required (%fmax)'));
        ylabel('force peak (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end    

    %%
fig = figure; set(fig,'Name','force5');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        selection = selectSubject & selectionChosen ;
        
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        reward = data.chosenCardinalReward(selection);
        side = data.sideChoice(selection);
        
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(force,{effort,side,treatment},@nanmean);
        err =   tapply(force,{effort,side,treatment},@sem);
        
        
        for s=1:2
            subplot(2,numel(subjectList),(s-1)*2+iSub);hold on
            for i = 1:2
                ind = [1:numel(unique(effort))];
                errbar( x(ind) ,mu(ind,s,i), err(ind,s,i) ,...
                         'Color',col{i},'LineWidth',2);
                h(i) = scatter(x(ind),mu(ind,s,i),'MarkerFaceColor',col{i},'MarkerEdgeColor',col{i},'LineWidth',1);
                plot(x(ind),mu(ind,s,i),'Color',col{i},'LineWidth',1);
            end
            plot(x(ind),x(ind),'--k','LineWidth',1);
            legend([h(1) h(2)],{'placebo','atx'});

            ax = gca;
            ax.XTick = x;
    %         ax.XLim = [0 1];
    %         ax.YLim = [0 1];
            xlabel(texlabel('effort required (%fmax)'));
            ylabel([sideList(s) ' force peak (%fmax)']); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
        end

    end   
    
%% 
fig = figure; set(fig,'Name','force25');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        selection = selectSubject & selectionChosen ;
        
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        reward = data.chosenCardinalReward(selection);
        
        x = tapply(effort,{effort,treatment},@nanmean);
        mu =  tapply(force,{effort,treatment},@nanmean);
        err =   tapply(force,{effort,treatment},@sem);
        
        
            subplot(1,numel(subjectList),iSub);hold on
            plot([0 1],[0 1],'--k','LineWidth',1);
            for i = 1:2
              [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
            end
            legend([h(1) h(2)],{'placebo','atx'});

            ax = gca;
            ax.XTick = x(:,1);
            ax.XLim = [0 1];
            ax.YLim = [0 1];
            xlabel(texlabel('effort required (%fmax)'));
            ylabel(['force peak (%fmax)']); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
        end

    
    
%%    
fig = figure; set(fig,'Name','force10');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        selection = selectSubject & selectionChosen ;
        
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        reward = data.chosenCardinalReward(selection);
        side = data.sideChoice(selection);
        
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(force,{effort,side,treatment},@nanmean);
        err =   tapply(force,{effort,side,treatment},@sem);
        
        
        for s = unique(side)'
            ie=0;
            for e = unique(effort)'
                ie = ie+1;
                subplot(4,numel(subjectList)*2,(ie-1)*4 +(iSub-1)*2 + s);hold on
                for i = 1:2
                    x = force(side==s & effort==e & treatment==trtList(i));
                    h(i) = histogram(x,[0:0.05:1],'Normalization','Probability',...
                        'FaceColor',col{i},'EdgeColor','none','FaceAlpha',0.4);
                    hold on;
                end
                legend([h(1) h(2)],{'placebo','atx'});

                ax = gca;
                ylabel(texlabel(['effort = ' num2str(round(e*100)) '(%fmax)']));
                xlabel([sideList(s)]); 
                yy =ylim;
                t = title(['monkey ' sub(1) ]); 
            end
        end

    end   
  
    
    
%%
fig = figure; set(fig,'Name','force11');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        selection = selectSubject & selectionChosen ;
%         selection = selectSubject & selectionChosen & ~selectionRepetition;

        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        correct = selectionCorrect(selection);
        side = data.sideChoice(selection);
        
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(correct,{effort,side,treatment},@nanmean);
        err =   tapply(correct,{effort,side,treatment},@sem);
        
        
        for s=1:2
            subplot(2,numel(subjectList),(s-1)*2+iSub);hold on
            for i = 1:2
                ind = [1:numel(unique(effort))];
                errbar( x(ind) ,mu(ind,s,i), err(ind,s,i) ,...
                         'Color',col{i},'LineWidth',2);
                h(i) = scatter(x(ind),mu(ind,s,i),'MarkerFaceColor',col{i},'MarkerEdgeColor',col{i},'LineWidth',1);
                plot(x(ind),mu(ind,s,i),'Color',col{i},'LineWidth',1);
            end
            legend([h(1) h(2)],{'placebo','atx'});

            ax = gca;
            ax.XTick = x;
    %         ax.XLim = [0 1];
    %         ax.YLim = [0 1];
            xlabel([ sideList(s) ' effort (%fmax)']);
            ylabel(['correct execution (%)']); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
        end

    end   
    
    
%%
fig = figure; set(fig,'Name','force12');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        selection = selectSubject & selectionChosen ;
%         selection = selectSubject & selectionChosen & ~selectionRepetition;

        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        error = abs(data.perf(selection)-data.chosenCardinalEffort(selection));
        side = data.sideChoice(selection);
        
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(error,{effort,side,treatment},@nanmean);
        err =   tapply(error,{effort,side,treatment},@sem);
        
        
        for s=1:2
            subplot(2,numel(subjectList),(s-1)*2+iSub);hold on
            for i = 1:2
                ind = [1:numel(unique(effort))];
                errbar( x(ind) ,mu(ind,s,i), err(ind,s,i) ,...
                         'Color',col{i},'LineWidth',2);
                h(i) = scatter(x(ind),mu(ind,s,i),'MarkerFaceColor',col{i},'MarkerEdgeColor',col{i},'LineWidth',1);
                plot(x(ind),mu(ind,s,i),'Color',col{i},'LineWidth',1);
            end
            legend([h(1) h(2)],{'placebo','atx'});

            ax = gca;
            ax.XTick = x;
    %         ax.XLim = [0 1];
    %         ax.YLim = [0 1];
            xlabel([ sideList(s) ' effort (%fmax)']);
            ylabel(['absolute error (%fmax)']); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
        end

    end   
    
%%
fig = figure; set(fig,'Name','force12');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
%         selection = selectSubject & selectionChosen ;
        selection = selectSubject & selectionChosen;

        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        reward =  data.chosenCardinalReward(selection);
        side = data.sideChoice(selection);
        
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(force,{effort,reward,treatment},@nanmean);
        err =   tapply(force,{effort,reward,treatment},@sem);
        
        
            for i = 1:2
                subplot(2,numel(subjectList),(i-1)*2+iSub);hold on

                for r=1:4
                    ind = [1:numel(unique(effort))];
                    errbar( x(ind) ,mu(ind,r,i), err(ind,r,i) ,...
                             'Color',col{i},'LineWidth',2);
                    scatter(x(ind),mu(ind,r,i),'MarkerFaceColor',col{i},'MarkerEdgeColor',col{i},'LineWidth',1);
                    h(r) = plot(x(ind),mu(ind,r,i),'Color',col{i},'LineWidth',r);
%                     if mod(r,2)==1;h(i).LineStyle='--';end
                end
                plot(x(ind),x(ind),'--k','LineWidth',0.5);
                legend([h(1) h(2) h(3) h(4)],{'R=1','R=2','R=3','R=4'});

            ax = gca;
            ax.XTick = x;
    %         ax.XLim = [0 1];
    %         ax.YLim = [0 1];
            xlabel([ ' effort (%fmax)']);
            ylabel(['force peak (%fmax)']); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            end

        end

    
%%
fig = figure; set(fig,'Name','force6');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        force = data.perf(selectSubject & selectionChosen &  ~selectionRepetition);
        effort = data.chosenCardinalEffort(selectSubject & selectionChosen & ~selectionRepetition);
        reward = data.chosenCardinalReward(selectSubject & selectionChosen & ~selectionRepetition);
        force = force - effort;
        
        
        x = tapply(reward,{reward},@nanmean);
        mu =  tapply(force,{reward,treatment},@nanmean);
        err =   tapply(force,{reward,treatment},@sem);
        
        for i = 1:2
            ind = [1:4];
            h(i) = errorbar( x ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = x;
%         ax.XTickLabels = {'left','right'};
        xlabel(texlabel('chosen reward '));
        ylabel('force exerted - required (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    
    
fig = figure; set(fig,'Name','force7');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
                selection = selectSubject & selectionChosen & ~selectionRepetition;
%         selection = selectSubject & selectionChosen ;
        
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        reward = data.chosenCardinalReward(selection);
        correct = selectionCorrect(selection);
        novelty = (~selectionRepetition(selection));

        
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(force,{effort,treatment},@nanmean);
        err =   tapply(force,{effort,treatment},@sem);
        
%         for i = 1:2
            ind = [1:numel(x)];
            plot(x,x,'--','Color',col{1},'LineWidth',1.5);
            h = plotSpread(force(treatment==trtList(1)),'distributionIdx',effort(treatment==trtList(1)),'distributionColors','k',...
                        'categoryIdx',treatment(treatment==trtList(1)),'categoryColors',col);
            h = plotSpread(force(treatment==trtList(2)),'distributionIdx',effort(treatment==trtList(2)),'distributionColors','k',...
                        'categoryIdx',treatment(treatment==trtList(2)),'categoryColors',col);
            h = errbar( x ,mu(ind), 0.1*ones(1,4) ,...
                     'Color',[0 0 0],'LineWidth',2,...
                     'horiz');
%         end

        ax = gca;
%         ax.XTick = x;
%         ax.XLim = x;
        ax.XLim = [0 1];
        ob = findobj('-property','MarkerSize');
        for io = 1:numel(ob); ob(io).MarkerSize = 10 ; end
        
%         ax.XTickLabels = {'left','right'};
        xlabel(texlabel('effort required (%fmax)'));
        ylabel('force peak (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end   
    
    
%% force distribution    
fig = figure; set(fig,'Name','force8');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        force = data.perf(selectSubject & selectionChosen &  ~selectionRepetition);
        effort = data.chosenCardinalEffort(selectSubject & selectionChosen & ~selectionRepetition);
        reward = data.chosenCardinalReward(selectSubject & selectionChosen & ~selectionRepetition);
        
        
        x = [1:10];
        binmax = 10;
        minF =  tapply(force,{force,treatment},@nanmin,{'continuous','discrete'},[binmax 2]);
        maxF =  tapply(force,{force,treatment},@nanmax,{'continuous','discrete'},[binmax 2]);

        mu =  tapply(force,{force,treatment},@nanmean,{'continuous','discrete'},[binmax 2]);
        err =   tapply(force,{force,treatment},@std,{'continuous','discrete'},[binmax 2]);
        
        for i = 1:2
%             ind = x;
%             h(i) = errorbar( x ,mu(ind,i), err(ind,i) ,...
%                      'Color',col{i},'LineWidth',2);
          for ind = x
                dist = force(force>=minF(ind,i) & force<maxF(ind,i) & treatment==trtList(i));
                hb = boxplot(dist,...
                    'boxstyle','filled','colors',col{i},...
                    'labels',ind,...
                    'positions',ind+ (i-1)./2,...
                    'medianstyle','line',...
                    'symbol','',...
                    'whisker',1,...
                    'widths',0.5);
           end
        end
%         legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XLim = [0 11];
        ax.YLim = [0 1.2];
        ax.XTick = x;
        ax.XTickLabels = num2cell(x);
        xlabel(texlabel('deciles of peak forces '));
        ylabel('peak force (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end   
%%
fig = figure; set(fig,'Name','force9');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        force = data.perf(selectSubject & selectionChosen &  ~selectionRepetition);
        effort = data.chosenCardinalEffort(selectSubject & selectionChosen & ~selectionRepetition);
        reward = data.chosenCardinalReward(selectSubject & selectionChosen & ~selectionRepetition);
        rt = data.reactionTime(selectSubject & selectionChosen &  ~selectionRepetition);

        
        x = [1:10];
        binmax = 10;
        minF =  tapply(force,{force,treatment},@nanmin,{'continuous','discrete'},[binmax 2]);
        maxF =  tapply(force,{force,treatment},@nanmax,{'continuous','discrete'},[binmax 2]);

        mu =  tapply(force,{force,treatment},@nanmean,{'continuous','discrete'},[binmax 2]);
        err =   tapply(force,{force,treatment},@nanstd,{'continuous','discrete'},[binmax 2]);
        
        for i = 1:2
            ind = x;
            h(i) = errorbar( mu(ind,i) ,err(ind,i), zeros(1,numel(ind)) ,...
                     'Color',col{i},'LineWidth',2,'Marker','o');
%           for ind = x
%                 dist = force(force>=minF(ind,i) & force<maxF(ind,i) & treatment==trtList(i));
%                 hb = boxplot(dist,...
%                     'boxstyle','filled','colors',col{i},...
%                     'labels',ind,...
%                     'positions',ind+ (i-1)./2,...
%                     'medianstyle','line',...
%                     'symbol','',...
%                     'whisker',1,...
%                     'widths',0.5);
%            end
        end
%         legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
%         ax.XLim = [0 11];
%         ax.YLim = [0 1.2];
        ax.XTick = [0:0.1:1];
%         ax.XTickLabels = num2cell( round(mean(mu,2),2) );
        xlabel(texlabel('mean peak forces (%fmax) '));
        ylabel('peak force variability (sem of %fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end   
 
     %%
    fig = figure; set(fig,'Name','force_error');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        side = data.sideChoice(selection);
        effort = data.chosenOrdinalEffort(selection);
        reward = data.chosenOrdinalReward(selection);
        rt = data.reactionTime(selection);
        
        error = abs(force-effort);
        for i = 1:2
            for e=1:4
                for r=1:4
                       sel = (treatment==trtList(i) & effort==e & reward==r);
                       sel = (treatment==trtList(i) & effort==e);
                       error(sel) = abs(errorscore(force(sel)));
                end
            end
        end
        
      % compute stats

%         x = tapply(force,{effort,reward,treatment},@nanmean);
%         y =  tapply(error,{effort,reward,treatment},@nanmean);
%         z =   tapply(error,{effort,reward,treatment},@sem);

        y =  tapply(error,{session,treatment,effort},@nanmean);
        y = nanmean(y,3);
        y =  tapply(error,{session,treatment},@nanmean);
        [h,p] = ttest2(exnan(y(:,1)),exnan(y(:,2)));
        z =   sem(y,1);
        y  =   nanmean(y,1);
        
%         y =  tapply(error,{treatment},@nanmean);
%         z =   tapply(error,{treatment},@sem);
%         
        % display
            subplot(numel(subjectList),1,iSub);hold on
            for i = 1:2
              [ b(i),~ ] = barplot( 1- 0.25 + 0.5*(i-1) ,y(i), z(i) , col{i} );
              b(i).BarWidth = 0.3;
            end
            if p<=0.05 ; s = sigstar({[0.75,1.25]},p); end
        
%         for i = 1:2
%             for r=1:4
%                   [~,~,h(i)] = errorscat(x(:,r,i),y(:,r,i), z(:,r,i)*0,col{i});
%             end
%             
%         end
%         for i =1:2
%             x = error(treatment==trtList(i));
%             h(i) = histogram(x,'BinWidth',0.02,'Normalization','Probability','FaceColor',col{i},'EdgeColor',col{i},'FaceAlpha',0.5);
%             xlabel('force residuals (%fmax)');
%             ylabel('frequency');
% 
%         end
        
%         legend([h(1) h(2)],{'placebo','atx'});
        legend([b(1) b(2)],{'placebo','atx'});

        ax = gca;
%         xlabel('force exerted (%fmax)'); 
        ylabel('force error (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
 
    %%
    fig = figure; set(fig,'Name','force_distribution');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        force = data.perf(selectSubject & selectionChosen &  ~selectionRepetition);
        effort = data.chosenOrdinalEffort(selectSubject & selectionChosen & ~selectionRepetition);

        
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(force,{effort,treatment},@nanmean);
        var =  tapply(force,{effort,treatment},@nanvar);
        err =   tapply(force,{effort,treatment},@sem);
        
        for i = 1:2
            subplot(2,numel(subjectList),iSub + 2*(i-1));hold on
            for ind = [1:4]
                dist = force(effort==ind & treatment==trtList(i));
                h(i) = histogram( dist,[0:0.01:1],...
                         'EdgeColor',col{i},'FaceColor','none','Normalization','countdensity');
            end
%             for ind = [1:4]
%                 dist = force(effort==ind & treatment==trtList(i));
%                 hb = boxplot(dist,...
%                     'boxstyle','filled','colors',col{i},...
%                     'labels',ind,...
%                     'positions',ind+ (i-1)./2,...
%                     'medianstyle','line',...
%                     'symbol','',...
%                     'whisker',1,...
%                     'widths',0.5);
%             end
        end
        legend([h(1) h(2)],{'placebo','atx'});
% 
%         for i = 1:2
%             ind = [1:4];
%             h(i) = errorbar( x ,var(ind,i), zeros(1,4) ,...
%                      'Color',col{i},'LineWidth',2);
%         end
%         legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = x;
%         ax.XLim = [0 5];

%         ax.XTickLabels = {'left','right'};
        xlabel(texlabel('effort required (%fmax)'));
        ylabel('force peak (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end   
    
  %%  
fig = figure; set(fig,'Name','dynamic');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
%         subplot(1,numel(subjectList),iSub);hold on
        
        y = []; indMax = 125;
        y = data.chosenPerfRT(selectSubject & selectionChosen  &  ~selectionRepetition,:);
        dd = nan(numel(y),1);
        for i=1:size(y,1)
           ind = min([indMax,numel(y{i,:})]);
           dd(i,1:ind) = y{i,:}(1:ind);
        end
%         data.chosenPerfRT(selectSubject & selectionChosen  &  ~selectionRepetition,:) = y;
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
%         dynamic = cell2mat(data.chosenPerfRT(selectSubject & selectionChosen  &  ~selectionRepetition,:));
        dynamic = dd;
        force = data.perf(selectSubject & selectionChosen &  ~selectionRepetition);
        effort = data.chosenOrdinalEffort(selectSubject & selectionChosen & ~selectionRepetition);
        reward = data.chosenOrdinalReward(selectSubject & selectionChosen & ~selectionRepetition);
        
        
        it=0; freq = forceSamplingFreq;
        for ie = 1:4
            for it = 1:2
                subplot(2,numel(subjectList),2*(iSub-1)+it);hold on

                
                sel = (treatment==trtList(it) & effort==ie);
                d = nanmean(dynamic(sel,:),1);
                e = sem(dynamic(sel,:),1);
                tmax = find(isnan(d),1,'first')-30;
                if isempty(tmax); tmax = numel(d);end;
                x = [1:numel(d(1:tmax))].*(1/freq);
                [h(it),hp] = boundedline( x , d(1:tmax) , e(1:tmax) , 'alpha' );   
                cl = col{it};
                set(h(it),'Color',cl,'LineWidth',2);set(hp,'FaceColor',cl);
                h(it).LineStyle = '-'; 

                ax = gca;
                ax.XTick = x;
                ax.XLim = [0 5];
        %         ax.XTickLabels = {'left','right'};
                xlabel(texlabel('time (sec) '));
                ylabel('force dynamic (%fmax)'); 
                yy =ylim;
                
                
            end
        end
        
%         legend([h(1) h(2) h(3) h(4)],'10%','40%','70%','100%');
        

        ax = gca;
        ax.XTick = x;
        ax.XLim = [0 5];
%         ax.XTickLabels = {'left','right'};
        xlabel(texlabel('time (sec) '));
        ylabel('force dynamic (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end   
    
    
    
%% 5) choice
   fig = figure; set(fig,'Name','choiceR');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen & data.deltaOrdinalReward~=0 ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            choice = (data.sideChoice(selection)-1)*2-1;
            dR = round(data.deltaOrdinalReward(selection));
            choiceR = (sign(choice)==sign(dR));
        
        % compute stats
            mu =  tapply(choiceR,{session,treatment},@nanmean);
            [h,p] = ttest2(exnan(mu(:,1)),exnan(mu(:,2)));
            err =   sem(mu,1);
            mu  =   nanmean(mu,1);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [ b(i),~ ] = barplot( 1- 0.25 + 0.5*(i-1) ,mu(i), err(i) , col{i} );
              b(i).BarWidth = 0.3;
            end
            if p<=0.05 ; s = sigstar({[0.75,1.25]},p); end
            
            
            
        % legending
            legend([b(1) b(2)],{'placebo','atx'});
            ax = gca; 
            ax.YLim = [0.5 1];
            ax.TickLength = [0 0];
            ax.XTick = [];

            ylabel('choice = max(reward) (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
%%
fig = figure; set(fig,'Name','choiceE');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen & data.deltaOrdinalEffort~=0 ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            choice = (data.sideChoice(selection)-1)*2-1;
            dE = round(data.deltaOrdinalEffort(selection));
            choiceE = (sign(choice)~=sign(dE));
        
        % compute stats
            mu =  tapply(choiceE,{session,treatment},@nanmean);
            [h,p] = ttest2(exnan(mu(:,1)),exnan(mu(:,2)));
            err =   sem(mu,1);
            mu  =   nanmean(mu,1);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [ b(i),~ ] = barplot( 1- 0.25 + 0.5*(i-1) ,mu(i), err(i) , col{i} );
              b(i).BarWidth = 0.3;
            end
            if p<=0.05 ; s = sigstar({[0.75,1.25]},p); end
            
            
            
        % legending
            legend([b(1) b(2)],{'placebo','atx'});
            ax = gca; 
            ax.YLim = [0.5 1];
            ax.TickLength = [0 0];
            ax.XTick = [];

            ylabel('choice = min(effort) (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
%%

fig = figure; set(fig,'Name','choice');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        force = data.perf(selectSubject & selectionChosen &  ~selectionRepetition);
        choice = data.sideChoice(selectSubject & selectionChosen & ~selectionRepetition)-1;
        dR = (data.deltaCardinalReward(selectSubject & selectionChosen & ~selectionRepetition));
        dR = round(dR,3);
        
        x = tapply(dR,{dR},@nanmean);
        mu =  tapply(choice,{dR,treatment},@nanmean);
        err =   tapply(choice,{dR,treatment},@sem);
        
        for i = 1:2
            ind = [1:numel(unique(dR))];
            h(i) = errorbar( x ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel(texlabel('Delta reward'));
        ylabel('right choice (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end

%%
fig = figure; set(fig,'Name','choicedE');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            choice = (data.sideChoice(selection)-1);
            dE = round(data.deltaOrdinalEffort(selection));
        
        % compute stats
            x = tapply(dE,{dE,treatment},@nanmean);
            y =  tapply(choice,{dE,treatment},@nanmean);
            z =   tapply(choice,{dE,treatment},@sem);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            xl = [min(dE)-1 max(dE)+1];
            plot([1 1]*0,[0 1],'k--','LineWidth',0.5);
            plot(xl,[1 1]*0.5,'k--','LineWidth',0.5);
            for i = 1:2
              [~,~,h(i)] = errorscat(x(:,i) ,y(:,i), z(:,i),col{i});
            end

            
            
        % legending
            legend([h(1) h(2)],{'placebo','atx'});
            ax = gca; 
%             ax.TickLength = [0 0];
%             ax.XTick = [];
            ax.XLim = xl;
            xlabel(texlabel('E_right - E_left')); 
            ylabel('right choice (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
%%
    fig = figure; set(fig,'Name','choicedR');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            choice = (data.sideChoice(selection)-1);
            dR = round(data.deltaOrdinalReward(selection));
        
        % compute stats
            x = tapply(dR,{dR,treatment},@nanmean);
            y =  tapply(choice,{dR,treatment},@nanmean);
            z =   tapply(choice,{dR,treatment},@sem);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            xl = [min(dR)-1 max(dR)+1];
            plot([1 1]*0,[0 1],'k--','LineWidth',0.5);
            plot(xl,[1 1]*0.5,'k--','LineWidth',0.5);
            for i = 1:2
              [~,~,h(i)] = errorscat(x(:,i) ,y(:,i), z(:,i),col{i});
            end

            
            
        % legending
            legend([h(1) h(2)],{'placebo','atx'});
            ax = gca; 
%             ax.TickLength = [0 0];
%             ax.XTick = [];
            ax.XLim = xl;
            xlabel(texlabel('R_right - R_left')); 
            ylabel('right choice (%)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
  
%%
fig = figure; set(fig,'Name','choice2');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen;
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        choice = data.sideChoice(selection)-1;
        dE = (data.deltaCardinalEffort(selection));
        dE = round(dE,3);
        dR = (data.deltaCardinalReward(selection));
        dR = round(dR,3);
        dE3 = data.cardinalEffortRight(selection).*data.calibRight(selection) - ...
                data.cardinalEffortLeft(selection).*data.calibLeft(selection);
        E2 = data.cardinalEffortRight(selection);
        E1 = data.cardinalEffortLeft(selection);
        
        % fit
        yy = nan(numel(treatment),1);
        dE2 = nan(numel(treatment),1);
        for i=1:2
            sel = treatment==trtList(i);
            [beta,~,stat] = glmfit([dR(sel),E2(sel), E1(sel)],choice(sel),'binomial','link','logit','constant','off');
            yy(sel) = glmval(beta,[dR(sel),E2(sel), E1(sel)],'logit','Constant','off');
%             dE2(sel) = -beta(3).*([E2(sel), E1(sel)]);
        end
        
%         dE3 = nan(numel(treatment),1);
%         [beta,~,stat] = glmfit([r1(sel),E(sel)],choice(sel),'binomial','link','logit');
%         yy(sel) = glmval(beta,[dR(sel),dE(sel)],'logit');
%         dE2(sel) = -beta(3).*dE(sel);
        

        delta = dE;
        x = tapply(delta,{dE,treatment},@nanmean);
        mu =  tapply(choice,{dE,treatment},@nanmean);
        err =   tapply(choice,{dE,treatment},@sem);
        mu2 =   tapply(yy,{dE,treatment},@nanmean);

        for i = 1:2
            
            ind = [1:numel(unique(dE))];
            errbar( x(ind,i) ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',1);
            h(i) = scatter(x(ind,i),mu(ind,i),'MarkerFaceColor',col{i},'MarkerEdgeColor',col{i},'LineWidth',1);
            plot(x(ind,i),mu(ind,i),'Color',col{i},'LineWidth',1);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel(texlabel('Delta effort'));
        ylabel('right choice (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
   
%%
fig = figure; set(fig,'Name','choice5');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        choice = data.sideChoice(selection)-1;
        dR = data.deltaCardinalReward(selection);
        dE = data.deltaCardinalEffort(selection);
        dE = round(dE,3);
        dR = round(dR,3);
        congruent = ((dR.*(-dE))>0)*2;
        congruent((dE.*dR)==0) = 1;
        E2 = data.cardinalEffortRight(selection);
        E1 = data.cardinalEffortLeft(selection);
        
         % fit
             yy = nan(numel(treatment),1);
             predictor = table(rangescore(dR),rangescore(E2),rangescore(E1),(treatment) ,(choice) );
             predictor.Properties.VariableNames = {'dR';'E2';'E1';'atx';'choice'};
             formula = 'choice ~ -1 + dR + E2 + E1 + atx:(E2 + E1) ';
             glm = fitglm( predictor, formula,'Distribution','binomial','link','logit','CategoricalVars',[4]);
             yy = glm.Fitted.Response;
        
        
        x = tapply(dR,{dR},@nanmean);
        mu =  tapply(choice,{dR,congruent,treatment},@nanmean);
        mu2 =  tapply(yy,{dR,congruent,treatment},@nanmean);
        err =   tapply(choice,{dR,congruent,treatment},@sem);
        congList = {'conflicting','unidimensional','congruent'};
        
        for ig=1:3
            
            subplot(numel(subjectList),3,3*(2-iSub)+ig);hold on

            for i = 1:2
                ind = [1:numel(unique(dR))];
                errbar( x ,mu(ind,ig,i), err(ind,ig,i) ,...
                         'Color',col{i},'LineWidth',1);
                h(i) = scatter(x,mu(ind,ig,i),'MarkerFaceColor',col{i},'MarkerEdgeColor',col{i},'LineWidth',1);
                plot(x,mu2(ind,ig,i),'Color',col{i},'LineWidth',1);
            end
            
            ax = gca;
            xlabel(texlabel('Delta reward'));
            if ig==1  
                ylabel(['monkey ' sub(1) ]); 
            else
                
                ylabel('right choice (%)'); 
                if ig==3;legend([h(1) h(2)],{'placebo','atx'});end
            end
            title(congList(ig));
            yy =ylim;

        end


    end

    
fig = figure; set(fig,'Name','choice6');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
%         subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen & ~selectionRepetition;
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        choice = data.sideChoice(selection)-1;
        dR = data.deltaCardinalReward(selection);
        dE = data.deltaCardinalEffort(selection);
        dE = round(dE,3);
        dR = round(dR,3);
        congruent = ((dR.*(-dE))>0)*2;
        congruent((dE.*dR)==0) = 1;
        E2 = data.cardinalEffortRight(selection);
        E1 = data.cardinalEffortLeft(selection);
        
         % fit
             yy = nan(numel(treatment),1);
             predictor = table(rangescore(dR),rangescore(E2),rangescore(E1),(treatment) ,(choice) );
             predictor.Properties.VariableNames = {'dR';'E2';'E1';'atx';'choice'};
             formula = 'choice ~ -1 + dR + E2 + E1 + atx:(E2 + E1) ';
             glm = fitglm( predictor, formula,'Distribution','binomial','link','logit','CategoricalVars',[4]);
             yy = glm.Fitted.Response;
        
        x = tapply(dE,{dE},@nanmean);
        mu =  tapply(choice,{dE,congruent,treatment},@nanmean);
        mu2 =  tapply(yy,{dE,congruent,treatment},@nanmean);
        err =   tapply(choice,{dE,congruent,treatment},@sem);
        congList = {'conflicting','unidimensional','congruent'};
        
        for ig=1:3
            
            subplot(numel(subjectList),3,3*(2-iSub)+ig);hold on

            for i = 1:2
                ind = [1:numel(unique(dE))];
                errbar( x ,mu(ind,ig,i), err(ind,ig,i) ,...
                         'Color',col{i},'LineWidth',1);
                h(i) = scatter(x,mu(ind,ig,i),'MarkerFaceColor',col{i},'MarkerEdgeColor',col{i},'LineWidth',1);
                plot(x,mu2(ind,ig,i),'Color',col{i},'LineWidth',1);
            end
            
            ax = gca;
            xlabel(texlabel('Delta effort'));
            if ig==1  
                ylabel(['monkey ' sub(1) ]); 
            else
                
                ylabel('right choice (%)'); 
                if ig==3;legend([h(1) h(2)],{'placebo','atx'});end
            end
            title(congList(ig));
            yy =ylim;
            
        end

    end

    
%%
fig = figure; set(fig,'Name','choice3');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        force = data.perf(selectSubject & selectionChosen &  ~selectionRepetition);
        choice = data.sideChoice(selectSubject & selectionChosen & ~selectionRepetition)-1;
        dE = (data.deltaCardinalEffort(selectSubject & selectionChosen & ~selectionRepetition));
        dE = round(dE,3);
        e1 = (data.ordinalEffortLeft(selectSubject & selectionChosen & ~selectionRepetition));
        e2 = (data.ordinalEffortRight(selectSubject & selectionChosen & ~selectionRepetition));

        
        x1 = tapply(e1,{e1},@nanmean);
        x2 = tapply(e2,{e2},@nanmean);
        mu =  tapply(choice,{e1,e2,treatment},@nanmean);
        
        for i = 1:2
            
            subplot(2,numel(subjectList), iSub + (i-1)*2);hold on
            
            p = mu(:,:,i);
            cm = colormap('cool'); cm = flipud(cm);
            minCM = min(min(p))-0.3 ; maxCM = max(max(p))+0.3;
            minCM = min(min(p))*(1) ; maxCM = max(max(p))*(1);
            minCM = 0 ; maxCM = 1;

            heatmap(p,x2,x1,(p),'Colormap',cm,'NaNColor',[1 1 1],...
            'MaxColorValue',maxCM,'MinColorValue',minCM,'TextColor',[1 1 1]*0);
            colorbar;

            ax = gca;
            ax.XTick = [x2];
            ax.XTickLabel = num2cell(x2);
            ax.YTick = [x1];
            ax.YTickLabel = num2cell(x1);
            xlabel(texlabel(' effort (right)'));
            ylabel(texlabel(' effort (left)'));
            yy =ylim;
            if i==1; t = title(['monkey ' sub(1) ]); end
            if i==2; t = title(['right choice (%)' ]); end

        
        end
    end
    
    
fig = figure; set(fig,'Name','choice4');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        
        session = data.session(selectSubject & selectionChosen & ~selectionRepetition);
        treatment = data.treatment(selectSubject & selectionChosen & ~selectionRepetition);
        force = data.perf(selectSubject & selectionChosen &  ~selectionRepetition);
        choice = data.sideChoice(selectSubject & selectionChosen & ~selectionRepetition)-1;
        dE = (data.deltaCardinalEffort(selectSubject & selectionChosen & ~selectionRepetition));
        dE = round(dE,3);
        r1 = (data.ordinalRewardLeft(selectSubject & selectionChosen & ~selectionRepetition));
        r2 = (data.ordinalRewardRight(selectSubject & selectionChosen & ~selectionRepetition));

        
        x1 = tapply(r1,{r1},@nanmean);
        x2 = tapply(r2,{r2},@nanmean);
        mu =  tapply(choice,{r1,r2,treatment},@nanmean);
        
        for i = 1:2
            
            subplot(2,numel(subjectList), iSub + (i-1)*2);hold on
            
            p = mu(:,:,i);
            cm = colormap('cool'); cm = flipud(cm);
            minCM = min(min(p))-0.3 ; maxCM = max(max(p))+0.3;
            minCM = min(min(p))*(1) ; maxCM = max(max(p))*(1);
            minCM = 0 ; maxCM = 1;

            heatmap(p,x2,x1,(p),'Colormap',cm,'NaNColor',[1 1 1],...
            'MaxColorValue',maxCM,'MinColorValue',minCM,'TextColor',[1 1 1]*0);
            colorbar;

            ax = gca;
            ax.XTick = [x2];
            ax.XTickLabel = num2cell(x2);
            ax.YTick = [x1];
            ax.YTickLabel = num2cell(x1);
            xlabel(texlabel(' reward (right)'));
            ylabel(texlabel(' reward (left)'));
            yy =ylim;
            if i==1; t = title(['monkey ' sub(1) ]); end
            if i==2; t = title(['right choice (%)' ]); end

        
        end
    end

%%
    fig = figure; set(fig,'Name','choice8');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        selection = selectSubject & selectionChosen; 
        session = data.session(selection);
        treatment = data.treatment(selection);
        choice = data.sideChoice(selection)-1;
        dE = (data.deltaCardinalEffort(selection));
        dE = round(dE,3);
        dR = (data.deltaCardinalReward(selection));
        dR = round(dR,3);

        
        x1 = tapply(dR,{dR},@nanmean);
        x2 = tapply(dE,{dE},@nanmean);
        mu =  tapply(choice,{dR,dE,treatment},@nanmean);
        
        for i = 1:2
            
            subplot(2,numel(subjectList), iSub + (i-1)*2);hold on
            
            p = mu(:,:,i);
            cm = colormap('cool'); cm = flipud(cm);
            minCM = min(min(p))-0.3 ; maxCM = max(max(p))+0.3;
            minCM = min(min(p))*(1) ; maxCM = max(max(p))*(1);
            minCM = 0 ; maxCM = 1;

            heatmap(p,x2,x1,(p),'Colormap',cm,'NaNColor',[1 1 1],...
            'MaxColorValue',maxCM,'MinColorValue',minCM,'TextColor',[1 1 1]*0);
            colorbar;

            ax = gca;
            ax.XTick = [x2];
            ax.XTickLabel = num2cell(x2);
            ax.YTick = [x1];
            ax.YTickLabel = num2cell(x1);
            xlabel(texlabel(' Delta E'));
            ylabel(texlabel(' Delta R'));
            yy =ylim;
            if i==1; t = title(['monkey ' sub(1) ]); end
            if i==2; t = title(['right choice (%)' ]); end

        
        end
    end
    
%%
    fig = figure; set(fig,'Name','choice9'); 
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        subplot(1,numel(subjectList),iSub);hold on

        selection = selectSubject & selectionChosen; 
        session = data.session(selection);
        treatment = data.treatment(selection);
        choice = data.sideChoice(selection)-1;
        dE = (data.deltaCardinalEffort(selection));
        dE = round(dE,3);
        sumE = (data.cardinalEffortLeft(selection) + data.cardinalEffortRight(selection));
        sumE = round(sumE,3);

        
        x1 = tapply(dE,{dE,sumE},@nanmean);
        x2 = tapply(sumE,{dE,sumE},@nanmean);
        mu =  tapply(choice,{dE,sumE,treatment},@nanmean);
        err =  tapply(choice,{dE,sumE,treatment},@sem);
        
        for i = 1:2
          x = [1:4];
          ind = [1:2:7];
          b(i) = plot( x , mu(4,ind,i),'Color',col{i},'LineWidth',2) ;
          h = errbar( x ,mu(4,ind,i), err(4,ind,i) ,...
                 'Color',col{i},'LineWidth',2);
             
        end

        legend([b(1) b(2)],{'placebo','atx'});
        ax = gca; 
        xlabel('effort levels (E1 & E2)'); 
        ylabel('right choice | dE=0 (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
    fig = figure; set(fig,'Name','choice10'); 
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        

        selection = selectSubject & selectionChosen; 
        session = data.session(selection);
        treatment = data.treatment(selection);
        choice = data.sideChoice(selection)-1;
        dR = (data.deltaCardinalReward(selection));
        dR = round(dR,3);
        dE = (data.deltaCardinalEffort(selection));
        dE = round(dE,3);
        sumE = (data.cardinalEffortLeft(selection) + data.cardinalEffortRight(selection));
        sumE = round(sumE,3);

        x = tapply(dR,{dR},@nanmean);
        x1 = tapply(dR,{dR,dE,sumE},@nanmean);
        x2 = tapply(dE,{dR,dE,sumE},@nanmean);
        x3 = tapply(sumE,{dR,dE,sumE},@nanmean);
        mu =  tapply(choice,{dR,dE,sumE,treatment},@nanmean);
        err =  tapply(choice,{dR,dE,sumE,treatment},@sem);
        
        for i = 1:2
            
            subplot(2,numel(subjectList),(i-1)*2+iSub);hold on

            for ind = [1:2:7]
                  b(i) = plot( x , mu(:,4,ind,i),'Color',col{i},'LineWidth',ind/2) ;
                  h = errbar( x ,mu(:,4,ind,i), err(:,4,ind,i) ,...
                         'Color',col{i},'LineWidth',ind/2);
            end

        ax = gca; 
        xlabel('Delta R'); 
        ylabel('right choice | E1&E2 (%)'); 
        yy =ylim;
        t = title(['unidimensional: monkey ' sub(1) ]); 
        end
            legend([b(1) b(2)],{'placebo','atx'});

    end
%%    
fig = figure; set(fig,'Name','choice11'); 
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
            subplot(1,numel(subjectList),iSub);hold on

        selection = selectSubject & selectionChosen; 
        session = data.session(selection);
        treatment = data.treatment(selection);
        choice = data.sideChoice(selection)-1;
        r1 = (data.ordinalRewardLeft(selection));
        r2 = (data.ordinalRewardRight(selection));

        x1 = tapply(r1,{r1,r2},@nanmean);
        x2 = tapply(r2,{r1,r2},@nanmean);
        mu =  tapply(choice,{r1,r2,treatment},@nanmean);
        err =  tapply(choice,{r1,r2,treatment},@sem);
        
        for i = 1:2
          x = [1:12];
          y =[mu(4,1,i) mu(4,2,i) mu(3,1,i) mu(4,3,i) mu(3,2,i) mu(2,1,i),...
              mu(1,2,i) mu(2,3,i) mu(3,4,i) mu(1,3,i) mu(2,4,i) mu(1,4,i)];
          z =[err(4,1,i) err(4,2,i) err(3,1,i) err(4,3,i) err(3,2,i) err(2,1,i),...
              err(1,2,i) err(2,3,i) err(3,4,i) err(1,3,i) err(2,4,i) err(1,4,i)];

          b(i) = plot( x , y,'Color',col{i},'LineWidth',2) ;
          h = errbar( x ,y, z ,...
                 'Color',col{i},'LineWidth',2);

            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = {'4/1','4/2','3/1','4/3','3/2','2/1',...
                             '1/2','2/3','3/4','1/3','2/4','1/4'};
            xlabel('reward options (left | right) '); 
            ylabel('right choice  (%)'); 
            yy =ylim;
            t = title([' monkey ' sub(1) ]); 
        end
        legend([b(1) b(2)],{'placebo','atx'});

    end
    %%    
fig = figure; set(fig,'Name','choice12'); 
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
            subplot(1,numel(subjectList),iSub);hold on

        selection = selectSubject & selectionChosen; 
        session = data.session(selection);
        treatment = data.treatment(selection);
        choice = data.sideChoice(selection)-1;
        e1 = (data.ordinalEffortLeft(selection));
        e2 = (data.ordinalEffortRight(selection));

        x1 = tapply(e1,{e1,e2},@nanmean);
        x2 = tapply(e2,{e1,e2},@nanmean);
        mu =  tapply(choice,{e1,e2,treatment},@nanmean);
        err =  tapply(choice,{e1,e2,treatment},@sem);
        
        for i = 1:2
          x = [1:12];
          y =[mu(4,1,i) mu(4,2,i) mu(3,1,i) mu(4,3,i) mu(3,2,i) mu(2,1,i),...
              mu(1,2,i) mu(2,3,i) mu(3,4,i) mu(1,3,i) mu(2,4,i) mu(1,4,i)];
          z =[err(4,1,i) err(4,2,i) err(3,1,i) err(4,3,i) err(3,2,i) err(2,1,i),...
              err(1,2,i) err(2,3,i) err(3,4,i) err(1,3,i) err(2,4,i) err(1,4,i)];

          b(i) = plot( x , y,'Color',col{i},'LineWidth',2) ;
          h = errbar( x ,y, z ,...
                 'Color',col{i},'LineWidth',2);

            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = {'4/1','4/2','3/1','4/3','3/2','2/1',...
                             '1/2','2/3','3/4','1/3','2/4','1/4'};
            xlabel('effort options (left | right) '); 
            ylabel('right choice  (%)'); 
            yy =ylim;
            t = title([' monkey ' sub(1) ]); 
        end
        legend([b(1) b(2)],{'placebo','atx'});

    end
%% 6) response time
fig = figure; set(fig,'Name','rt');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            rt = data.responseTime(selection);
            value = data.participationValue(selection);
        
        % compute stats
            nbin = [6 2 2];
%             mu =  tapply(rt,{value,session,treatment},@nanmean,{'continuous','discrete','discrete'},nbin);
%             mu = nanmean(mu,1);
%             [h,p] = ttest2(exnan(mu(:,:,1)),exnan(mu(:,:,2)));
%             err =   sem(mu,2);
%             mu  =   nanmean(mu,2);
            
            mu =  tapply(rt,{session,treatment},@nanmean);
            [h,p] = ttest2(exnan(mu(:,1)),exnan(mu(:,2)));
            err =   sem(mu,1);
            mu  =   nanmean(mu,1);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [ b(i),~ ] = barplot( 1- 0.25 + 0.5*(i-1) ,mu(i), err(i) , col{i} );
              b(i).BarWidth = 0.3;
            end
            if p<=0.05 ; s = sigstar({[0.75,1.25]},p); end
            
            
            
        % legending
            legend([b(1) b(2)],{'placebo','atx'});
            ax = gca; 
%             ax.YLim = [0.5 1];
            ax.TickLength = [0 0];
            ax.XTick = [];

            ylabel('response time (s.)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
    
%%
fig = figure; set(fig,'Name','rt_var');
    for iSub = 1:numel(subjectList)
        % select
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            selection = selectSubject & selectionChosen ;
       
        % variables
            session = data.session(selection);
            treatment = data.treatment(selection);
            rt = data.responseTime(selection);
        
        % compute stats
            mu =  tapply(rt,{session,treatment},@nanvar);

            [h,p] = ttest2(exnan(mu(:,1)),exnan(mu(:,2)));
            err =   sem(mu,1);
            mu  =   nanmean(mu,1);

        % display
            subplot(1,numel(subjectList),iSub);hold on
            for i = 1:2
              [ b(i),~ ] = barplot( 1- 0.25 + 0.5*(i-1) ,mu(i), err(i) , col{i} );
              b(i).BarWidth = 0.3;
            end
            if p<=0.05 ; s = sigstar({[0.75,1.25]},p); end
            
            
            
        % legending
            legend([b(1) b(2)],{'placebo','atx'});
            ax = gca; 
%             ax.YLim = [0.5 1];
            ax.TickLength = [0 0];
            ax.XTick = [];

            ylabel('response time variance (s.)'); 
            yy =ylim;
            t = title(['monkey ' sub(1) ]); 
            
    end
    
%%
fig = figure; set(fig,'Name','rt');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        
        mu =  tapply(rt,{treatment},@nanmean);
        err =   tapply(rt,{treatment},@sem);
        
        
        for i = 1:2
          b(i) = bar( 1- 0.25 + 0.5*(i-1) , mu(i)) ;
          h = errbar( 1- 0.25 + 0.5*(i-1) ,mu(i), err(i) ,...
                 'Color',[0.5 0.5 0.5],'LineWidth',2);
             b(i).FaceColor = col{i};
             b(i).EdgeColor = 'none';
             b(i).BarWidth = 0.3;
             
        end

        legend([b(1) b(2)],{'placebo','atx'});
        ax = gca; 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    
fig = figure; set(fig,'Name','rt2');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
  
        mu =  tapply(rt,{session},@nanmean);
        err =   tapply(rt,{session},@sem);
        
        for i = 1:2
            ind = trtSessions{i};
            h(i) = errorbar( ind ,mu(ind), err(ind) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = [1:numel(sessionList)];
        ax.XTickLabels = cellstr(sessionList');
        ax.XTickLabelRotation = 45;
        xlabel('session date');
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    
%%
fig = figure; set(fig,'Name','rt3');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        dR = data.deltaCardinalReward(selection);
        dE = data.deltaCardinalEffort(selection);
        dE = round(dE,3);
        dR = round(dR,3);
  
        x = tapply(dR,{dR},@nanmean);
        mu =  tapply(rt,{dR,treatment},@nanmean);
        err =   tapply(rt,{dR,treatment},@sem);
        
       for i = 1:2
            ind = [1:numel(unique(dR))];
            h(i) = errorbar( x ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel(texlabel('Delta reward'));
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end

fig = figure; set(fig,'Name','rt4');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        dR = data.deltaCardinalReward(selection);
        dE = data.deltaCardinalEffort(selection);
        dE = round(dE,3);
        dR = round(dR,3);
        
          x = tapply(dE,{dE},@nanmean);
        mu =  tapply(rt,{dE,treatment},@nanmean);
        err =   tapply(rt,{dE,treatment},@sem);
        
       for i = 1:2
            ind = [1:numel(unique(dE))];
            h(i) = errorbar( x ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel(texlabel('Delta effort'));
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
fig = figure; set(fig,'Name','rt5');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
%         selection = selectSubject & selectionChosen & ~selectionRepetition;
        session = data.session(selection);
        dR = data.deltaCardinalReward(selection);
        dE = data.deltaCardinalEffort(selection);
        dE = round(dE,3);
        dR = round(dR,3);
        congruent = ((dR.*(-dE))>0)*2;
        congruent((dE.*dR)==0) = 1;

        treatment = data.treatment(selection);
%         rt = data.reactionTime(selection)-data.offerTime(selection);        
        rt = data.responseTime(selection);        
  
        mu =  tapply(rt,{congruent,treatment},@nanmean);
        err =   tapply(rt,{congruent,treatment},@sem);
        
       for i = 1:2
            ind = [1:3];
            x= ind;
            [~,~,h(i)] = errorscat(x ,mu(ind,i), err(ind,i),col{i});

        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = [1:3];
        ax.XTickLabel = {'conflicting','unidimensional','congruent'};
        xlabel('dimension congruence'); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    %%
    fig = figure; set(fig,'Name','rt5');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
%         selection = selectSubject & selectionChosen & ~selectionRepetition;
        session = data.session(selection);
        choice = data.sideChoice(selection)-1;
        r1 = data.ordinalRewardLeft(selection);
        r2 = data.ordinalRewardRight(selection);
        e1 = data.ordinalEffortLeft(selection);
        e2 = data.ordinalEffortRight(selection);
        nt = data.trialNumber(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        
         dv = nan(numel(rt),1);
         conflict = nan(numel(rt),1);
         % fit choice
         for i = 1:2
             sel = (treatment==trtList(i) );
             predictor = table(rangescore(r2(sel)),rangescore(r1(sel)),rangescore(e2(sel)),rangescore(e1(sel)),rangescore(nt(sel)),choice(sel) );
             predictor.Properties.VariableNames = {'R2';'R1';'E2';'E1';'ntrial';'choice'};
             formula = 'choice ~ -1 + R1 + R2 + E1 + E2';
             glm = fitglm( predictor,formula,'Distribution','binomial','link','logit');
             dv(sel) = glm.Fitted.LinearPredictor;
             beta = glm.Coefficients.Estimate;
             desmat = [rangescore(r2(sel)),rangescore(r1(sel)),rangescore(e2(sel)),rangescore(e1(sel))];
             dvR = (((beta(1:2))')*(desmat(:,1:2)'))';
             dvE = (((beta(3:4))')*(desmat(:,3:4)'))';
             vR = ((abs(beta(1:2))')*(desmat(:,1:2)'))';
             vE = ((abs(beta(3:4))')*(desmat(:,3:4)'))';
             conflictR = abs(dvR);
             conflictE = abs(dvE);
             
             conflict(sel) = abs(dv(sel));
%              conflict(sel) = conflictR;
%              conflict(sel) = conflictE;
%              conflict(sel) = conflictR+conflictE;
%              conflict(sel) = dvR.*(-dvE);
%              conflict(sel) = abs(vR-vE);
%              conflict(sel) = vR+vE;


         end
  
        nbin = 10;
        x = tapply(conflict,{conflict,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        mu =  tapply(rt,{conflict,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        err =   tapply(rt,{conflict,treatment},@sem,{'continuous','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('conflict'); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    %%
    fig = figure; set(fig,'Name','rt8');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
%         selection = selectSubject & selectionChosen & ~selectionRepetition;
        session = data.session(selection);
        choice = data.sideChoice(selection)-1;
        e = data.chosenCardinalEffort(selection);
        r = data.chosenCardinalReward(selection);
        nt = data.trialNumber(selection);
        treatment = data.treatment(selection);
%         rt = data.reactionTime(selection)-data.offerTime(selection);
        rt = data.responseTime(selection);

  
        nbin = 10;
        x = tapply(r,{r,treatment},@nanmean);
        mu =  tapply(rt,{r,treatment},@nanmean);
        err =   tapply(rt,{r,treatment},@sem);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('chosen reward'); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
    fig = figure; set(fig,'Name','rt9');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
%         selection = selectSubject & selectionChosen & ~selectionRepetition;
        session = data.session(selection);
        choice = data.sideChoice(selection)-1;
        e = data.chosenCardinalEffort(selection);
        r = data.chosenCardinalReward(selection);
        nt = data.trialNumber(selection);
        treatment = data.treatment(selection);
%         rt = data.reactionTime(selection)-data.offerTime(selection);
        rt = data.responseTime(selection);        
  
        nbin = 10;
        x = tapply(e,{e,treatment},@nanmean);
        mu =  tapply(rt,{e,treatment},@nanmean);
        err =   tapply(rt,{e,treatment},@sem);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('chosen effort'); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
%%
    fig = figure; set(fig,'Name','rt10');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
%         selection = selectSubject & selectionChosen & ~selectionRepetition;
        session = data.session(selection);
        choice = data.sideChoice(selection)-1;
        force = data.perf(selection);
        r = data.chosenCardinalReward(selection);
        nt = data.trialNumber(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        
  
        nbin = 10;
        x = tapply(force,{force,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        mu =  tapply(rt,{force,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        err =   tapply(rt,{force,treatment},@sem,{'continuous','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('exerted force (%fmax)'); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
             
    
%%
fig = figure; set(fig,'Name','rt11');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        sv = data.participationValue(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        
  
        nbin = 6;
        x = tapply(sv,{sv,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        mu =  tapply(rt,{sv,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        err =   tapply(rt,{sv,treatment},@sem,{'continuous','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('offer value '); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    
%%
fig = figure; set(fig,'Name','rt12');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        v = round(data.ordinalRewardLeft(selection)+data.ordinalRewardRight(selection));
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        
  
        nbin = 10;
        x = tapply(v,{v,treatment},@nanmean,{'discrete','discrete'},[nbin 2]);
        mu =  tapply(rt,{v,treatment},@nanmean,{'discrete','discrete'},[nbin 2]);
        err =   tapply(rt,{v,treatment},@sem,{'discrete','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('R1 + R2'); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    
%%
fig = figure; set(fig,'Name','rt13');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        v = round(data.ordinalEffortLeft(selection)+data.ordinalEffortRight(selection));
%         v = round(data.ordinalEffortLeft(selection));
%         v = round(data.ordinalEffortRight(selection));

        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        
  
        nbin = 10;
        x = tapply(v,{v,treatment},@nanmean,{'discrete','discrete'},[nbin 2]);
        mu =  tapply(rt,{v,treatment},@nanmean,{'discrete','discrete'},[nbin 2]);
        err =   tapply(rt,{v,treatment},@sem,{'discrete','discrete'},[nbin 2]);
        
        for i = 1:2
           [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('E1 + E2'); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
%%
    fig = figure; set(fig,'Name','rt13');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        nt = data.trialNumber(selection);

        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        
  
        nbin = 10;
        x = tapply(nt,{nt,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        mu =  tapply(rt,{nt,treatment},@nanmean,{'continuous','discrete'},[nbin 2]);
        err =   tapply(rt,{nt,treatment},@sem,{'continuous','discrete'},[nbin 2]);
        
        for i = 1:2
%            [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});
          [~,~,h(i)] = errorscat(x(:,i) ,err(:,i), err(:,i)*0,col{i});

        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('session progression (total trial)'); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
             
  %%  
    fig = figure; set(fig,'Name','rt14');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        sv = data.participationValue(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
        r1 = data.ordinalRewardLeft(selection);
        r2 = data.ordinalRewardRight(selection);
        e1 = data.ordinalEffortLeft(selection);
        e2 = data.ordinalEffortRight(selection);
        dr = round(r2-r1);
        de = round(e2-e1);
        congruent = ((dr.*(-de))>0)*2;
        congruent((de.*dr)==0) = 1;
        time = data.offerTime(selection);

  
        nbin = [6 2 2];
        x = tapply(sv,{sv,congruent,treatment},@nanmean,{'continuous','discrete','discrete'},nbin);
        mu =  tapply(rt,{sv,congruent,treatment},@nanmean,{'continuous','discrete','discrete'},nbin);
        err =   tapply(rt,{sv,congruent,treatment},@sem,{'continuous','discrete','discrete'},nbin);
        x = tapply(sv,{sv,time,treatment},@nanmean,{'continuous','continuous','discrete'},nbin);
        mu =  tapply(rt,{sv,time,treatment},@nanmean,{'continuous','continuous','discrete'},nbin);
        err =   tapply(rt,{sv,time,treatment},@sem,{'continuous','continuous','discrete'},nbin);     
        
        for i = 1:2
           [~,~,a] = errorscat(x(:,1,i) ,mu(:,1,i), err(:,1,i),col{i}); a.LineStyle = '-.';
           [~,~,h(i)] = errorscat(x(:,2,i) ,mu(:,2,i), err(:,2,i),col{i});
%            [~,~,b] = errorscat(x(:,3,i) ,mu(:,3,i), err(:,3,i),col{i}); b.LineStyle = '--';
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('offer value '); 
        ylabel('response time (sec.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end
    
    
%%
fig = figure; set(fig,'Name','rt7');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
%         selection = selectSubject & selectionChosen & ~selectionRepetition;
        session = data.session(selection);
%         time = data.intertrialTime(selection);
        time = data.offerTime(selection);

        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
  
        nbin = [2 2];
        x =  tapply(time,{time,treatment},@nanmean,{'continuous','discrete'},nbin);
        mu =  tapply(rt,{time,treatment},@nanmean,{'continuous','discrete'},nbin);
        err =   tapply(rt,{time,treatment},@sem,{'continuous','discrete'},nbin);
        
       for i = 1:2
            [~,~,h(i)] = errorscat(x(:,i) ,mu(:,i), err(:,i),col{i});

       end
%         plot(x(ind,i),x(ind,i),'k--');
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        xlabel('offer time (s.)'); 
        ylabel('response time (s.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end

%%
fig = figure; set(fig,'Name','rt8');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
%         time = data.intertrialTime(selection);
        time = data.offerTime(selection);
        r1 = data.ordinalRewardLeft(selection);
        r2 = data.ordinalRewardRight(selection);
        e1 = data.ordinalEffortLeft(selection);
        e2 = data.ordinalEffortRight(selection);
        dr = round(r2-r1);
        de = round(e2-e1);
        congruent = ((dr.*(-de))>0)*2;
        congruent((de.*dr)==0) = 1;
        
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);        
  
        nbin = [3 2 2];
        x =  tapply(congruent,{congruent,time,treatment},@nanmean,{'discrete','continuous','discrete'},nbin);
        mu =  tapply(rt,{congruent,time,treatment},@nanmean,{'discrete','continuous','discrete'},nbin);
        err =   tapply(rt,{congruent,time,treatment},@sem,{'discrete','continuous','discrete'},nbin);
        
       for i = 1:2
            [~,~,h(i)] = errorscat(x(:,1,i) ,mu(:,1,i), err(:,1,i),col{i}); h(i).LineStyle = '--';
            [~,~,h(i)] = errorscat(x(:,2,i) ,mu(:,2,i), err(:,2,i),col{i}); h(i).LineStyle = '-';

       end
%         plot(x(ind,i),x(ind,i),'k--');
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = x(:,1,i);
        ax.XTickLabel = {'conflicting','unidimensional','congruent'};
        ax.XLim = [-1 3];
        xlabel('dimension congruence'); 
        ylabel('response time (s.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end

%%
fig = figure; set(fig,'Name','rt8');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject & selectionChosen ;
        session = data.session(selection);
        novelty = (~selectionRepetition(selection));

        treatment = data.treatment(selection);
        rt = data.reactionTime(selection);
  
        mu =  tapply(rt,{novelty,treatment},@nanmean);
        err =   tapply(rt,{novelty,treatment},@sem);
        
       for i = 1:2
            ind = [1:2];
%             h(i) = scatter( x(ind,i) ,mu(ind,i),10,...
%                      'MarkerFaceColor',col{i});
             h(i) = plot( ind ,mu(ind,i),...
                     'Color',col{i},'LineWidth',2);
             errbar( ind ,mu(ind,i), err(ind,i) ,...
                     'Color',col{i},'LineWidth',2);
       end
%         plot(x(ind,i),x(ind,i),'k--');
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca;
        ax.XTick = [1 2];
        ax.XLim = [0 3];
        ax.XTickLabel = {'repetition','novel'};
        ylabel('response time (s.)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end

    
%% regressions analysis

%% 1) force
% select
 selection = selectionChosen & ~selectionRepetition;
%  selection = selectionChosen & ~selectionRepetition & data.treatment=='placebo';

 y = data.perf(selection);        
%  y = data.cumulativePerf(selection);        

% reward effect
     x3 = data.chosenCardinalReward(selection);


% effort effect
     x1 = data.sideChoice(selection);
     x2 = data.chosenCardinalEffort(selection);


% dimension integration
%     x8 = x6.*x7;
%     x8 = x2.*(-x1);

     
% dynamical effect
     % short-term
%      x4 = y(1:end-1); x4 = [NaN; x4]; % previousForce
%      x8 = data.ordinalRewardOutcome(selection); x8 = x8(1:end-1); x8 = [NaN; x8]; % previousOutcome
%      x8 = data.offerTime(selection);
     x8 = data.reactionTime(selection) - data.offerTime(selection);

     
     % long-term
      x4 = data.trialNumber(selection);
%       x8 = data.cumulativeReward(selection);
%       x4 = data.cumulativeEffort(selection);
%       x8 = data.rewardRate(selection);
%       x4 = data.effortRate(selection);



 % treatment effect
    x5 = data.treatment(selection);

  % fit
  f = figure; hold on;
%   f.Units = 'normalized';
%   f.Position = [0.05 0.05 0.6 0.9];
  for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);

         predictor = table(nanzscore(x1),nanzscore(x2),nanzscore(x3),nanzscore(x4),nanzscore(x8),(x5) ,(y) );

         predictor.Properties.VariableNames = {'side';'E';'R';'ntrial';'rt';'atx';'force'};
         predictor = predictor(selectSubject(selection),:);

         formula = 'force ~ E + side + ntrial + rt:E';
         formula = 'force ~ E + side + ntrial + rt:E + atx:(1 +  side)';


         glm = fitglm( predictor,...
            formula,...
            'Distribution','normal','link','identity','CategoricalVars',[6]);

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
  
  %% 2 ) choice
% select
 selection = selectionChosen;
%   selection =  selectionChosen & ~selectionRepetition & data.treatment=='placebo';

 y = data.sideChoice(selection)-1;     
 
% reward effect
    x1 = data.deltaOrdinalReward(selection);
% effort effect
    x2 = data.cardinalEffortRight(selection);
    x3 = data.cardinalEffortLeft(selection);
    
% dimension integration
%     x8 = x6.*x7;
%     x8 = x2.*(-x1);

% dynamical effect
     % short-term
     x4 = y(1:end-1); x4 = [NaN; x4]; % previousChoice
%      x8 = data.ordinalRewardOutcome(selection); x8 = x8(1:end-1); x8 = [NaN; x8]; % previousOutcome
%      x8 = data.offerTime(selection);

     
     % long-term
      x5 = data.trialNumber(selection);
%       x8 = data.cumulativeReward(selection);
%       x8 = data.cumulativeEffort(selection);
%       x8 = data.rewardRate(selection);
%       x8 = data.effortRate(selection);


 % treatment effect
    x10 = data.treatment(selection);

  % fit
  f = figure; hold on;
%   f.Units = 'normalized';
%   f.Position = [0.05 0.05 0.6 0.9];
  for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);

         predictor = table(rangescore(x1),rangescore(x2),rangescore(x3),nanzscore(x4),(x10) , (y) );
         predictor.Properties.VariableNames = {'dR';'E2';'E1';'choice_t_1';'atx';'choice'};
%          predictor = predictor(selectSubject(selection),:);
%          formula = 'choice ~ dR + dE + sumR:dR + sumE:dE + sumE:dR + sumR:dE ';
%          formula = 'choice ~ dR + dE + choice_t_1';
%          formula = 'choice ~ dR + dE + choice_t_1 + atx + atx:dE + atx:dR + atx:choice_t_1  ';
         formula = 'choice ~ -1 + dR + E2 + E1 + atx:( dR + E2 + E1) ';


         glm = fitglm( predictor,...
            formula,...
            'Distribution','binomial','link','logit','CategoricalVars',[5]);

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
  
  
%% 3 ) participation model
% select

 selection =  ~selectionRepetition |  selectionRepetition;
%  selection =  selectionAll ;
%  selection =  selectionAll & data.treatment=='1_placebo';
 c = col{1};
%  selection =  selectionAll & data.treatment=='2_atx';
%  c = col{2};

 y = selectionParticipate(selection);       
 
 % reward effect
     x1 = data.ordinalRewardLeft(selection);
     x2 = data.ordinalRewardRight(selection);
     x3 = x1 + x2;
%      x3 = max([x1,x2],[],2);
% effort effect
     x4 = data.ordinalEffortLeft(selection);
     x5 = data.ordinalEffortRight(selection);
     x6 = x4 + x5;
%      x6 = min([x4,x5],[],2);
% dimension integration
%      x10 = abs(x1-x2).*abs(x3-x4);
%      x10 = (x2-x1).*(x3-x4);
%      x10 = x3.*x6;
     
% dynamical effect
     % short-term
%          x10 = data.ordinalRewardOutcome(selection); x7 = x7(1:end-1); x7 = [NaN; x7]; % previousOutcome
         x7 = y(1:end-1); x7 = [NaN; x7]; % previousParticipation
%          x10 = data.offerTime(selection);
%          x10 = data.intertrialTime(selection);

     % long-term
         x8 = data.trialNumber(selection);
    %      x8 = data.cumulativeReward(selection);
%          x8 = data.cumulativeEffort(selection);
             x10 = data.rewardRate(selection);

 % treatment effect
    x9 = data.treatment(selection);


%   % fit
%   f1 = figure; hold on;
%   f2 = figure; hold on;
%   f2.Name = 'fit';
% %   f.Units = 'normalized';
% %   f.Position = [0.05 0.05 0.6 0.9];
%   for iSub = 1:numel(subjectList)
%         sub = subjectList{iSub};
%         selectSubject = ismember(data.subject(selection), sub);
% 
%          predictor = table(nanzscore(x1),nanzscore(x2),nanzscore(x3),nanzscore(x4),nanzscore(x5),nanzscore(x6),nanzscore(x7), nanzscore(x8) , nanzscore(x10) , (x9) , (y) );
% %          predictor = table(rangescore(x1),rangescore(x2),rangescore(x3),rangescore(x4),rangescore(x5),rangescore(x6),rangescore(x7), rangescore(x8) , rangescore(x10) , (x9) , (y) );
%          predictor.Properties.VariableNames = {'R1';'R2';'sumR';'E1';'E2';'sumE';'participation_t_1';'ntrial';'rewardrate';'atx';'participation'};
%          formula = 'participation ~ 1 + R1 + R2 + E1 + E2 + ntrial + atx:( 1 + R1 + R2 + E1 + E2 + ntrial )';
% 
%          
%          predictor = predictor(selectSubject(selection),:);
%          glm = fitglm( predictor,...
%             formula,...
%             'Distribution','binomial','link','logit','CategoricalVars',[10]);
% 
%         figure(f1);
%         subplot(numel(subjectList),1,iSub);hold on
%         ax = gca;
%         option.pValue=0;
%         f = displayGLMcoeff( glm , f1 , ax, option );
%         ax.FontSize = 14;
% %         ax.Color = c;
%         if iSub~= numel(subjectList)
%             ax.XTick = [];
%             ax.XLabel = [];
%         end
%         ylabel(ax,'% of variance');
%         t = title(ax,['monkey ' sub(1) ]); 
%         
%         figure(f2);
%         subplot(numel(subjectList),1,iSub);hold on
%         nbin = 10;
%         x = glm.Fitted.LinearPredictor ;
%         y1 = predictor.participation;
%         f  = glm.Link.Inverse;
%         xx = tools.tapply(x,{x},@nanmean,{'continuous'},nbin);
%         yy = tools.tapply(y1,{x},@nanmean,{'continuous'},nbin);
%         zz = tools.tapply(y1,{x},@sem,{'continuous'},nbin);
%         plot(sort(x),f(sort(x)),'-','Color',c);
%         scatter(xx,yy,'MarkerEdgeColor',c,'MarkerFaceColor',c,'LineWidth',2);
% %         errbar(xx,yy,zz,'Color',c,'LineWidth',2);
%         ax = gca;
%         xlabel('decision value'); 
%         ylabel('participation (%)');
%         t = title(ax,['monkey ' sub(1) ]); 
% 
%   end
  
  % model selection
  for iSub = 1:numel(subjectList)
      
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject(selection), sub);
        
         predictor = table( rangescore(x1),rangescore(x2),rangescore(x3),rangescore(x4),rangescore(x5),rangescore(x6),rangescore(x7), rangescore(x8)  , (x9) , (y) );
%          predictor = table(nanzscore(x1),nanzscore(x2),nanzscore(x3),nanzscore(x4),nanzscore(x5),nanzscore(x6),nanzscore(x7), nanzscore(x8)  , (x9) , (y) );
         predictor.Properties.VariableNames = {'R1';'R2';'sumR';'E1';'E2';'sumE';'participation_t_1';'ntrial';'atx';'participation'};
         predictor = predictor(selectSubject(selection),:);
         
         predictor.atx = zscore(double(predictor.atx));
         x = predictor{:,[1 2 4 5 8 9]};
         for ix =  [1 2 4 5 8]
            x = [x , predictor{:,ix}.*predictor{:,9} ];
         end
          x = [x , predictor{:,10} ];
         x = array2table(x,'VariableNames',{'R1';'R2';'E1';'E2';'ntrial';...
                                            'atx';'atx_R1';'atx_R2';'atx_E1';'atx_E2';'atx_ntrial';...
                                            'participation'});
         c = [1 1 1 1 1 1  1 1 1 1 1 1 ;
              1 1 1 1 1 1  1 0 0 0 0 0 ;
              1 1 1 1 1 1  0 1 0 0 0 0 ;
              1 1 1 1 1 1  0 0 0 1 0 0 ;
              1 1 1 1 1 1  0 0 0 0 0 1 ;
              1 1 1 1 1 1  0 0 0 0 1 0 ;
              1 1 1 1 1 1  1 0 0 0 1 0 ];
          
        [b,stat] = VBA_logit_ATX1(x,[],[],c);
         
     end
  
  y = beta(1) + (beta(2:6)')*u(1:5); % placebo effects
y = y +  safepos(beta(7))*u(6) ;  % main effect
y = y +  safepos(beta(8))*(beta(2)*u(1)+beta(3)*u(2));
y = y +  safepos(beta(10))*(beta(4)*u(3)+beta(5)*u(4));
y = y +  safepos(beta(12))*(beta(6)*u(5));
%% 4) rt
% select
 selection = selectionChosen & ~selectionRepetition;
%  selection = selectionChosen & ~selectionRepetition & data.treatment=='placebo';

 y = data.reactionTime(selection);        
%  y = data.cumulativePerf(selection);        

% reward effect
     x3 = data.chosenCardinalReward(selection);
%      x3 = data.ordinalRewardLeft(selection) + data.ordinalRewardRight(selection);

% effort effect
     x1 = data.sideChoice(selection);
     x2 = data.chosenCardinalEffort(selection);
%      x2 = data.ordinalEffortLeft(selection) + data.ordinalEffortRight(selection);



% dimension integration
%     x6 = x2.*x3;
%     x6 = abs(data.deltaOrdinalReward(selection)).*abs(-data.deltaOrdinalEffort(selection));
%     x6 = abs(data.deltaOrdinalEffort(selection));
    x6 = (data.deltaOrdinalReward(selection)).*(-data.deltaOrdinalEffort(selection))>0;

     
% dynamical effect
     % short-term
%      x4 = y(1:end-1); x4 = [NaN; x4]; % previousForce
%      x8 = data.ordinalRewardOutcome(selection); x8 = x8(1:end-1); x8 = [NaN; x8]; % previousOutcome
     x8 = data.offerTime(selection);
%      x8 = data.reactionTime(selection) - data.offerTime(selection);

     
     % long-term
%       x4 = data.trialNumber(selection);
%       x4 = data.cumulativeReward(selection);
%       x4 = data.cumulativeEffort(selection);
%       x4 = data.rewardRate(selection);
      x4 = data.effortRate(selection);



 % treatment effect
    x5 = data.treatment(selection);

  % fit
  f = figure; hold on;
%   f.Units = 'normalized';
%   f.Position = [0.05 0.05 0.6 0.9];
  for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);

         predictor = table(rangescore(x6),rangescore(x2),rangescore(x3),rangescore(x4),rangescore(x8),(x5) , rangescore(y) );

         predictor.Properties.VariableNames = {'x6';'E';'R';'ntrial';'offertime';'atx';'rt'};
         predictor = predictor(selectSubject(selection),:);

         formula = 'rt ~ offertime + R ';
         formula = 'rt ~ offertime + R  + atx:(1 + offertime + R)';


         glm = fitglm( predictor,...
            formula,...
            'Distribution','gamma','link','identity','CategoricalVars',[6]);

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
  
%% 5) noise correlation

% force noise
% select
 selection = selectionChosen & ~selectionRepetition;
 y = data.perf(selection);         
 x2 = data.chosenCardinalEffort(selection);
 e = y-x2;
 
 y2 = data.sideChoice(selection)-1;       
 x1 = data.deltaOrdinalEffort(selection);
 x2 = data.deltaOrdinalReward(selection);
 x3 = data.reactionTime(selection)- data.offerTime(selection);
 trt = data.treatment(selection);
 x4 = data.chosenCardinalEffort(selection);

  % fit
  f = figure; hold on;
  for iSub = 1:numel(subjectList)
              
        subplot(1,numel(subjectList),iSub);hold on
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);

         predictor = table(zscore(x1),zscore(x2), (y2) );
         predictor.Properties.VariableNames = {'dE';'dR';'choice'};
         predictor = predictor(selectSubject(selection),:);
         formula = 'choice ~ -1 + dE + dR  ';
         glm = fitglm( predictor,...
            formula,...
            'Distribution','binomial','link','logit');
        
        e3 = glm.Residuals.Raw;
        e2 = e(selectSubject(selection));
        rt = x3(selectSubject(selection));
        x = x4(selectSubject(selection));
        
        for i = 1:2
             select = (trt(selectSubject(selection)) == trtList(i));
             h(i) = scatter((abs(e3(select))),(e2(select)),'filled','MarkerFaceColor',col{i},'SizeData',30);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca; 
        xlabel(ax,'choice error');
        ylabel(ax,'force error');

        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
        
  end
  
   % fit
  f = figure; hold on;
  for iSub = 1:numel(subjectList)
              
        subplot(1,numel(subjectList),iSub);hold on
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);

         predictor = table(zscore(x1),zscore(x2), (y2) );
         predictor.Properties.VariableNames = {'dE';'dR';'choice'};
         predictor = predictor(selectSubject(selection),:);
         formula = 'choice ~ -1 + dE + dR  ';
         glm = fitglm( predictor,...
            formula,...
            'Distribution','binomial','link','logit');
        
        e3 = glm.Residuals.Raw;
        e2 = e(selectSubject(selection));
        rt = x3(selectSubject(selection));
        x = x4(selectSubject(selection));
        
        for i = 1:2
             select = (trt(selectSubject(selection)) == trtList(i));
             h(i) = scatter((rt(select)),(e2(select)),'filled','MarkerFaceColor',col{i},'SizeData',30);
        end
        legend([h(1) h(2)],{'placebo','atx'});

        ax = gca; 
        xlabel(ax,'rt');
        ylabel(ax,'force error');
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
        
  end
  
%% Session analysis
fig = figure; set(fig,'Name','sessionEffect_1');
for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        
        selection = selectSubject ;
        
        session = data.session(selection);
        treatment = data.treatment(selection);
        
        participation = selectionChosen(selection);
        force = data.perf(selection);
        choice = data.sideChoice(selection)-1;
        rt = data.responseTime(selection);
        
        nt = data.trialNumber(selection);
        r1 = data.ordinalRewardLeft(selection);
        r2 = data.ordinalRewardRight(selection);
        e1 = data.ordinalEffortLeft(selection);
        e2 = data.ordinalEffortRight(selection);
        e = data.chosenCardinalEffort(selection);
        r = data.chosenCardinalReward(selection);

        % inti
        beta1 = nan(numel(sessionList),1);
        beta2 = nan(numel(sessionList),1);
        beta3 = nan(numel(sessionList),1);
        beta4 = nan(numel(sessionList),1); 
        beta5 = nan(numel(sessionList),1);
        beta6 = nan(numel(sessionList),1);
        beta7 = nan(numel(sessionList),1);
        beta8 = nan(numel(sessionList),1);   
        
        dev = nan(numel(sessionList),1);  
        accuracy = nan(numel(sessionList),1);  
        p0 = nan(numel(sessionList),1);  
        f0 = nan(numel(sessionList),1);  
        rt0 = nan(numel(sessionList),1);  
        
        % session loop
        for is = 1:numel(sessionList)
            
         % fit participation
             sel = (session==sessionList(is));
             predictor = table(zscore(r1(sel)),zscore(r2(sel)),zscore(e1(sel)),zscore(e2(sel)),zscore(nt(sel)),participation(sel) );
             predictor.Properties.VariableNames = {'R1';'R2';'E1';'E2';'ntrial';'participation'};
             formula = 'participation ~ 1 + R1 + R2 + E1 + E2 + ntrial';
             glm = fitglm( predictor,formula,'Distribution','binomial','link','logit');
             beta1(is) = glm.Coefficients.Estimate(1);
             beta2(is) = sum(glm.Coefficients.Estimate([2 3]));
             beta3(is) = sum(glm.Coefficients.Estimate([4 5]));
             beta4(is) = glm.Coefficients.Estimate(6);
             p0(is) = mean(participation(sel));

          % fit choice
             sel = (session==sessionList(is) & participation==1 );
             predictor = table(zscore(r2(sel)),zscore(r1(sel)),zscore(e2(sel)),zscore(e1(sel)),zscore(nt(sel)),choice(sel) );
             predictor.Properties.VariableNames = {'R2';'R1';'E2';'E1';'ntrial';'choice'};
             formula = 'choice ~ -1 + R2 + R1 + E1 + E2';
             glm = fitglm( predictor,formula,'Distribution','binomial','link','logit');
             beta5(is) = glm.Coefficients.Estimate(1);
             beta6(is) = sum(glm.Coefficients.Estimate([2 3]));
             accuracy(is) = mean((abs(glm.Residuals.Raw)<0.5));

             
          % fit force
          
             sel = (session==sessionList(is) & participation==1);
             predictor = table(zscore(e(sel)),zscore(choice(sel)),zscore(r(sel)),zscore(nt(sel)),(force(sel)) );
             predictor.Properties.VariableNames = {'E';'side';'R';'ntrial';'force'};
             formula = 'force ~ 1 + E*R + side + ntrial';
             glm = fitglm( predictor,formula,'Distribution','gamma','link','identity');
             beta7(is) = glm.Coefficients.Estimate(1);
             beta8(is) = sum(glm.Coefficients.Estimate([2]));
             dev(is) = mean(abs(glm.Residuals.Raw));
             f0(is) = mean(force(sel));
             
             rt0(is) = nanmean(rt(sel));


        end
        
        % concatenate
        metric = [p0 accuracy f0 dev rt0];
        metricName = {'participation','accuracy','force','variability','rt'};
        [cormat,pmat] = corr(metric);
        
        % plot
            figure(fig);
            subplot(1,numel(subjectList),iSub);hold on
            
            
            p = cormat.*(pmat<=0.05);
            p2 = cormat;col
%             cm = colormap('cool');
            cm = colormap_italian;
%             cm = flipud(cm);
            minCM = min(min(p))-0 ; maxCM = max(max(p))+0;
            minCM = min(min(p))*(1) ; maxCM = max(max(p))*(1);
            minCM = -1 ; maxCM = 1;

            heatmap(flipud(p),metricName,fliplr(metricName),flipud(p2),'Colormap',cm,'NaNColor',[1 1 1],...
            'MaxColorValue',maxCM,'MinColorValue',minCM,'TextColor',[1 1 1]*0);
            colorbar;

            ax = gca;
            ax.XTick = [1:size(p,1)];
            ax.XTickLabel = metricName;
            ax.YTick = [1:size(p,1)];
            ax.YTickLabel = metricName;
            yy =ylim;
            t = title(['monkey ' sub(1) ]);
            
            fig2 = figure; set(fig2,'Name',['sessionEffect_' sub(1)]);
            for i = 1:numel(metricName)
                for j = 1:numel(metricName)
                    
                    subplot(numel(metricName),numel(metricName),(i-1)*numel(metricName)+j);
                    hold on;
                    for t = 1:2
                        ind = trtSessions{t};
                        x = metric(ind,i);
                        y = metric(ind,j);
                        scatter(x,y,'MarkerFaceColor',col{t},'MarkerEdgeColor',col{t});
                        xlabel(metricName{i});
                        ylabel(metricName{j});
                    end
                    l = lsline;
                end
            end
end
        



