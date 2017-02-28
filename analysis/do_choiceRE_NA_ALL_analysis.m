%% do_choiceRE_NA_ALL_analysis

clc;
clear all;
close all;

%%
analysisname = ['monkeydrug_NA_dataset.mat'];
load(analysisname);

%% local variables
    % data acquisition parameters
        forceSamplingFreq = 25; % Hz.
        rewardUnit2Volume = 0.9; % conversion rate between valve opening time (sec) & water volume (ml)
        
        col = {[0.54 0.27 0.07],...
               [1 1 1]*0.5,...
               [1 0.5 0]};
        
        subjectList = {'aliosha','bob','esmeralda'};
        nsub = numel(subjectList);
        trtList = {'clonidine','placebo','atomoxetine'};
        ntrt = numel(trtList);
        sessionList = unique(data.session);
        maxSess = max(design.sessionNumber);
        sideList = {'left','right'};
        
        
%% reformat
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
        dR =  round(data.deltaOrdinalReward);
        dE =  round(data.deltaOrdinalEffort);
        selectionCongruent = (sign(dR)~=sign(dE));
        selectionCongruent(sign(dR)==0 & sign(dE)==0) = 0;
        
%% data preprocessing/completion
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
        data.successRate = ones(height(data),1);
        data.sessionNumber = zeros(height(data),1);
        data.successRate_left = zeros(height(data),1);
        data.successRate_right = zeros(height(data),1);

%         data.force_byEffort = zeros(height(data),1);


        for iSub = 1:numel(subjectList)
            sub = subjectList{iSub};
            selectSubject = ismember(data.subject, sub);
            sessList = unique(data.session(selectSubject));
            for iSess = 1:numel(sessList)
                
                % display
                fprintf('subject %d ; session %d \n',iSub,iSess);
                
                selectSession = (data.session==sessList(iSess)) ;
                select = selectSubject & selectSession;
                
                % re-assign treatment per session
                treatment = design.treatment(ismember(design.session,sessList(iSess)) & ismember(design.subject,sub) );
                data.treatment(select) = repmat(treatment,numel(data.treatment(select)),1);
                sessionNumber = design.sessionNumber(ismember(design.session,sessList(iSess)) & ismember(design.subject,sub) );
                data.sessionNumber(select) = repmat(sessionNumber,numel(data.treatment(select)),1);

                cumRew =   nancumsum(data.volumeRewardOutcome(select)); 
%                 cumRew(2:end) = cumRew(1:end-1);cumRew(1)= 0;
                data.cumulativeReward(select) = cumRew;
                
                cumPerf =  nancumsum(data.perf(select)); 
%                 cumPerf(2:end) = cumPerf(1:end-1);cumPerf(1)= 0;
                data.cumulativeEffort(select) = cumPerf;
                
                success = data.ordinalRewardOutcome(select)>0;
                cumSucess = nancumsum(success); 
                cumPart = nancumsum(selectionParticipate(select)); 
                chosenEffort = data.chosenOrdinalEffort(select);
                nt = data.trialNumber(select);
                data.successRate(select) = cumSucess./cumPart;
                p1 = nan(max(nt),1); p2 = nan(max(nt),1);
                for it=1:max(nt)
                    e1 = data.ordinalEffortLeft(it);
                    e2 = data.ordinalEffortRight(it);
                    p1(it) = nanmean(success(chosenEffort==e1 & nt<=it));
                    p2(it) = nanmean(success(chosenEffort==e2 & nt<=it));
                end
                data.successRate_left(select) = p1;
                data.successRate_right(select) = p2;

                
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

% save
save(analysisname,'data','-append');
     

% add session-number info
ethogram.sessionNumber = zeros(height(ethogram),1);
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(ethogram.subject, sub);
        sessList = unique(ethogram.session(selectSubject));
        for iSess = 1:numel(sessList)
            selectSession = (ethogram.session==sessList(iSess)) ;
            select = selectSubject & selectSession;
            sessionNumber = design.sessionNumber(ismember(design.session,sessList(iSess)) & ismember(design.subject,sub) );
            ethogram.sessionNumber(select) = repmat(sessionNumber,numel(ethogram.treatment(select)),1);
        end
    end
save(analysisname,'ethogram','-append');


%% aggregate descriptive statistics: session-level

varnames = {'session','session_number','weekday','subject','treatment'...
            'session_length','session_time','session_reward',...
            'participation_rate','participation_ntrial','participation_inertia','omission_inertia','premature_rate'...
            'choice_R','choice_E','accuracy_R','accuracy_E','accuracy_RE','choice_HRHE','choice_consistency','choice_repetition',...
            'peak_force','varcoef_force','peak_velocity','time2peak','success_force',...
            'responsetime','varcoef_responsetime'};
nvar = numel(varnames);
sessiondata = design(:,1:5);
sessiondata = [sessiondata , array2table(nan(height(design),nvar-5))];
sessiondata.Properties.VariableNames = varnames;

for isub = 1:nsub
    
    % select
    indsub = ismember(design.subject,subjectList{isub}) ;
    indsession = design.sessionNumber(indsub);
    selectSubject = ismember(data.subject, subjectList{isub});   
    
    selection = selectSubject ;
    
    % session info
    session = data.sessionNumber(selection);
    subsession = indsession;
    
    % participation
    participation = selectionParticipate(selection);
    nt = data.trialNumber(selection);
    sessiondata.session_length(indsub) = tapply(nt,{session},@max)';
    sessiondata.participation_rate(indsub) = tapply(participation,{session},@nanmean)';
    sessiondata.participation_ntrial(indsub) = tapply(participation,{session},@nansum)';
    
    time = data.cumulativeTime(selection);
    sessiondata.session_time(indsub) = tapply(time,{session},@nanmax)';
    
    reward = data.cumulativeReward(selection);
    sessiondata.session_reward(indsub) = tapply(reward,{session},@nanmax)';
    
    % premature
    premature = selectionPremature(selection);
    sessiondata.premature_rate(indsub) = tapply(premature,{session},@nanmean)';
    
    % consecutive states
    dwork = nan(numel(subsession),1);
    dpause = nan(numel(subsession),1);
    for is = [1:numel(subsession)]
        selection = (session==subsession(is));
        partsession = participation(selection); 
        nextpartsession = [ partsession(2:end); NaN ];  
        dwork(is) = mean(nextpartsession(partsession==1)==1); 
        dpause(is) = mean(nextpartsession(partsession==0)==0);
    end
    sessiondata.participation_inertia(indsub) = dwork;
    sessiondata.omission_inertia(indsub) = dpause;
    
    
    % select
    selection = selectSubject & selectionChosen;  
    % session info
    session = data.sessionNumber(selection);
    
    % choice
    choice = (data.sideChoice(selection)-1)*2-1;
    nnt = data.normalizedTrialNumber(selection);
    nextchoice = [ choice(2:end); NaN ];  
    nextchoice(nnt==1) = NaN; 
    repetition = (choice==nextchoice);
    dR = round(data.deltaOrdinalReward(selection));
    choiceR = (sign(choice)==sign(dR));
    dE = round(data.deltaOrdinalEffort(selection));
    choiceE = (sign(choice)==-sign(dE));
    accuracy = (choiceR | choiceE) ; 
    congruent = selectionCongruent(selection);
    sessiondata.choice_R(indsub) = tapply(choiceR,{session},@nanmean)';
    sessiondata.choice_E(indsub) = tapply((1-choiceE),{session},@nanmean)';
    sessiondata.accuracy_R(indsub) = tapply(choiceR(dR~=0 & dE==0),{session(dR~=0 & dE==0)},@nanmean)';
    sessiondata.accuracy_E(indsub) = tapply(choiceE(dE~=0 & dR==0),{session(dE~=0 & dR==0)},@nanmean)';
    sessiondata.accuracy_RE(indsub) = tapply(choiceR(congruent & dE~=0 & dR~=0),{session(congruent & dE~=0 & dR~=0)},@nanmean)';
    sessiondata.choice_HRHE(indsub) = tapply(choiceR(~congruent & dE~=0 & dR~=0),{session(~congruent & dE~=0 & dR~=0)},@nanmean)';
    sessiondata.choice_repetition(indsub) = tapply(repetition,{session},@nanmean)';

    % force
    force = data.perf(selection);
    success_force = (data.volumeRewardOutcome(selection)>0);
    t2p = data.time2peak(selection);
    velocity = data.yankPeak(selection);
    sessiondata.peak_force(indsub) = tapply(force,{session},@nanmean)';
    sessiondata.peak_velocity(indsub) = tapply(velocity,{session},@nanmean)';
    sessiondata.time2peak(indsub) = tapply(t2p,{session},@nanmean)';
    sessiondata.success_force(indsub) = tapply(success_force,{session},@nanmean)';

    varcoef  = @(x) nanvar(x)/nanmean(x);
    rchoice = data.chosenCardinalReward(selection);
    echoice = data.chosenCardinalEffort(selection);
    ysub = tapply(force,{session,rchoice,echoice},varcoef);
    sessiondata.varcoef_force(indsub) = nanmean(nanmean(ysub,3),2);
    
    % rt
    rt = data.responseTime(selection);
    sessiondata.responsetime(indsub) = tapply(rt,{session},@nanmean)';
    ysub = tapply(rt,{session,rchoice,echoice},varcoef);
    sessiondata.varcoef_responsetime(indsub) = nanmean(nanmean(ysub,3),2);

    
end


save(analysisname,'sessiondata','-append');



%% aggregate inferential statistics, subject-level

varnames = {'session','session_number','weekday','subject','treatment'...
            'session_length','session_time','session_reward',...
            'participation_rate','participation_ntrial','participation_inertia','omission_inertia','premature_rate'...
            'accuracy_R','accuracy_E','accuracy_RE','choice_HRHE','choice_consistency','choice_repetition',...
            'peak_force','varcoef_force','peak_velocity','time2peak',...
            'responsetime','varcoef_responsetime'};
nvar = numel(varnames);
aggregateStat = array2table(nan(nsub*2*3,nvar-5));
aggregateStat.Properties.VariableNames = varnames([6:end]);
aggregateStat = [ cell2table(repmat({''},nsub*2*3,3)) , aggregateStat];
stat_varnames = [ varnames(4) , {'effect'} , {'stat'} , varnames([6:end])];
aggregateStat.Properties.VariableNames = stat_varnames;

sessiondata_2 = sessiondata;

for isub = 1:nsub
    
    % select
    indsub = ismember(design.subject,subjectList{isub}) ;
    selection = indsub ;
    ind = [1:6] + (isub-1)*6;

    aggregateStat.subject(ind) = repmat({subjectList{isub}},6,1);
    aggregateStat.effect(ind) =  [repmat({'clonidine'},3,1);repmat({'atomoxetine'},3,1)];
    aggregateStat.stat(ind) =  repmat([{'beta'};{'p'};{'h'}],2,1);


    % glm
    sessiondata.clonidine = (sessiondata.treatment=='clonidine');
    sessiondata.atomoxetine = (sessiondata.treatment=='atomoxetine');

    for ivar = [1:nvar-5]
        try
            var = varnames{ivar+5};
            tabstat = sessiondata(indsub,:);
            
            model = fitglm(tabstat,...
                    [ var ' ~ 1 +  session_number + weekday' ]);
%             disp(model);
            tabstat.(var) = model.Residuals.Raw + model.Coefficients.Estimate(1);
            sessiondata_2.(var)(selection) = tabstat.(var);
            model = fitglm(sessiondata(indsub,:),...
                    [ var ' ~ 1 + clonidine + atomoxetine ' ]);
%             disp(model);
            beta = model.Coefficients.Estimate(2:3);
            p = model.Coefficients.pValue(2:3);
            h = (p<=0.05);
            aggregateStat.(var)(ind) = [ beta(1) ; p(1) ; h(1) ; 
                                         beta(2) ; p(2) ; h(2) ];
        end
    end

end


%% 0.1) FLexible plot to compare dependant variables btw treatments
     
    % variables
    vartext = 'participation_rate';
    trt2comp = [1 2];
    ytext = 'participation_rate (%)';
    ylimits = [];
    
    nbin = [ntrt,nsub];
    Y = nan(nbin);
    Z = nan(nbin);

    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(sessiondata.subject, sub);   
        selection = selectSubject ;

        % variables
        tab = sessiondata;
        session = tab.session_number(selection);
        treatment = tab.treatment(selection);
        participation_rate = tab.participation_rate(selection);
        accuracy_R = tab.accuracy_R(selection);
        accuracy_E = 1 - tab.accuracy_E(selection);
        choice_HRHE = tab.choice_HRHE(selection);
        accuracy_RE = tab.accuracy_RE(selection);
        
        varcoef_force = tab.varcoef_force(selection);
        varcoef_responsetime = tab.varcoef_responsetime(selection);
        peak_force = tab.peak_force(selection);
        responsetime = tab.responsetime(selection);

        eval(['vary = ' vartext ]);

        
        % subject-display
        if isub==1; fig = figure; set(fig,'Name','vary_comp_sub'); end
        subplot(1,nsub,isub);hold on ; clear h
        for it=trt2comp
            % inputs
            x = [1];
            xx = x + (it-2)*0.25;
            indtrt = ismember(treatment,trtList(it));
            y = nanmean(vary(indtrt));  Y(it,isub) = y;
            z = sem(vary(indtrt));      Z(it,isub) = z;
            % plot
            [ h(it)  ] = barplot( xx ,y,z, col{it} );
            h(it).BarWidth = 0.25;
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend(trtList{trt2comp});
        end
        ax = gca;
        if ~isempty(ylimits); ax.YLim = ylimits; end
        ax.TickLength = [0 0];
        ax.XTick = []; 
        if isub==1 ; ylabel(ytext); end
        t = title(['monkey ' sub(1) ]); 
        
    end

     % group-display
        fig = figure; set(fig,'Name','vary_comp_group');
        hold on ; clear h
        for it=trt2comp
            % inputs
            x = [1];
            xx = x + (it-2)*0.25;
            y = nanmean(Y(it,:));  
            z = nanmean(Z(it,:));       
            % plot
            [ h(it)  ] = barplot( xx ,y,z, col{it} );
            h(it).BarWidth = 0.25;
        end
        % legending
        legend(trtList{trt2comp});
        ax = gca;
        if ~isempty(ylimits); ax.YLim = ylimits; end
        ax.TickLength = [0 0];
        ax.XTick = []; 
        ylabel(ytext); 
        
        % stat
        predictor = [];
        predictor = [predictor ,(tab.treatment=='clonidine')];
        predictor = [predictor ,(tab.treatment=='placebo')];
        predictor = [predictor , (tab.treatment=='atomoxetine')];
        y = tab.(vartext);
       
        contrast = [0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
%         contrast = contrast'; [pv,stat,df,all] = GLM_contrast(predictor,y,contrast,'t',1);
        model = fitglm(predictor,y,'linear','Intercept',false);
        
%         predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),tab.subject,y,...
%              'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
%         model2 = fitglme(predictor2,[ vartext ' ~ clonidine + atomoxetine + (clonidine|subject) + (atomoxetine|subject)']);
%         disp(model2.Coefficients);
        
        [p,F] = coefTest(model,contrast);
        xx = x + (trt2comp(trt2comp~=2)-2)*0.25;
        [s] = sigstar( num2cell(xx),p);
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
%% 0.1.2) FLexible plot to correlate cross-session measures
     
    % variables
    varxtext = 'success_force';
    varytext = 'choice_E';
    trt2comp = [1 2];
    xtext = 'correct force (%)';
    ytext = 'choice = HE (%)';
    xlimits = [0 1];
    ylimits = [0 1];
    xbin = [1:nbin(1)];
    xticks = [];
    xlabels = [];
    BETA = []; T = [];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(sessiondata.subject, sub);   
        selection = selectSubject ;

        % variables
        tab = sessiondata;
        session = tab.session_number(selection);
        treatment = tab.treatment(selection);
        participation_rate = tab.participation_rate(selection);
        accuracy_R = tab.accuracy_R(selection);
        accuracy_E = 1 - tab.accuracy_E(selection);
        choice_HRHE = tab.choice_HRHE(selection);
        accuracy_RE = tab.accuracy_RE(selection);
        varcoef_force = tab.varcoef_force(selection);
        varcoef_responsetime = tab.varcoef_responsetime(selection);
        peak_force = tab.peak_force(selection);
        responsetime = tab.responsetime(selection);
        success_force = tab.success_force(selection);
        choice_E = tab.choice_E(selection);
        
        eval(['varx = ' varxtext ]);
        eval(['vary = ' varytext ]);
        
         T = [T ; treatment ];
         BETA = [BETA ; [varx,vary] ];
         
        % subject-display
        if isub==1; fig = figure; set(fig,'Name','var_corr_sub'); end
        subplot(1,nsub,isub);hold on ; clear h
        for it=trt2comp
            % inputs
            indtrt = ismember(treatment,trtList(it));
            x = (varx(indtrt));  
            y = (vary(indtrt));  
            % plot
            [~,h(it),s] = errorscat(x,y,y.*0,col{it});
            s.LineStyle = 'none';
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend(trtList{trt2comp});
        end
        ax = gca;
        if ~isempty(ylimits); ax.YLim = ylimits; end
        ax.TickLength = [0 0];
        if isub==1 ; ylabel(ytext); end
        xlabel(xtext); 
        t = title(['monkey ' sub(1) ]); 
        
    end

     % group-display
        fig = figure; set(fig,'Name','vary_comp_group');
        hold on ; clear h
        for it=trt2comp
            % inputs
            indtrt = ismember(T,trtList(it));
            x = (BETA(indtrt,1));  
            y = (BETA(indtrt,2));  
            % plot
            [~,h(it),s] = errorscat(x,y,y.*0,col{it});
            s.LineStyle = 'none';
            l = lsline; 
        end
        
        % legending
        legend(trtList{trt2comp});
        for it=trt2comp
            l(it).Color = col{it};
        end
        ax = gca;
        if ~isempty(ylimits); ax.YLim = ylimits; end
        ax.TickLength = [0 0];
        xlabel(xtext); 
        ylabel(ytext); 
        
        % stat
        predictor = [];
        predictor = [predictor ,(tab.treatment=='clonidine')];
        predictor = [predictor ,(tab.treatment=='placebo')];
        predictor = [predictor , (tab.treatment=='atomoxetine')];
        predictor = [predictor , tab.(varxtext)];
        y = tab.(varytext);
       
        contrast = [0 0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
%         contrast = contrast'; [pv,stat,df,all] = GLM_contrast(predictor,y,contrast,'t',1);
        model = fitglm(predictor,y,'linear','Intercept',false);
        
%         predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),predictor(:,4),tab.subject,y,...
%              'VariableNames',{'clonidine','placebo','atomoxetine','confound','subject',vartext});
%         model2 = fitglme(predictor2,[ vartext ' ~ clonidine + atomoxetine + (clonidine|subject) + (atomoxetine|subject) + (confound|subject)']);
%         disp(model2.Coefficients);
        
        [p,F] = coefTest(model,contrast);
        xx = x + (trt2comp(trt2comp~=2)-2)*0.25;
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
    
%% 0.2) FLexible plot to compare dependant variables btw treatments

    % variables
    vartext = 'peak_force';
    binfactor = [];
    factor = 'rchoice';
    trt2comp = [2 3];
    xtext = texlabel('R_chosen');
    ytext = 'force peak (%)';
    xlimits = [0.5 4.5];
    ylimits = [0 1];
    xxplot = 1;
    nqtle = 4;
    nbin = [ nqtle , ntrt , maxSess ];
    nbin2 = [ nqtle ,ntrt,nsub ];
    xbin = [1:nbin(1)];
    xticks = [];
    xlabels = [];
%     xlabels = {'incorrect','correct'};
    fstat = @nanmean;
    dataselect = 'selectSubject & selectionChosen';
    
    X = nan(nbin2);
    Y = nan(nbin2);
    Z = nan(nbin2);

    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        eval(['selection = ' dataselect ';']);
        
        % variables
        xsub = nan(nbin);
        ysub = nan(nbin);
        
        session = data.sessionNumber(selection);
        treatment = data.treatment(selection);
        
        participation_rate = selectionParticipate(selection);
%         participation_prediction = data.participation_prediction(selection);
%         participation_predictionbin = quantileranks(participation_prediction,nbin(1)); 
%         participation_predictionbin(participation_predictionbin==0)=NaN;
        
        choice = (data.sideChoice(selection)-1)*2-1;
        sdR = (round(data.deltaOrdinalReward(selection)));
        choiceR = (sign(choice)==sign(sdR));
        dR = abs(sdR);
        sdE = (round(data.deltaOrdinalEffort(selection)));
        choiceE = (sign(choice)==sign(sdE));
        dE = abs(sdE);
        r1 = round(data.ordinalRewardLeft(selection));
        r2 = round(data.ordinalRewardRight(selection));
        e1 = round(data.ordinalEffortLeft(selection));
        e2 = round(data.ordinalEffortRight(selection));
        sum_r = r1+r2;
        sum_e = e1+e2;
        congruency = double(sign(sdR)==-sign(sdE)) ;
        dimensionality = (sdR~=0) + (sdE~=0) ;
        congruency(dimensionality~=2)=NaN;
        
        offertime = data.intertrialTime(selection);
        offertimebin = quantileranks(offertime,nbin(1)); 

        
        accuracy = (choiceR | ~choiceE) ;
        responsetime = data.responseTime(selection);

        responsetimebin = quantileranks(responsetime,nbin(1)); 
        responsetimebin(responsetimebin==0)=NaN;
        echoice = data.chosenCardinalEffort(selection);
        rchoice = data.chosenCardinalReward(selection);
        peak_force = data.perf(selection);   
        peak_force_bin = quantileranks(peak_force,nbin(1)); 
        yank = data.yankPeak(selection);   
        time2peak = data.time2peak(selection);   
        nnt = data.normalizedTrialNumber(selection);
        nnt = quantileranks(nnt,nbin(1)); 
        nt = data.trialNumber(selection);
        tmax=1; p_t = nan(numel(participation_rate),tmax);
        for it=1:tmax
            p_t(1+it:end,it) = participation_rate(1:end-it);  p_t(nt<=it,it) = NaN;
        end
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);
        
        eval(['varx = ' factor ';' ]); 
        eval(['vary = ' vartext ';' ]); 
        if ~isempty(binfactor) ; 
        eval(['varx2 = ' binfactor ';' ]); 
        else
            varx2 = varx;
        end
        x = tapply(varx,{varx2,treatment,session},@nanmean);
        xsub(:,subtrt,subsess) =  x;
        y = tapply(vary,{varx2,treatment,session},fstat);
        ysub(:,subtrt,subsess) =  y;
            

        % second-level stat
         dim=3;
         x = nanmean(xsub,dim); X(:,subtrt,isub) = x;
         y = nanmean(ysub,dim); Y(:,subtrt,isub) = y;
         z = sem(ysub,dim);     Z(:,subtrt,isub) = z;
        
        % display
        if isub==1; fig = figure; set(fig,'Name','bivar_comp_sub'); end
        subplot(1,nsub,isub);hold on ; clear h
        for it=trt2comp
            xx = xbin ;
            [~,~,h(it)] = errorscat( x(xx,it) ,y(xx,it), z(xx,it),col{it});
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            h2leg = []; for it=trt2comp ; h2leg = [h2leg, h(it)];end;
            legend(h2leg,trtList{trt2comp});
        end
        ax = gca;
        if ~isempty(ylimits); ax.YLim = ylimits; end
        ax.TickLength = [0 0];
        ax.XTick = []; if ~isempty(xticks);  ax.XTick = xticks; end
        if ~isempty(xlabels);  ax.XTickLabel = xlabels; end
        xlabel(xtext);
        ylabel(ytext); 
        t = title(['monkey ' sub(1) ]); 
        
    end

    % group-display
    fig = figure; set(fig,'Name','bivar_comp_group');
    hold on ; clear h;  dim=3;
    if xxplot ; plot(x(xx),x(xx),'--k'); end
    for it=trt2comp
        % inputs
        x = nanmean(X(:,it,:),dim);  
        y = nanmean(Y(:,it,:),dim);  
        z = nanmean(Z(:,it,:),dim);       
        % plot
        xx = xbin ;
        [~,~,h(it)] = errorscat( x(xx) ,y(xx), z(xx),col{it});
    end
    
    % legending
    h2leg = []; for it=trt2comp ; h2leg = [h2leg, h(it)];end;
    legend(h2leg,trtList{trt2comp});
    ax = gca;
    if ~isempty(xlimits); ax.XLim = xlimits; end
    if ~isempty(ylimits); ax.YLim = ylimits; end
    ax.TickLength = [0 0];
    ax.XTick = []; if ~isempty(xticks);  ax.XTick = xticks; end
    if ~isempty(xlabels);  ax.XTickLabel = xlabels; end
    xlabel(xtext);
    ylabel(ytext); 
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);    
    
%% 0.3) FLexible plot to compare GLM parameters between treatments
% variables
     trt2comp = [2 3];
     vartext = 'force';
     formula = ' force ~ 1 + side + rchoice*echoice + ntrial  ';
     varnames = {'side';'rchoice';'echoice';'ntrial';'Rtot';'Etot';'force'};
     effectName = {texlabel('beta_0'),'side','rchoice','echoice','ntrial','rchoice:echoice'};
     nregressor = 6;
     fname= 'identity';  distname = 'normal';
     dataselect = 'selectSubject & selectionChosen';
     predictor_select = '[ side*2-1 , rchoice, echoice , nt , rr , ee ];';


     nfactor = numel(varnames)-1;
     nbin = [ nregressor , ntrt , maxSess ];
     nbin2 = [ nregressor ,ntrt,nsub ];
     ysub = nan(nbin);
     Y = nan(nbin2);
     Z = nan(nbin2);
     BETA = []; T = [];

     data.([vartext '_prediction']) = nan(height(data),1);
     
    for isub = 1:nsub
                
         % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        eval(['selection = ' dataselect ';']);

        % variables
        ysub = nan(nbin);
        session = data.sessionNumber(selection);
        participation = selectionParticipate(selection);
        participation_prediction = data.participation_prediction(selection);
        treatment = data.treatment(selection);
        
        r1 = round(data.ordinalRewardLeft(selection));
        r2 = round(data.ordinalRewardRight(selection));
        sum_r = r1+r2;
        max_r = max(r1,r2);
        e1 = round(data.ordinalEffortLeft(selection));
        e2 = round(data.ordinalEffortRight(selection));
        sum_e = e1+e2;
        min_e = min(e1,e2);
        nt = data.trialNumber(selection);
        nnt = data.normalizedTrialNumber(selection);
        ee =  data.cumulativeEffort(selection);
        [b,~,stat] = glmfit(nnt,ee,'normal');
        ee = stat.resid + b(1);
        rr =  data.cumulativeReward(selection);
        [b,~,stat] = glmfit(nnt,rr,'normal');
        rr = stat.resid + b(1);
        rchoice = data.chosenCardinalReward(selection);
        echoice = data.chosenCardinalEffort(selection);
        side = data.sideChoice(selection)-1;   
        
        choice = (data.sideChoice(selection)-1)*2-1;
        sdR = (round(data.deltaOrdinalReward(selection)));
        choiceR = (sign(choice)==sign(sdR));
        dR = abs(sdR);
        sdE = (round(data.deltaOrdinalEffort(selection)));
        choiceE = (sign(choice)==sign(sdE));
        dE = abs(sdE);
        r1 = round(data.ordinalRewardLeft(selection));
        r2 = round(data.ordinalRewardRight(selection));
        e1 = round(data.ordinalEffortLeft(selection));
        e2 = round(data.ordinalEffortRight(selection));
        sum_r = r1+r2;
        sum_e = e1+e2;
        congruency = double(sign(sdR)==-sign(sdE)) ;
        dimensionality = (sdR~=0) + (sdE~=0) ;
        congruency(dimensionality~=2)=1;
        offertime = data.offerTime(selection);
        offertimebin = quantileranks(offertime,nbin(1)); 
        
        responsetime = data.responseTime(selection)+eps;
        force = data.perf(selection);   
        cum_part = nan(size(participation));

        
        % previousParticipation
        tmax=1;
        p_t = nan(numel(participation),tmax);
        for it=1:tmax
            p_t(1+it:end,it) = participation(1:end-it);
            p_t(nt<=it,it) = NaN;
        end
        
        cumul = data.stateDuration(selection);
        cumul(participation==0) = 0;
        cumul2 = data.stateDuration(selection);
        cumul2(participation==1) = 0;
        tmax=1;
        work_t = nan(numel(cumul),tmax);
        pause_t = nan(numel(cumul),tmax);
        for it=1:tmax
            work_t(1+it:end,it) = cumul(1:end-it);
            work_t(nt<=it,it) = NaN;
            pause_t(1+it:end,it) = cumul2(1:end-it);
            pause_t(nt<=it,it) = NaN;
        end
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);

        for is = subsess'
            cum_part(session==(is)) = cumsum(participation(session==(is)));
            
            itrt = double(unique(treatment(session==(is))));
            eval(['predictor = ' predictor_select]);
            predictor = (predictor(session==(is),:));
            predictor = nanzscore(predictor);
            eval(['y = ' vartext '(session==is);']);
            stat = fitglm(predictor,y,...
                           formula,'VarNames',varnames,'Distribution',distname,'link',fname);
            coef = stat.Coefficients;
            disp(coef);
            ysub(:,itrt,is) = coef.Estimate;
            data.([vartext '_prediction'])(selection & data.sessionNumber==is) = stat.Fitted.LinearPredictor;
        end

           
        % second-level stat
         dim=3;
         T = [T ; nominal(repmat({'clonidine','placebo','atomoxetine'}',maxSess,1)) ];
         BETA = [BETA ; reshape(ysub,[nregressor,ntrt*maxSess])' ];
         y = nanmean(ysub,dim); Y(:,subtrt,isub) = y;
         z = sem(ysub,dim); Z(:,subtrt,isub) = z;
        
        % display
        if isub==1; fig = figure; set(fig,'Name','glm_participation'); end
        subplot(nsub,1,isub);hold on ; clear h ; h = gobjects(ntrt,1);
        for it=trt2comp
            x = [1:nregressor];
            xx = x + (it-2)*0.25;
            [ h(it)  ] = barplot( xx ,y(:,it),z(:,it), col{it} );
            h(it).BarWidth = 0.25;
            
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            h2leg = []; for it=trt2comp ; h2leg = [h2leg, h(it)];end;
%             legend(h2leg,trtList(trt2comp));
        end
        ax = gca;
%         ax.YLim = [0.5 1];
        ax.TickLength = [0 0];
        ax.XTick = []; 
        ax.XTick = x;
        ax.XTickLabel = effectName(x);
        ax.XTickLabelRotation = 45;
        ax.XLim = [0 x(end)+1];
        ylabel('regression coefficients '); 
        t = title(['monkey ' sub(1) ]); 
        
    end
% group-display
    fig = figure; set(fig,'Name','glm_comp_group');
    hold on ; clear h; h = gobjects(ntrt,1);  dim=3;
    for it=trt2comp
        % inputs
        y = nanmean(Y(:,it,:),dim);  
        z = nanmean(Z(:,it,:),dim);       
        % plot
        x = [1:nregressor];
        xx = x + (it-2)*0.25;
        [ h(it)  ] = barplot( xx ,y,z, col{it} );
        h(it).BarWidth = 0.25;
    end
    
    % legending
%     legend(h2leg,trtList(trt2comp));
    ax = gca;
    ax.TickLength = [0 0];
    ax.XTick = []; 
    ax.XTick = x;
    ax.XTickLabel = effectName(x);
    ax.XTickLabelRotation = 45;
    ax.XLim = [0 x(end)+1];
    ylabel('regression coefficients ');     
    
     % stat
    T = T(~isnan(BETA(:,1)),:);
    predictor = [];
    predictor = [predictor ,(T=='clonidine')];
    predictor = [predictor ,(T=='placebo')];
    predictor = [predictor , (T=='atomoxetine')];
    BETA = BETA(~isnan(BETA(:,1)),:);
    p = nan(nregressor,1);
    for i = 1:nregressor
        contrast = [0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
        y = BETA(:,i);
%         model = fitglm(predictor,y,'linear','Intercept',false);
%         [p(i),F] = coefTest(model,contrast);

        predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),tab.subject,y,...
             'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
        model2 = fitglme(predictor2,[ vartext ' ~ clonidine + atomoxetine + (clonidine|subject) + (atomoxetine|subject)']);
%         disp(model2.Coefficients);
        ind = max(trt2comp);
        p(i) =  model2.Coefficients.pValue(ind);
    end
    xx = x + (trt2comp(trt2comp~=2)-2)*0.25;
    [s] = sigstar( num2cell(xx),p);
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);  


%% 1) Trial classification

    % variables
        nbin = [ 5 , ntrt , maxSess ];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        selection = selectSubject ;

        % variables
        Y = nan(nbin);
        session = data.session(selection);
        trialtype = data.errorType(selection);
        treatment = data.treatment(selection);
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        [~,subsess] = ismember(unique(session),sessionList);
        subtype = unique(trialtype);
        subtype = subtype+1;
        ysub = tapply(trialtype,{trialtype,treatment,session},@numel);
        ysub(ysub==0) = NaN;
        Y(subtype,subtrt,subsess) =  ysub;

        % averaging
         dim=3;
         y = nanmean(Y,dim);
         z = sem(Y,dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','n_trials'); end
        subplot(1,nsub,isub);hold on
        for it=1:3
            x = [1 2 3 5];
            xx = x + (it-2)*0.2;
            [ h(it)  ] = barplot( xx ,y(x,it),z(x,it), col{it} );
            h(it).BarWidth = 0.2;
        end
        
        % legending
        legend([h(1) h(2) h(3)],trtList);
        labels = {'chosen(correct)';'chosen(incorrect)';'unchosen(omitted)';'';'unchosen(premature)'};
        ax = gca; 
        ax.XTick = x;
        ax.XTickLabel = labels(x);
        ax.XTickLabelRotation = 45;
        xlabel('trial classification'); 
        ylabel('number of trials');         
        t = title(['monkey ' sub(1) ]); 
    end
        
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
%% 2) Participation analysis

Participation_analysis;
     
%% 3) Choice
     
Choice_analysis;
    
%% 4) Force

Force_analysis;    
    
%% 5) RT
    
    

% variables
     trt2comp = [2 3];
     vartext = 'responsetime';
     formula = ' responsetime ~ 1 + offervalue + offertime  + dimensionality + congruency';
     varnames = {'offervalue';'offertime';'p_t';'dimensionality';'congruency';'responsetime'};
     effectName = {texlabel('beta_0'),'offervalue','offertime','dimensionality','congruency',};
     nregressor = 5;
     fname= 'identity';  distname = 'normal';
     dataselect = 'selectSubject';
     predictor_select = '[ participation_prediction , offertime , p_t , dimensionality , congruency ];';


     nfactor = numel(varnames)-1;
     nbin = [ nregressor , ntrt , maxSess ];
     nbin2 = [ nregressor ,ntrt,nsub ];
     ysub = nan(nbin);
     Y = nan(nbin2);
     Z = nan(nbin2);
     BETA = []; T = [];

     data.([vartext '_prediction']) = nan(height(data),1);
     
    for isub = 1:nsub
                
         % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        eval(['selection = ' dataselect ';']);

        % variables
        ysub = nan(nbin);
        session = data.sessionNumber(selection);
        participation_prediction = data.participation_prediction(selection);
        participation = selectionParticipate(selection);
        treatment = data.treatment(selection);
        
        r1 = round(data.ordinalRewardLeft(selection));
        r2 = round(data.ordinalRewardRight(selection));
        sum_r = r1+r2;
        max_r = max(r1,r2);
        e1 = round(data.ordinalEffortLeft(selection));
        e2 = round(data.ordinalEffortRight(selection));
        sum_e = e1+e2;
        min_e = min(e1,e2);
        nt = data.trialNumber(selection);
        nnt = data.normalizedTrialNumber(selection);
        ee =  data.cumulativeEffort(selection);
        rr =  data.cumulativeReward(selection);
        
        choice = (data.sideChoice(selection)-1)*2-1;
        sdR = (round(data.deltaOrdinalReward(selection)));
        choiceR = (sign(choice)==sign(sdR));
        dR = abs(sdR);
        sdE = (round(data.deltaOrdinalEffort(selection)));
        choiceE = (sign(choice)==sign(sdE));
        dE = abs(sdE);
        r1 = round(data.ordinalRewardLeft(selection));
        r2 = round(data.ordinalRewardRight(selection));
        e1 = round(data.ordinalEffortLeft(selection));
        e2 = round(data.ordinalEffortRight(selection));
        sum_r = r1+r2;
        sum_e = e1+e2;
        congruency = double(sign(sdR)==-sign(sdE)) ;
        dimensionality = (sdR~=0) + (sdE~=0) ;
        congruency(dimensionality~=2)=1;
        offertime = data.offerTime(selection);
        offertimebin = quantileranks(offertime,nbin(1)); 
        
        responsetime = data.responseTime(selection)+eps;

        
        % previousParticipation
        tmax=1;
        p_t = nan(numel(participation),tmax);
        for it=1:tmax
            p_t(1+it:end,it) = participation(1:end-it);
            p_t(nt<=it,it) = NaN;
        end
        
        cumul = data.stateDuration(selection);
        cumul(participation==0) = 0;
        cumul2 = data.stateDuration(selection);
        cumul2(participation==1) = 0;
        tmax=1;
        work_t = nan(numel(cumul),tmax);
        pause_t = nan(numel(cumul),tmax);
        for it=1:tmax
            work_t(1+it:end,it) = cumul(1:end-it);
            work_t(nt<=it,it) = NaN;
            pause_t(1+it:end,it) = cumul2(1:end-it);
            pause_t(nt<=it,it) = NaN;
        end
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);

        for is = subsess'
            itrt = double(unique(treatment(session==(is))));
            eval(['predictor = ' predictor_select]);
            predictor = (predictor(session==(is),:));
            predictor = nanzscore(predictor);
            eval(['y = ' vartext '(session==is);']);
            stat = fitglm(predictor,y,...
                           formula,'VarNames',varnames,'Distribution',distname,'link',fname);
            coef = stat.Coefficients;
            disp(coef);
            ysub(:,itrt,is) = coef.Estimate;
            data.([vartext '_prediction'])(selection & data.sessionNumber==is) = stat.Fitted.LinearPredictor;
        end

           
        % second-level stat
         dim=3;
         T = [T ; nominal(repmat({'clonidine','placebo','atomoxetine'}',maxSess,1)) ];
         BETA = [BETA ; reshape(ysub,[nregressor,ntrt*maxSess])' ];
         y = nanmean(ysub,dim); Y(:,subtrt,isub) = y;
         z = sem(ysub,dim); Z(:,subtrt,isub) = z;
        
        % display
        if isub==1; fig = figure; set(fig,'Name','glm_participation'); end
        subplot(nsub,1,isub);hold on ; clear h ; h = gobjects(ntrt,1);
        for it=trt2comp
            x = [1:nregressor];
            xx = x + (it-2)*0.25;
            [ h(it)  ] = barplot( xx ,y(:,it),z(:,it), col{it} );
            h(it).BarWidth = 0.25;
            
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            h2leg = []; for it=trt2comp ; h2leg = [h2leg, h(it)];end;
%             legend(h2leg,trtList(trt2comp));
        end
        ax = gca;
%         ax.YLim = [0.5 1];
        ax.TickLength = [0 0];
        ax.XTick = []; 
        ax.XTick = x;
        ax.XTickLabel = effectName(x);
        ax.XTickLabelRotation = 45;
        ax.XLim = [0 x(end)+1];
        ylabel('regression coefficients '); 
        t = title(['monkey ' sub(1) ]); 
        
    end
% group-display
    fig = figure; set(fig,'Name','glm_comp_group');
    hold on ; clear h; h = gobjects(ntrt,1);  dim=3;
    for it=trt2comp
        % inputs
        y = nanmean(Y(:,it,:),dim);  
        z = nanmean(Z(:,it,:),dim);       
        % plot
        x = [1:nregressor];
        xx = x + (it-2)*0.25;
        [ h(it)  ] = barplot( xx ,y,z, col{it} );
        h(it).BarWidth = 0.25;
    end
    
    % legending
%     legend(h2leg,trtList(trt2comp));
    ax = gca;
    ax.TickLength = [0 0];
    ax.XTick = []; 
    ax.XTick = x;
    ax.XTickLabel = effectName(x);
    ax.XTickLabelRotation = 45;
    ax.XLim = [0 x(end)+1];
    ylabel('regression coefficients ');     
    
     % stat
    T = T(~isnan(BETA(:,1)),:);
    predictor = [];
    predictor = [predictor ,(T=='clonidine')];
    predictor = [predictor ,(T=='placebo')];
    predictor = [predictor , (T=='atomoxetine')];
    BETA = BETA(~isnan(BETA(:,1)),:);
    p = nan(nregressor,1);
    for i = 1:nregressor
        contrast = [0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
        model = fitglm(predictor,BETA(:,i),'linear','Intercept',false);
        [p(i),F] = coefTest(model,contrast);

    %         predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),tab.subject,y,...
    %              'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
    %         model2 = fitglme(predictor2,[ vartext ' ~ clonidine + atomoxetine + (clonidine|subject) + (atomoxetine|subject)']);
    %         disp(model2.Coefficients);
    end
    xx = x + (trt2comp(trt2comp~=2)-2)*0.25;
    [s] = sigstar( num2cell(xx),p);
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);  

%% rt

    % variables
        nbin = [ 4 , ntrt , maxSess ];

    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        selection = selectSubject & selectionChosen ;

        % variables
        Y = nan(nbin);
        session = data.session(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection)+0.1;        
        log_rt = log(rt);    
        
        R_options = data.deltaOrdinalReward(selection)~=0 & data.deltaOrdinalEffort(selection)==0  ;
        E_options = data.deltaOrdinalEffort(selection)~=0 & data.deltaOrdinalReward(selection)==0  ;
        congruent_options = selectionCongruent(selection) & data.deltaOrdinalReward(selection)~=0 & data.deltaOrdinalEffort(selection)~=0 ;
        incongruent_options = ~selectionCongruent(selection) & data.deltaOrdinalReward(selection)~=0 & data.deltaOrdinalEffort(selection)~=0 ;

        
        % stats
        var = rt;
%         var = log_rt;
        varcoef  = @(x) nanvar(x)/nanmean(x);
        statistic = @nanmean;
%         statistic = @nanvar;
%         statistic = varcoef;

        
        [~,subtrt] = ismember(unique(treatment),trtList);
        [~,subsess] = ismember(unique(session),sessionList);
        
        ysub = tapply(var,{R_options,treatment,session},statistic);
        Y(1,subtrt,subsess) =  ysub(2,:,:);
        ysub = tapply(var,{E_options,treatment,session},statistic);
        Y(2,subtrt,subsess) =  ysub(2,:,:);
        ysub = tapply(var,{congruent_options,treatment,session},statistic);
        Y(3,subtrt,subsess) =  ysub(2,:,:);
        ysub = tapply(var,{incongruent_options,treatment,session},statistic);
        Y(4,subtrt,subsess) =  ysub(2,:,:);
        
        % averaging
         dim=3;
         y = nanmean(Y,dim);     
         z = sem(Y,dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','rt'); end
        subplot(1,nsub,isub);hold on ; clear h ;
        for it=1:3
            x = [1:nbin(1)];
            xx = x ;
            [~,~,h(it)] = errorscat( xx ,y(x,it), z(x,it),col{it});
        end        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend([h(1) h(2) h(3)],trtList);
        end
        conditionNames = {'1D:reward','1D:effort','2D:congruent','2D:incongruent'};
        ax = gca; 
%         ax.YLim = [0.5 1];
        ax.TickLength = [0 0];
        ax.XTick = x;
        ax.XTickLabel = conditionNames(x);
        ylabel(texlabel('response time (sec)'));
        t = title(['monkey ' sub(1) ]); 
        
    end
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2)
    
%% 5 ) Full Model 

 % variables
     nparam = 13;
     nbin = [ nparam , ntrt , maxSess ];
     Y = nan(nbin);

    for isub = 1:nsub
                
         % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        selection = selectSubject ;

        % variables
        Y = nan(nbin);
        
        session = data.session(selection);
        treatment = data.treatment(selection);
        
        participation = selectionParticipate(selection);
        choice = data.sideChoice(selection)-1;   
        force = data.perf(selection);   
        rt = data.responseTime(selection)+0.1;        
        log_rt = log(rt); 
        
        side = data.sideChoice(selection)-1; 
        side = side*2-1;
        r1 = round(data.ordinalRewardLeft(selection));
        r2 = round(data.ordinalRewardRight(selection));
        sum_r = r1+r2;
        dr = r2-r1;
        e1 = round(data.ordinalEffortLeft(selection));
        e2 = round(data.ordinalEffortRight(selection));
        sum_e = e1+e2;  
        de = e2-e1;
        nt = data.trialNumber(selection);
        nnt = data.normalizedTrialNumber(selection);
        ee =  data.cumulativeEffort(selection);
        rr =  data.cumulativeReward(selection);
        % previous behavior
        tmax=1;
        p_t = nan(numel(participation),tmax);
        choice_t = nan(numel(choice),tmax);
        for it=1:tmax
            p_t(1+it:end,it) = participation(1:end-it);
            p_t(nt<=it,it) = NaN;
            choice_t(1+it:end,it) = choice(1:end-it);
            choice_t(nt<=it,it) = NaN;
        end
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        [~,subsess] = ismember(unique(session),sessionList);

        for is = subsess'
            
            % input/output
            input = ([ r1, r2 , e1 , e2 , nt , side , p_t*2-1, choice_t*2-1 , rr , ee ]);
            input = nanzscore(input(session==sessionList(is),:))';
            observation = ([ participation , choice , force , log_rt ]);
            observation = observation(session==sessionList(is),:)';
            
            % formula
            g_fname = @g_choice_RE;
            
            % priors
            nphi = 13 ;
            param = struct;
            param.prior.mu =    [ 1 , 0  , 1 , 0 , 1  , 0 , 0 , 0 , 0 , ...
                                  0.1 , 0.1 , ...
                                  0 , 1 ];
            param.prior.sigma = [ 1 , 1  , 1 , 1 , 1  , 1 , 1 , 1 , 1 , ...
                                  1 , 1 , ...
                                  1 , 1 ];
            param.type = repmat({'Phi'},1,dim.n_phi);
            param.labels = {'kr','kr2','ke','ke2','kn','krt','ket','kt','k0',...
                            's0','se',...
                            't0','gamma'};
            param.transform.direct = { @identity ,@identity ,@identity ,@identity ,@identity ,@identity ,@identity,@identity ,@identity,...
                                       @safepos ,@safepos ,...
                                       @safepos ,@safepos};
            inG.transform = param.transform.direct(ismember(param.type,'Phi'));
            opt.display=0;
            [priors]  = setParam(param,opt);

            % hyperpriors
            
            % dimensions
            dim = struct('n',0,... % number of hidden states
                    'n_theta',0 ,...   % number of evolution parameters
                    'n_phi',nphi,...     % number of observation parameters
                    'p',4,...          % output (data) dimension
                    'n_t',size(observation,2));   % number of time samples or trials
            options.dim             = dim;
    
            % options
            options.isYout  =  zeros(size(observation,1),size(observation,2));   % data exclusion
            options.isYout([2 3 4],participation==0)   =  1;   
            observation(2,participation==0) = -1 ;
            observation([3 4],participation==0) = 0 ;
            options.priors          = priors;   % include priors in options structure
            options.inG             = inG;      % input structure (grid)
            options.dim             = dim;
            options.DisplayWin      = 0;
            options.verbose         = 1;
            options.updateHP = 1 ;
            options.kernelSize = 1 ;
            options.extended=1;
            options.sources(1).out  = [1,2];    % participation , choice
            options.sources(1).type = 1;        % binomial data
            options.sources(2).out  = 3 ;     % force
            options.sources(2).type = 0 ;      % Normal continous data
            options.sources(3).out  = 4 ;       % response time
            options.sources(3).type = 0;        % Normal continous data
            
            % fit
            [posterior,inversion] = VBA_NLStateSpaceModel(observation,input,[],g_fname,dim,options);
            extract_model_output;
            disp(output.parameters);
            itrt = double(unique(treatment(session==sessionList(is))));
            Y(:,itrt,is) = param.estimate;
        end
           
%         % second-level stat
%          dim=3;
%          y = nanmean(Y,dim);
%          z = sem(Y,dim);
%         
%         % display
%         if isub==1; fig = figure; set(fig,'Name','glm_participation'); end
%         subplot(nsub,1,isub);hold on
%         for it=1:3
%             x = [1:nregressor];
%             xx = x + (it-2)*0.25;
%             [ h(it)  ] = barplot( xx ,y(:,it),z(:,it), col{it} );
%             h(it).BarWidth = 0.25;
%             
% %             h = plotSpread(Y','distributionColors',col);
% %             m = findobj('-property','MarkerFaceColor');
% %             set(m,'MarkerSize',18);
% %             xx = double(unique(sessionList));
% %             h(it) = scatter(xx,Y(it,:),40,col{it},'filled');
%             
%         end
%         
%         % legending
%         if isa(h,'matlab.graphics.Graphics') && isub==nsub
%             legend([h(1) h(2) h(3)],trtList);
%         end
%         effectName = varnames(1:end-1);
%         effectName = {'k_0', effectName{:}};
%         ax = gca;
% %         ax.YLim = [0.5 1];
%         ax.TickLength = [0 0];
%         ax.XTick = []; 
%         ax.XTick = x;
%         ax.XTickLabel = effectName(x);
%         ax.XLim = [0 x(end)+1];
%         ylabel('regression coefficients '); 
%         t = title(['monkey ' sub(1) ]); 
        
    end

%     % reformat
%     setFigProper('FontSize',20,'LineWidth',2);


%% 6) Ethogram



     % variables
    trt2comp = [2 3];
    ytext = 'number';
    vartext = 'diversity';
    xtext = vartext;
    aggregate = 1;
    nbehavior = 12;
    nbin = [ nbehavior, ntrt , maxSess ];
    if aggregate 
        nbin = [ ntrt , maxSess ]; 
        nbehavior = 1;
    end
    nbin2 = [ nbehavior, ntrt , nsub ];
    behaviorList = getlevels(ethogram.behavior);
    xlist =    {'AFFILIATIVE','AGGRESSIVE','ALERT','AVOIDANCE','FOOD CONSUMPTION','FOOD SEARCH',...
                'OBJECT MANIPULATION','STEREOTYPING','RESTING','SELF INJURING','SELF PROTECTING','SUBMISSION'};

    Y = nan(nbin2);
    Z = nan(nbin2);
    BETA = []; T = [];
        
    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(ethogram.subject, sub);  
        selection = selectSubject   ;

        % variables
        ysub = nan(nbin);
        session = ethogram.sessionNumber(selection);
        treatment = ethogram.treatment(selection);
        behavior = ethogram.behavior(selection);
        time = ethogram.time(selection);
        valence = double( behavior=='AFFILIATIVE' |...
                    behavior=='FOOD_CONSUMPTION' |...
                    behavior=='OBJECT_MANIPULATION' |...
                    behavior=='SELF_PROTECTION' |...
                    behavior=='SUBMISSION' ) ;
        valence(valence==0) = -1;
        valence(behavior=='RESTING') = 0;

        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);
        
        for is = subsess'

            % metric
            freq = histcounts(double(behavior(session==(is))),'Normalization','Probability','BinLimits',[1,nbehavior]);
            totaltime = tapply(time(session==(is)),{behavior(session==(is))},@nansum);
            diversity = numel(unique(behavior(session==(is))));
            stability = nanmean(time(session==(is)));
            positivity  = nanmean(time( session==(is) & valence==1 ));
            activity = mean(~(behavior(session==(is))=='RESTING'));
            
            ib = double(unique(behavior(session==(is))));
            itrt = double(unique(treatment(session==(is))));
%             Y(ib,itrt,is) = totaltime;  
            if aggregate 
                eval(['vary = ' vartext ';']);
                ysub(itrt,is) = vary;
            else
            ysub(:,itrt,is) = freq;     
            end     
        end


        % averaging
         dim=3;
         if aggregate; dim=2; end
         y = nanmean(ysub,dim);     Y(:,subtrt,isub) = y;
         z = sem(ysub,dim);         Z(:,subtrt,isub) = z;
         T = [T ; nominal(repmat({'clonidine','placebo','atomoxetine'}',maxSess,1)) ];
         BETA = [BETA ; reshape(ysub,[nbehavior,ntrt*maxSess])' ];
        
        % display
        if isub==1; fig = figure; set(fig,'Name','ethogram'); end
        subplot(1,nsub,isub);hold on ; clear h ;
        for it=trt2comp
            x = [1:nbehavior];
            if aggregate; x = [1]; end
            xx = x + (it-2)*0.25;
            if aggregate
                [ h(it)  ] = barplot( xx ,y(it),z(it), col{it} );
            else
                [ h(it)  ] = barplot( xx ,y(:,it),z(:,it), col{it} );
            end
            h(it).BarWidth = 0.15;
            
        end
        
        % legending
%         if isa(h,'matlab.graphics.Graphics') && isub==nsub
%             legend([h(1) h(2) h(3)],trtList);
%         end
        ax = gca; 
        ax.TickLength = [0 0];
        if ~aggregate
            ax.XTick = x;
            ax.XTickLabel = cellstr(xlist);
        end
        ax.XTickLabelRotation = 45; 
        ylabel(ytext); 
        if aggregate; xlabel(xtext);  end
        t = title(['monkey ' sub(1) ]); 
        
    end
    
     % group-display
    fig = figure; set(fig,'Name','bivar_comp_group');
    hold on ; clear h;  dim=3;
    for it=trt2comp
        % inputs
        y = nanmean(Y(:,it,:),dim);  
        z = nanmean(Z(:,it,:),dim);       
        % plot
        x = [1:nbehavior];
        xx = x + (it-2)*0.25;
        [ h(it)  ] = barplot( xx ,y,z, col{it} );
        h(it).BarWidth = 0.15;
    end
    
    % legending
%     legend(h2leg,trtList{trt2comp});
    ax = gca;
    ax.TickLength = [0 0];
    if ~aggregate
        ax.XTick = x;
        ax.XTickLabel = cellstr(behaviorList);
    else
        ax.XTick = [];
    end
    ax.XTickLabelRotation = 45; 
    ylabel(ytext); 
    if aggregate; xlabel(xtext);  end
    
    % stat
    T = T(~isnan(BETA(:,1)),:);
    predictor = [];
    predictor = [predictor ,(T=='clonidine')];
    predictor = [predictor ,(T=='placebo')];
    predictor = [predictor , (T=='atomoxetine')];
    BETA = BETA(~isnan(BETA(:,1)),:);
    p = nan(nbehavior,1);
    for i = 1:nbehavior
        contrast = [0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
        model = fitglm(predictor,BETA(:,i),'linear','Intercept',false);
        [p(i),F] = coefTest(model,contrast);

    %         predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),tab.subject,y,...
    %              'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
    %         model2 = fitglme(predictor2,[ vartext ' ~ clonidine + atomoxetine + (clonidine|subject) + (atomoxetine|subject)']);
    %         disp(model2.Coefficients);
    end
    xx = x + (trt2comp(trt2comp~=2)-2)*0.25;
%     [s] = sigstar( num2cell(xx),p);

    % reformat
    setFigProper('FontSize',20,'LineWidth',2);

