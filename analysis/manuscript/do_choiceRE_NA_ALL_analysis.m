%% do_choiceRE_NA_ALL_analysis

clc;
clear all;
close all;

%%
datadir = 'B:\nicolas.borderies\projets\monkey_pharma_choice\data\dataset_behavior_ethogram';
codedir = 'B:\nicolas.borderies\projets\monkey_pharma_choice\code\analysis';
cd(datadir);
analysisname = ['monkeydrug_NA_dataset.mat'];
load(analysisname);
cd(codedir);


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
        
%% reformat
    ind = find(ismember(data.Properties.VariableNames,'gain'));
    data.Properties.VariableNames{ind} = 'ordinalRewardOutcome';
        
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


%% 0.1) FLexible plots to compare dependant variables btw treatments
     
%  compare dependant variables btw treatments
[stat,model]  = barcomp_monkeydrug(subjectList,sessiondata,...
                            'participation_rate','participation rate (%)',[1 2],'RFX');


    
%% 0.1.2) FLexible plot to correlate cross-session measures
     
[stat] = scattercomp_monkeydrug(design,sessiondata,'success_force','correct force (%)','choice_E','choice = HE (%)',...
                                                      'betweenConditionConfound',...
                                                      [0 1],[0 1],[1 2]); 
    
    
%% 0.2) FLexible plot to compare dependant variables btw treatments

stat = linecomp_monkeydrug(design,data,'peak_force','peak force(%)','peak_velocity','peak force(%)','peak_force_bin','glm',...
                            10,0,[0 1],[],'selectSubject & selectionChosen',[1 2]);   
    
%% 0.3) FLexible plot to compare GLM parameters between treatments

formula = 'force ~ 1 + side + side:echoice + rchoice*echoice + ntrial';
effectName = {'k_0','rchoice','echoice','side','ntrial','rchoice*echoice ','echoice*side'};
[data,stat] =  do_glm_monkeydrug(design,data,'force',...
                                formula,effectName,...
                                'identity','selectSubject & selectionChosen',[1 2],'RFX');


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
   
Responsetime_analysis;
    
%% 5 ) Complete computational model (participation + choice + force)

 % variables
     nparam = 11;
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
            
            % dimensions
            nphi = nparam ;
            dim = struct('n',0,... % number of hidden states
                    'n_theta',0 ,...   % number of evolution parameters
                    'n_phi',nphi,...     % number of observation parameters
                    'p',4,...          % output (data) dimension
                    'n_t',size(observation,2));   % number of time samples or trials
            options.dim             = dim;
            
            % priors
            param = struct;
            param.prior.mu =    [ 0.3 , -0.1  , 0.1 , 0.5 , 1 ,...
                                  0.8 , 0 , 0.5 , ...
                                  0.2 ,...
                                  0 , 1 ];
            param.prior.sigma = [ 1 , 1 , 1 , 1 , 1 , ...
                                  1 , 1 , 1 , ...
                                  1 , ...
                                  0 , 0 ];
            param.type = repmat({'Phi'},1,dim.n_phi);
            param.labels = {'kr','ks','ke','kf','k0',...
                            'beta','bm','bp',...
                            'se',...
                            't0','gamma'};
            param.transform.direct = { @identity ,@identity ,@identity ,@identity ,@identity ,...
                                       @safepos ,@identity,@identity ,...
                                       @safepos ,...
                                       @safepos ,@safepos};
            inG.transform = param.transform.direct(ismember(param.type,'Phi'));
            opt = struct('display',0);
            [priors]  = setParam(param,opt);

            % hyperpriors
            
    
            % options
            %%% data exclusion
            options.isYout  =  zeros(size(observation,1),size(observation,2));   
            options.isYout([2 3 4],participation(session==sessionList(is))==0)   =  1;    % unparticipated trials
            options.isYout(3,force(session==sessionList(is))==0)   =  1;    % performance not available
            options.isYout(4,:)   =  1;    % response-time
            observation(options.isYout==1) = 0 ; % excluded data set to zero

            input([7,8],1) = 0; % set to zero previous participation & choice for first trial
            
            options.priors          = priors;   % include priors in options structure
            options.inG             = inG;      % input structure (grid)
            options.dim             = dim;
            options.DisplayWin      = 1;
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

Ethogram_analysis;


