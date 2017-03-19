function [dataset,stat] = do_glm_monkeydrug(designtab,dataset,yvar,formula,paramNames,fname,dataselection,trt2comp,modelname,glmtype)  

% variables
 vartext = yvar;
 nregressor = numel(paramNames);
 switch fname
     case 'identity'
     distname = 'normal';
     case 'logit'
     distname = 'binomial';
 end
     
 % parameters
 %%% display option
 col = {[0.54 0.27 0.07],...
    [1 1 1]*0.5,...
    [1 0.5 0]};
 %%% variables selection
 varnames = {'sumR';'sumE';'dR';'dE';'rchoice';'echoice';
            'side';'ntrial';'Rtot';'Etot';
            'part_t';'choice_t';yvar};
%%% design information
 subjectList = {'aliosha','bob','esmeralda'};
 nsub = numel(subjectList);
 trtList = {'clonidine','placebo','atomoxetine'};
 ntrt = numel(trtList);
 sessionList = unique(dataset.session);
 maxSess = max(designtab.sessionNumber);
%%% data selection pointers
selectionAll =  ones(height(dataset),1) ;
selectionRepetition = (dataset.isRepeatedTrial==1) ;
selectionMissed = ( dataset.errorType==1);
selectionPremature = ( dataset.errorType==2);
selectionSwitch = ( dataset.errorType==3);
selectionIncorrect = ( dataset.errorType==4);
selectionParticipate = ~selectionMissed ;
selectionResponseTime =  ~selectionMissed & ~isnan(dataset.intertrialTime) ;
selectionChosen = ( ~selectionMissed & ~selectionPremature & ~selectionSwitch  );
selectionCorrect = ( selectionChosen  & ~selectionIncorrect );
dR =  round(dataset.deltaOrdinalReward);
dE =  round(dataset.deltaOrdinalEffort);
selectionCongruent = (sign(dR)~=sign(dE));
selectionCongruent(sign(dR)==0 & sign(dE)==0) = 0;
     

% data preparation
 nfactor = numel(varnames)-1;
 nbin = [ nregressor , ntrt , maxSess ];
 nbin2 = [ nregressor ,ntrt,nsub ];
 ysub = nan(nbin);
 Y = nan(nbin2);
 Z = nan(nbin2);
 BETA = []; T = [];
 accuracy = [];
 predictor = [];
 dataset.([vartext '_prediction']) = nan(height(dataset),1);
     
    for isub = 1:nsub
                
         % select
        sub = subjectList{isub};
        selectSubject = ismember(dataset.subject, sub);   
        eval(['selection = ' dataselection ';']);

        % variables
        ysub = nan(nbin);
        session = dataset.sessionNumber(selection);
        participation = selectionParticipate(selection);
%         participation_prediction = dataset.participation_prediction(selection);
        treatment = dataset.treatment(selection);
        
        r1 = round(dataset.ordinalRewardLeft(selection));
        r2 = round(dataset.ordinalRewardRight(selection));
        sumR = r1+r2;
        max_r = max(r1,r2);
        e1 = round(dataset.ordinalEffortLeft(selection));
        e2 = round(dataset.ordinalEffortRight(selection));
        sumE = e1+e2;
        min_e = min(e1,e2);
        ntrial = dataset.trialNumber(selection);
        nnt = dataset.normalizedTrialNumber(selection);
        Etot =  dataset.cumulativeEffort(selection);
        [b,~,stat] = glmfit(nnt,Etot,'normal'); % orthogonalize
        Etot = stat.resid + b(1);
        Rtot =  dataset.cumulativeReward(selection);
        [b,~,stat] = glmfit(nnt,Rtot,'normal');
        Rtot = stat.resid + b(1);
        rchoice = dataset.chosenCardinalReward(selection);
        echoice = dataset.chosenCardinalEffort(selection);
        choice = dataset.sideChoice(selection)-1;   
        side = (choice)*2-1;
        dR = (round(dataset.deltaOrdinalReward(selection)));
        choiceR = (sign(choice)==sign(dR));
        dE = (round(dataset.deltaOrdinalEffort(selection)));
        choiceE = (sign(choice)==sign(dE));
        congruency = double(sign(dR)==-sign(dE)) ;
        dimensionality = (dR~=0) + (dE~=0) ;
        congruency(dimensionality~=2)=1;
        offertime = dataset.offerTime(selection);
        offertimebin = quantileranks(offertime,nbin(1)); 
        responsetime = dataset.responseTime(selection)+eps;
        force = dataset.perf(selection);   
        cum_part = nan(size(participation));
        % previousParticipation
        tmax=1;
        part_t = nan(numel(participation),tmax);
        for it=1:tmax
            part_t(1+it:end,it) = participation(1:end-it);
            part_t(ntrial<=it,it) = NaN;
        end
        cumul = dataset.stateDuration(selection);
        cumul(participation==0) = 0;
        cumul2 = dataset.stateDuration(selection);
        cumul2(participation==1) = 0;
        tmax=1;
        work_t = nan(numel(cumul),tmax);
        pause_t = nan(numel(cumul),tmax);
        for it=1:tmax
            work_t(1+it:end,it) = cumul(1:end-it);
            work_t(ntrial<=it,it) = NaN;
            pause_t(1+it:end,it) = cumul2(1:end-it);
            pause_t(ntrial<=it,it) = NaN;
        end
        % previousChoice
        tmax=1;
        choice_t = nan(numel(choice),tmax);
        for it=1:tmax
            choice_t(1+it:end,it) = choice(1:end-it);
            choice_t(ntrial<=it,it) = NaN;
        end
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);

        for is = subsess'
            cum_part(session==(is)) = cumsum(participation(session==(is)));
            
            itrt = double(unique(treatment(session==(is))));
            predictor=[];
            for iv = 1:(numel(varnames)-1)
                eval([ ' predictor = [predictor , ' varnames{iv} '];']);
            end
            predictor = (predictor(session==(is),:));
            predictor = nanzscore(predictor);
            eval(['y = ' vartext '(session==is);']);
            
            % fit
            switch glmtype
                case 'MLH' % maximum-likelihood
                    stat = fitglm(predictor,y,...
                                   formula,'VarNames',varnames','Distribution',distname,'link',fname);
                   % extract coefficient
                    coef = stat.Coefficients;
                    disp(coef);
                    beta = coef.Estimate;
                    % extract prediction
                    prediction = stat.Fitted.LinearPredictor;
                    % extract model accuracy
                    switch fname
                        case 'identity'
                            accuracy = [accuracy; mean(stat.Rsquared.Ordinary) ];
                        case 'logit'
                            acc1 =  mean(stat.Fitted.Response(stat.Variables.(vartext)==0)<=0.5);
                            acc2 = mean(stat.Fitted.Response(stat.Variables.(vartext)==1)>0.5);
                            accuracy = [accuracy; mean([acc1 acc2]) ];
                    end
                case 'LASSO'
                case 'MAP' % maximum-a-posteriori
                    % construct predictor matrix
                    stat = fitglm(predictor,y,...
                                   formula,'VarNames',varnames','Distribution',distname,'link',fname);
                    terms = stat.Formula.Terms(:,1:end-1);
                    predictormat = [];
                    for iterm = 1:size(terms,1)
                        ind = find(terms(iterm,:)==1);
                        predictorvec = prod(predictor(:,ind),2);
                        predictorvec = predictorvec + (predictorvec==0);
                        predictormat = [predictormat , predictorvec  ];
                    end
                    % fit
                    [beta,L,posterior,out] = nanglm(predictormat,y,'logit',1);
                    disp(beta);
                    accuracy = [accuracy; out.fit.acc];
                    prediction = out.suffStat.gx;
            end
            ysub(:,itrt,is) = beta;
            dataset.([vartext '_prediction'])(selection & dataset.sessionNumber==is) = prediction;
        end



           
        % second-level stat
         dim=3;
         T = [T ; nominal(repmat({'clonidine','placebo','atomoxetine'}',maxSess,1)) ];
         BETA = [BETA ; reshape(ysub,[nregressor,ntrt*maxSess])' ];
         y = nanmean(ysub,dim); Y(:,subtrt,isub) = y;
         z = sem(ysub,dim); Z(:,subtrt,isub) = z;
         
        % display
        if isub==1; fig = figure; set(fig,'Name','glm_parameters_subjects'); end
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
        ax.XTickLabel = paramNames(x);
        ax.XTickLabelRotation = 45;
        ax.XLim = [0 x(end)+1];
        ylabel('regression coefficients '); 
        t = title(['monkey ' sub(1) ]); 
         
         
    end
    
    % group-display
    fig = figure; set(fig,'Name','glm_parameters_subjects');
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
    ax.XTickLabel = paramNames(x);
    ax.XTickLabelRotation = 45;
    ax.XLim = [0 x(end)+1];
    ylabel('regression coefficients ');  
    
    
    % stat
    T = T(~isnan(BETA(:,1)),:);
    predictor = [];
    predictor = [predictor ,(T=='clonidine')];
    predictor = [predictor ,(T=='placebo')];
    predictor = [predictor , (T=='atomoxetine')];
    subject = designtab.subject;
    BETA = BETA(~isnan(BETA(:,1)),:);
    stat=struct;
    stat.p = nan(1,nregressor);    stat.F = nan(1,nregressor);
    stat.pRFX = nan(3,nregressor); stat.F_RFX = nan(3,nregressor); 
    for i = 1:nregressor
        contrast = [0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
        y = BETA(:,i);

        switch modelname
            case 'FFX'
                % fixed-effect model
                model = fitglm(predictor,y,'linear','Intercept',false);
                [p(i),F(i)] = coefTest(model,contrast);
                [~,~,fixedstat] = fixedEffects(model);
                nffx = numel(fixedstat.Estimate);
                fixedstat.Name = repmat(paramNames(i),nffx,1) ;
            case 'RFX'
                % random-effect model
                if numel(trt2comp)>1
                    %%% treatment effect
                    predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),subject,y,...
                         'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
                    model = fitglme(predictor2,[ vartext '  ~ -1 + clonidine + placebo +  atomoxetine + ( -1 + clonidine + placebo +  atomoxetine | subject)']);
                    [p(i),F(i)] = coefTest(model,contrast);
                   [~,~,randomstat] = randomEffects(model);
                    randomstat.Name = repmat(paramNames(i),3,1) ;
                    for isub=1:3
                        contrastsub = zeros(1,9);
                        contrastsub([1:3]+(isub-1)*3) = contrast;
                       [stat.pRFX(isub,i),stat.F_RFX(isub,i)] = coefTest(model,contrast,0,'REContrast',contrastsub);
                    end

                else
                    %%% average effect
                    predictor2 = table(subject,y,...
                         'VariableNames',{'subject',vartext});
                    model = fitglme(predictor2,[ vartext ' ~ 1 + (1|subject)']);
                   [p(i),F(i)] = model.Coefficients.pValue;
                   [~,~,randomstat] = randomEffects(model);
                    randomstat.Name = repmat(paramNames(i),3,1) ;
                end
                disp(model.Coefficients);
                [~,~,fixedstat] = fixedEffects(model);
                nffx = numel(fixedstat.Estimate);
                fixedstat.Name = repmat(paramNames(i),nffx,1) ;
                if i==1
                    coef = fixedstat;
                    if exist('randomstat'); rcoef = randomstat;end
                else
                    coef = vertcat(coef,fixedstat);
                    if exist('randomstat');rcoef = vertcat(rcoef,randomstat);end
                end

        end
    end

    
    % accuracy comparison
    if numel(trt2comp)>1
        y = accuracy;
        predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),subject,y,...
             'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
        model = fitglme(predictor2,[ vartext '  ~ -1 + clonidine + placebo +  atomoxetine + ( -1 + clonidine + placebo +  atomoxetine | subject)']);
       [accup,accuF] = coefTest(model,contrast);
        [~,~,accutab] = fixedEffects(model);
        nffx = numel(accutab.Estimate);
        accutab.Name = repmat('accuracy',nffx,1) ;
    end
    
switch modelname
    case 'MANOVA'
        y = BETA; group = T;
        y = y(ismember(T,trtList(trt2comp)),:);
        group = group(ismember(T,trtList(trt2comp)));
        [stat.d,stat.p,manovastat] = manova1(y,group);
        stat.eigenvec = manovastat.eigenvec;
    otherwise
        if numel(trt2comp)>1

        else
            coef = dataset2table(coef);
            if exist('randomstat')
                rcoef = dataset2table(rcoef);
                firstletter = @(s) s(1) ;
                rcoef.Level = cellfun(firstletter,rcoef.Level);
            end
        end

        if numel(trt2comp)>1
            xx = x + (trt2comp(trt2comp~=2)-2)*0.25; % select treatment except placebo
        else
            xx = x + (trt2comp-2)*0.25; % select treatment
        end
        [s] = sigstar( num2cell(xx),p);

        % reformat
        setFigProper('FontSize',20,'LineWidth',2);    


         % stores
         stat.coef = coef;
         stat.accuracy.coef = accutab;
         stat.accuracy.F = accuF;
         stat.accuracy.p = accup;
         stat.paramNames = paramNames;
         stat.parameters = BETA;
         stat.p = p;
         stat.F = F;
         if exist('randomstat'); stat.rcoef = rcoef;end
end
     
end

         
        
