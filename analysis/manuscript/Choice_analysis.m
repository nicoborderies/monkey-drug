%% Choice_analysis
%   this script execute the following analysis:
%    - display between-treatments comparisons of choices (HR/HE/ correct choice)
%    - display between-treatments comparisons of choices dependance
%      to experimental factors 
%    - glmm modeling of treatment effect on choice
%    - between-treatments comparisons of glmm coefficients
%    - model selection between different glmm models
%

% Requirements:
%   Script: do_choiceRE_NA_ALL_analysis
%   Subfunctions: 
%   Data-files: monkeydrug_NA_dataset
%   Matlab-version:
%   Matlab-toolbox: 
%
% See also: 
%
% Author: Nicolas Borderies
% email address: nico.borderies@gmail.com 
% February 2017; Last revision: 


%% CLONIDINE-PLACEBO comparisons
% -------------------------------

% between-treatments comparisons of choices (HR/HE/ correct choice)
% choice HE
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'accuracy_E','choice = HE (%)',[1 2],'RFX');
% choice HR
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'accuracy_R','choice = HR (%)',[1 2],'RFX');
                            
% choice HRHE
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'choice_HRHE','choice = HRHE (%)',[1 2],'RFX');
        
% participation
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'participation_rate','participation(%)',[1 2],'RFX');
    
%% GLM models between treatments
% -------------------------------

% logistic regression of choice + MFX comparison of treatments
formula = 'choice ~ 1 + dR + dE + dR:ntrial + dE:ntrial + choice_t';
effectName = {texlabel('b_right'),'dR','dE','choice_t','dR*ntrial','dE*ntrial'};
[data,stat] =  do_glm_monkeydrug(design,data,'choice',...
                                formula,effectName,...
                                'logit','selectSubject & selectionChosen',[1 2],'RFX');
                            
% manova on regression coefficients by treatment
[data,stat] =  do_glm_monkeydrug(design,data,'choice',...
                                formula,effectName,...
                                'logit','selectSubject & selectionChosen',[1 2],'MANOVA');
     
%% option value
%     factor = 'R';
    factor = 'E';
    trt2comp = [1 2];
    ytext = 'choice = HE (%)';
    ylimits = [0 0.5];

    % variables
        nbin = [ 4 , ntrt , maxSess ];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        selection = selectSubject & selectionChosen ;
        
        % variables
        Y = nan(nbin);
        
        session = data.sessionNumber(selection);
        treatment = data.treatment(selection);
        
        participation = selectionParticipate(selection);
        choice = (data.sideChoice(selection)-1)*2-1;
        dR = (round(data.deltaOrdinalReward(selection)));
        choiceR = (sign(choice)==sign(dR));
        dR = abs(dR);
        dE = (round(data.deltaOrdinalEffort(selection)));
        choiceE = (sign(choice)==sign(dE));
        dE = abs(dE);
        accuracy = (choiceR | choiceE) ; 
            
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);
        
            % mean participation
            if factor=='R'; varx = dR; vary = choiceR; 
            else  varx = dE;  vary = choiceE;  end
            ysub = tapply(vary,{varx,treatment,session},@nanmean);
            Y(:,subtrt,subsess) =  ysub;
            

        % second-level stat
         dim=3;
         y = nanmean(Y,dim);
         z = sem(Y,dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','choice'); end
        subplot(1,nsub,isub);hold on ; clear h
        for it=trt2comp
            x = [2:nbin(1)];
            xx = x ;
            [~,~,h(it)] = errorscat( xx ,y(x,it), z(x,it),col{it});
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend(trtList{trt2comp});
        end
        ax = gca;
        if factor=='R';ax.YLim = [0.5 1];
        elseif factor=='E';ax.YLim = [0 0.5];end
        ax.TickLength = [0 0];
        ax.XTick = []; 
        if factor=='R'; xlabel(texlabel('abs(Delta R)'));
        else            xlabel(texlabel('abs(Delta E)'));  end
        if isub==1 ; 
            if factor=='R'; ylabel('choice = HR (%)'); 
            else            ylabel('choice = HE (%)');   end
        end
        t = title(['monkey ' sub(1) ]); 
        
    end

    % reformat
    setFigProper('FontSize',20,'LineWidth',2);

    
%% option value
    % variables
        nbin = [ 6 , ntrt , maxSess ];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        selection = selectSubject & selectionChosen  & data.deltaOrdinalReward==0 & data.deltaOrdinalEffort~=0;
        selection = selectSubject & selectionChosen ;
 
        
        % variables
        X = nan(nbin);
        Y = nan(nbin);
        session = data.session(selection);
        treatment = data.treatment(selection);
        
        choice = (data.sideChoice(selection)-1)*2-1;
        sumR = round(data.ordinalRewardLeft(selection) + data.ordinalRewardRight(selection));
        dE = (round(data.deltaOrdinalEffort(selection)));
        choiceE = (sign(choice)==-sign(dE));
        dE = abs(dE);
        successRate = data.successRate(selection);
        successRate2 = quantileranks(successRate,nbin(1));
        successRate2(successRate2==0) = NaN;
        nnt = data.normalizedTrialNumber(selection); 
        nnt = quantileranks(nnt,nbin(1));

        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        [~,subsess] = ismember(unique(session),sessionList);
        
        varx = nnt;  vary = choiceE;
        xsub = tapply(varx,{varx,treatment,session},@nanmean);
        ysub = tapply(vary,{varx,treatment,session},@nanmean);
        X(:,subtrt,subsess) =  xsub;
        Y(:,subtrt,subsess) =  ysub;
            

        % second-level stat
         dim=3;
%          xlevel = unique(varx);
         xlevel = nanmean(X,dim);
         y = nanmean(Y,dim);
         z = sem(Y,dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','choice'); end
        subplot(1,nsub,isub);hold on ; clear h
        for it=1:3
%         for it=1:2
            x = [1:nbin(1)];
            xx = xlevel(x) ;
            [~,~,h(it)] = errorscat( xx ,y(x,it), z(x,it),col{it});
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend([h(1) h(2) h(3)],trtList);
        end
        ax = gca;
%         ax.XLim = [0 max(xx)*(1.1)];
        ax.YLim = [0 1];
%         ax.TickLength = [0 0];
        ax.XTick = sort(xx(~isnan(xx))) ; 
%         xlabel(texlabel('R_left + R_right'));
        xlabel(texlabel('success rate (%)'));
        if isub==1 ; 
          ylabel('choice = min(E) (%)'); 
        end
        t = title(['monkey ' sub(1) ]); 
        
    end

    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
%% GLMM models
 % variables
     trt2comp = [2];
     vartext = 'choice';
     formula = ' choice ~ 1 + dR + dE + dR:ntrial + dE:ntrial + choice_t_1 ';
     varnames = {'dR';'R2';'dE';'E2';'ntrial';'choice_t_1';'Rtot';'Etot';'choice'};
     effectName = {texlabel('b_right'),'dR','dE','choice_t_1','dR*ntrial','dE*ntrial'};
     nregressor = 6;
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
        selection = selectSubject & selectionChosen ;
        
        % variables
        ysub = nan(nbin);
        session = data.sessionNumber(selection);
        choice = data.sideChoice(selection)-1;   
        treatment = data.treatment(selection);
        
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
        % previousChoice
        tmax=1;
        choice_t = nan(numel(choice),tmax);
        for it=1:tmax
            choice_t(1+it:end,it) = choice(1:end-it);
            choice_t(nt<=it,it) = NaN;
        end
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);

        for is = subsess'
            
            itrt = double(unique(treatment(session==is)));
            predictor = [ dr, r2 , de , e2 , nt , choice_t*2-1 , rr , ee ];
            predictor = (predictor(session==(is),:));
            predictor = nanzscore(predictor);
            y = choice(session==(is));
            stat = fitglm(predictor,y,...
                           formula,'VarNames',varnames,'Distribution','binomial','link','logit');
            coef = stat.Coefficients;
            disp(coef);
            ysub(:,itrt,is) = coef.Estimate;
            data.([vartext '_prediction'])(selection & data.sessionNumber==is) = stat.Fitted.LinearPredictor;

            %             x = [ predictor(:,[1:4 6]) , predictor(:,1).*predictor(:,5) , predictor(:,3).*predictor(:,5) ];
%             [B,fit] = lasso(x,y);
%             Y(1,itrt,is) = fit.Intercept(1);
%             Y(2:end,itrt,is) =  B(:,1);

            
        end
           
        % second-level stat
         dim=3;
         T = [T ; nominal(repmat({'clonidine','placebo','atomoxetine'}',maxSess,1)) ];
         BETA = [BETA ; reshape(ysub,[nregressor,ntrt*maxSess])' ];
         y = nanmean(ysub,dim); Y(:,subtrt,isub) = y;
         z = sem(ysub,dim); Z(:,subtrt,isub) = z;
        
        % display
        if isub==1; fig = figure; set(fig,'Name','glm_choice'); end
        subplot(nsub,1,isub);hold on ; clear h ; h = gobjects(ntrt,1);
        for it=1:trt2comp
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
    subject = sessiondata.subject;
    BETA = BETA(~isnan(BETA(:,1)),:);
    p = nan(nregressor,1);
    for i = 1:nregressor
        contrast = [0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
        y = BETA(:,i);

        % fixed-effect model
%         model = fitglm(predictor,y,'linear','Intercept',false);
%         [p(i),F] = coefTest(model,contrast);
        % random-effect model
            % treatment effect
%             predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),subject,y,...
%                  'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
%             model2 = fitglme(predictor2,[ vartext ' ~ clonidine + atomoxetine + (clonidine|subject) + (atomoxetine|subject)']);
%             disp(model2.Coefficients);
            % average effect
            predictor2 = table(subject,y,...
                 'VariableNames',{'subject',vartext});
            model2 = fitglme(predictor2,[ vartext ' ~ 1 + (1|subject)']);
            [~,~,randomstat] = randomEffects(model2);
            randomstat.Name = repmat(effectName(i),3,1) ;
            disp(model2.Coefficients);
            if i==1
                coef = model2.Coefficients;
                rcoef = randomstat;
            else
                coef = vertcat(coef,model2.Coefficients);
                rcoef = vertcat(rcoef,randomstat);
            end
            p(i) = model2.Coefficients.pValue;
        
    end
    coef.Name = effectName';
    coef = dataset2table(coef);
    rcoef = dataset2table(rcoef);
    firstletter = @(s) s(1) ;
    rcoef.Level = cellfun(firstletter,rcoef.Level);
%     xx = x + (trt2comp(trt2comp~=2)-2)*0.25; % select treatment except placebo
    xx = x + (trt2comp-2)*0.25; % select treatment
    [s] = sigstar( num2cell(xx),p);
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);    
    
%% GLMM model selection
% variables
     trt2comp = [2];
     vartext = 'choice';
     formula = ' choice ~ 1 + dR + dE + dR:ntrial + dE:ntrial + choice_t_1 ';
     varnames = {'dR';'R2';'dE';'E2';'ntrial';'choice_t_1';'Rtot';'Etot';'choice'};
     effectName = {texlabel('b_right'),'dR','dE','choice_t_1','dR*ntrial','dE*ntrial'};
     nregressor = 6;
     fname= 'logit';  distname = 'binomial';
     dataselect = 'selectSubject & selectionChosen';
     predictor_list = {'dr';'r2';'de';'e2';'nt';
                        'nt.*dr';'nt.*de';'choice_t';'rr';'ee';
                        'dPsuccess'};
     npredictor = numel(predictor_list);
     model_matrix = {[1 3 6 7 8] ; 
                     [1 6 8 11 ]};
     nmodel=size(model_matrix,1);

     nfactor = numel(varnames)-1;
     nbin = [ nregressor , ntrt , maxSess ];
     nbin2 = [ nregressor ,ntrt,nsub ];
     ysub = nan(nbin);
     Y = nan(nbin2);
     Z = nan(nbin2);
     BETA = []; T = [];
     LPP = [];
     

     data.([vartext '_prediction']) = nan(height(data),1);
     
    for isub = 1:nsub
                
         % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        eval(['selection = ' dataselect ';']);

        % variables
          ysub = nan(nbin);
        session = data.sessionNumber(selection);
        choice = data.sideChoice(selection)-1;   
        treatment = data.treatment(selection);
        
        r1 = round(data.ordinalRewardLeft(selection));
        r2 = round(data.ordinalRewardRight(selection));
        sum_r = r1+r2;
        dr = r2-r1;
        e1 = round(data.ordinalEffortLeft(selection));
        e2 = round(data.ordinalEffortRight(selection));
        sum_e = e1+e2;
        de = e2-e1;
        dPsuccess = data.successRate_right(selection)-data.successRate_left(selection);
        nt = data.trialNumber(selection);
        nnt = data.normalizedTrialNumber(selection);
        ee =  data.cumulativeEffort(selection);
        rr =  data.cumulativeReward(selection);
        % previousChoice
        tmax=1;
        choice_t = nan(numel(choice),tmax);
        for it=1:tmax
            choice_t(1+it:end,it) = choice(1:end-it);
            choice_t(nt<=it,it) = NaN;
        end
        
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);

        subLPP = nan(numel(subsess),nmodel);
        for is = subsess'            
            itrt = double(unique(treatment(session==(is))));
            predictor = [];
            for ip=1:npredictor
                eval(['predictor = [ predictor , ' predictor_list{ip} ' ];']);
            end
            predictor = (predictor(session==(is),:));
            predictor = nanzscore(predictor);
            eval(['y = ' vartext '(session==is);']);
            
            
            for iM = 1:nmodel
                inputs = [ ones(numel(y),1) , predictor(:,model_matrix{iM}) ];
                [beta,L,posterior,out] = nanglm(inputs,y,'logit',1);
                 subLPP(is,iM) = L;
                 disp(L);
            end
        end

           
        % second-level stat
         dim=3;
         T = [T ; nominal(repmat({'clonidine','placebo','atomoxetine'}',maxSess,1)) ];
         BETA = [BETA ; reshape(ysub,[nregressor,ntrt*maxSess])' ];
         LPP = [LPP ; subLPP ];

        
    end
    
    
% model-selection
options.DisplayWin=0;
% between-session random-effect 
subject =  double(nominal(sessiondata.subject)) ;
[posterior,out] = VBA_groupBMC(LPP',options) ;
bms = table(mean(out.L,2),out.Ef,out.ep','VariableNames',{'expLPP','expFreq','excProba'});

% between-subject random-effect 
p = nan(3,1);
[h,p(1)] = VBA_groupBMC_btwGroups({LPP(subject==1,:)',LPP(subject==2,:)'} ,options);
[h,p(2)] = VBA_groupBMC_btwGroups({LPP(subject==2,:)',LPP(subject==3,:)'} ,options);
[h,p(3)] = VBA_groupBMC_btwGroups({LPP(subject==1,:)',LPP(subject==3,:)'} ,options);

% within-subject between-session random-effect 
subjectname = {'a','b','e'};
for isub=1:nsub
    subLPP = LPP(subject==isub,:)';
    options.DisplayWin=0;
   [posterior,out] = VBA_groupBMC(subLPP,options) ;
   if isub==1
       subbms = table(repmat(subjectname(isub),numel(out.Ef),1),mean(out.L,2),out.Ef,out.ep','VariableNames',{'Level','expLPP','expFreq','excProba'});
   else
       subbms = vertcat(subbms,table(repmat(subjectname(isub),numel(out.Ef),1),mean(out.L,2),out.Ef,out.ep','VariableNames',{'Level','expLPP','expFreq','excProba'}));
   end
end
fprintf('done\n');
    
%------------- END OF CODE --------------
