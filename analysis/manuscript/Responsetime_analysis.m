%% Responsetime_analysis
%   this script execute the following analysis:
%    - display between-treatments comparisons of rt
%    - display between-treatments comparisons of rt dependance
%      to experimental factors 
%    - glmm modeling of treatment effect on rt
%    - between-treatments comparisons of glmm coefficients
%    - model selection between different glmm models

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

%------------- BEGIN CODE --------------

%% average effect
    % variables
        nbin = [ ntrt , maxSess ];

    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);  
        selection = selectSubject & selectionChosen  ;

        % variables
        Y = nan(nbin);
        session = data.session(selection);
        treatment = data.treatment(selection);
        rt = data.responseTime(selection);
        
        % stats
        var = rt;
        [~,subtrt] = ismember(unique(treatment),trtList);
        [~,subsess] = ismember(unique(session),sessionList);
        ysub = tapply(var,{treatment,session},@nanmean);
        Y(subtrt,subsess) =  ysub;

        % averaging
         dim=2;
         y = nanmean(Y,dim);     
         z = sem(Y,dim);
%          z = nanstd(Y,[],dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','rt'); end
        subplot(1,nsub,isub);hold on ; clear h ;
        for it=1:3
            x = [1];
            xx = x + (it-2)*0.25;
            [ h(it)  ] = barplot( xx ,y(it),z(it), col{it} );
            h(it).BarWidth = 0.25;
%             
%             h = plotSpread(Y','distributionColors',col);
%             m = findobj('-property','MarkerFaceColor');
%             set(m,'MarkerSize',18);
            
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend([h(1) h(2) h(3)],trtList);
        end
        ax = gca; 
        ax.YLim = [0 1];
        ax.TickLength = [0 0];
        ax.XTick = []; 
        if isub==1 ; ylabel('response time (sec)'); end
        t = title(['monkey ' sub(1) ]); 
        
    end
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);

%% option value
%     factor = 'R';
%     factor = 'E';
    factor = 'T';

    % variables
        nbin = [ 7, ntrt , maxSess ];
%         nbin = [ 4, ntrt , maxSess ];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        selection = selectSubject & selectionChosen ;

           % variables
        Y = nan(nbin);
        session = data.session(selection);
        rt = data.responseTime(selection)+0.1;        
        log_rt = log(rt);        
        treatment = data.treatment(selection);
        
        r1 = round(data.ordinalRewardLeft(selection));
        r2 = round(data.ordinalRewardRight(selection));
        sumR = r1+r2;
        dr = r2-r1;
        e1 = round(data.ordinalEffortLeft(selection));
        e2 = round(data.ordinalEffortRight(selection));
        sumE = e1+e2;  
        de = e2-e1;
        rchoice = data.chosenCardinalReward(selection);
        echoice = data.chosenCardinalEffort(selection);
        nt = data.trialNumber(selection);
        nnt = data.normalizedTrialNumber(selection);
        nnt = quantileranks(nnt,nbin(1));  

        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        [~,subsess] = ismember(unique(session),sessionList);
        
            % mean participation
            if factor=='R'; vary = sumR;
            elseif factor=='E';  vary = sumE; 
            elseif factor=='T';  vary = nnt; 
            end
            ysub = tapply(rt,{vary,treatment,session},@nanmean);
            Y(:,subtrt,subsess) =  ysub;
            

        % second-level stat
         dim=3;
         y = nanmean(Y,dim);
         z = sem(Y,dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','rt_value'); end
        subplot(1,nsub,isub);hold on ; clear h
        for it=1:3
%         for it=1:2
            x = [1:nbin(1)];
            xx = x ;
            [~,~,h(it)] = errorscat( xx ,y(:,it), z(:,it),col{it});
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend([h(1) h(2) h(3)],trtList);
        end
        ax = gca;
%         ax.YLim = [0 1];
        ax.TickLength = [0 0];
        ax.XTick = []; 
        if factor=='R'; xlabel(texlabel('R_left + R_right'));
        elseif factor=='E'; xlabel(texlabel('E_left + E_right')); 
        elseif factor=='T'; xlabel(texlabel('session progression (%)'));  end
        if isub==1 ; ylabel('response time (sec)'); end
        t = title(['monkey ' sub(1) ]); 
        
    end

    % reformat
    setFigProper('FontSize',20,'LineWidth',2);

%% GLMM models

 % parameters
     trt2comp = [2];
     vartext = 'responsetime';
     formula = ' responsetime ~ 1 + sumR + sumE + offertime  + dimensionality + congruency';
     varnames = {'sumR';'sumE';'offertime';'dimensionality';'congruency';'responsetime'};
     effectName = {texlabel('beta_0'),'sumR','sumE','offertime','dimensionality','congruency',};
     fname= 'identity';  distname = 'normal';
     dataselect = 'selectSubject & selectionChosen';
     predictor_list = {'sum_r';'sum_e';'offertime';'dimensionality';'congruency'};
     npredictor = numel(predictor_list)+1;
     nmodel=2;
     nbin = [ npredictor , ntrt , maxSess ];
     nbin2 = [ npredictor ,ntrt,nsub ];
     
    % prepare variables
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
        
        p_choice = sig(data.participation_prediction(selection));
        p_choice2 = sig(data.choice_prediction(selection));
        p_choice1 = 1-p_choice2;
        uncertainty = - ( (1-p_choice).*log(1-p_choice) + (p_choice.*p_choice2).*log(p_choice.*p_choice2) + (p_choice.*p_choice1).*log(p_choice.*p_choice1) );

        responsetime = data.responseTime(selection)+eps;
        % previousParticipation
        tmax=1;
        p_t = nan(numel(participation),tmax);
        for it=1:tmax
            p_t(1+it:end,it) = participation(1:end-it);
            p_t(nt<=it,it) = NaN;
        end
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);

        for is = subsess'
            
            predictor = [];
            for ip=1:(npredictor-1)
                eval(['predictor = [ predictor , ' predictor_list{ip} ' ];']);
            end
            predictor = predictor(session==(is),:);
            predictor = nanzscore(predictor);
            eval(['y = ' vartext '(session==is);']);
            stat = fitglm(predictor,y,...
                           formula,'VarNames',varnames,'Distribution',distname,'link',fname);
            coef = stat.Coefficients;
            disp(coef);
            
            itrt = double(unique(treatment(session==is)));
            ysub(:,itrt,is) = coef.Estimate;
        end
           
        % second-level stat
         dim=3;
         T = [T ; nominal(repmat({'clonidine','placebo','atomoxetine'}',maxSess,1)) ];
         BETA = [BETA ; reshape(ysub,[npredictor,ntrt*maxSess])' ];
         y = nanmean(ysub,dim); Y(:,subtrt,isub) = y;
         z = sem(ysub,dim); Z(:,subtrt,isub) = z;
         
        % display
        if isub==1; fig = figure; set(fig,'Name','glm_choice'); end
        subplot(nsub,1,isub);hold on ; clear h ; h = gobjects(ntrt,1);
        for it=trt2comp
            x = [1:npredictor];
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
        x = [1:npredictor];
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
    p = nan(npredictor,1);
    for i = 1:npredictor
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
 % parameters
     trt2comp = [2];
     vartext = 'responsetime';
     formula = ' responsetime ~ 1 + sumR + sumE + offertime  + dimensionality + congruency';
     varnames = {'sumR';'sumE';'offertime';'dimensionality';'congruency';'responsetime'};
     effectName = {texlabel('beta_0'),'sumR','sumE','offertime','dimensionality','congruency',};
     fname= 'identity';  distname = 'normal';
     dataselect = 'selectSubject & selectionChosen';
     predictor_list = {'sum_r';'sum_e';'dr';'de';'offertime';
                        'dimensionality';'congruency';'uncertainty';'rchoice';'echoice'};
     npredictor = numel(predictor_list)+1;
     model_matrix = {[1 2 5 6 7] ; 
                     [3 4 5 6 7];
                     [5 6 7 8];
                     [5 6 7 9 10]};
     nmodel=size(model_matrix,1);
     nbin = [ npredictor , ntrt , maxSess ];
     nbin2 = [ npredictor ,ntrt,nsub ];
     
    % prepare variables
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
        dr = abs(sdR);
        sdE = (round(data.deltaOrdinalEffort(selection)));
        choiceE = (sign(choice)==sign(sdE));
        de = abs(sdE);
        r1 = round(data.ordinalRewardLeft(selection));
        r2 = round(data.ordinalRewardRight(selection));
        e1 = round(data.ordinalEffortLeft(selection));
        e2 = round(data.ordinalEffortRight(selection));
        sum_r = r1+r2;
        sum_e = e1+e2;
        rchoice = data.chosenCardinalReward(selection);
        echoice = data.chosenCardinalEffort(selection);
        congruency = double(sign(sdR)==-sign(sdE)) ;
        dimensionality = (sdR~=0) + (sdE~=0) ;
        congruency(dimensionality~=2)=1;
        offertime = data.offerTime(selection);
        offertimebin = quantileranks(offertime,nbin(1)); 
        
        p_choice = sig(data.participation_prediction(selection));
        p_choice2 = sig(data.choice_prediction(selection));
        p_choice1 = 1-p_choice2;
        uncertainty = - ( (1-p_choice).*log(1-p_choice) + (p_choice.*p_choice2).*log(p_choice.*p_choice2) + (p_choice.*p_choice1).*log(p_choice.*p_choice1) );
        
        responsetime = data.responseTime(selection)+eps;
        % previousParticipation
        tmax=1;
        p_t = nan(numel(participation),tmax);
        for it=1:tmax
            p_t(1+it:end,it) = participation(1:end-it);
            p_t(nt<=it,it) = NaN;
        end
        
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);

        subLPP = nan(numel(subsess),nmodel);
        for is = subsess'            
            itrt = double(unique(treatment(session==(is))));
            predictor = [];
            for ip=1:(npredictor-1)
                eval(['predictor = [ predictor , ' predictor_list{ip} ' ];']);
            end
            predictor = (predictor(session==(is),:));
            predictor = nanzscore(predictor);
            eval(['y = ' vartext '(session==is);']);
            
            
            for iM = 1:nmodel
                inputs = [ ones(numel(y),1) , predictor(:,model_matrix{iM}) ];
                [beta,L,posterior,out] = nanglm(inputs,y,'logit',0);
                 subLPP(is,iM) = L;
                 disp(L);
            end
        end

           
        % second-level stat
         dim=3;
         T = [T ; nominal(repmat({'clonidine','placebo','atomoxetine'}',maxSess,1)) ];
         BETA = [BETA ; reshape(ysub,[npredictor,ntrt*maxSess])' ];
         LPP = [LPP ; subLPP ];

        
    end
    
    
% model-selection
disp=0;
options.DisplayWin=disp;
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
    options.DisplayWin=disp;
   [posterior,out] = VBA_groupBMC(subLPP,options) ;
   if isub==1
       subbms = table(repmat(subjectname(isub),numel(out.Ef),1),mean(out.L,2),out.Ef,out.ep','VariableNames',{'Level','expLPP','expFreq','excProba'});
   else
       subbms = vertcat(subbms,table(repmat(subjectname(isub),numel(out.Ef),1),mean(out.L,2),out.Ef,out.ep','VariableNames',{'Level','expLPP','expFreq','excProba'}));
   end
end
fprintf('done\n');
    
%------------- END OF CODE --------------

