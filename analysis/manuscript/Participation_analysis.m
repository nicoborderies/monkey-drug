%% Participation_analysis
%   this script execute the following analysis:
%    - display between-treatments comparisons of participation
%    - display between-treatments comparisons of participation dependance
%      to experimental factors 
%    - glmm modeling of treatment effect on participation
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

%------------- BEGIN CODE --------------

%% average effect
     
    % variables


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(sessiondata.subject, sub);   
        selection = selectSubject ;

        % variables
        tab = sessiondata;
        session = tab.session_number(selection);
        treatment = tab.treatment(selection);
        participation = tab.participation_rate(selection);

        vary = participation;
        
        % display
        if isub==1; fig = figure; set(fig,'Name','participation'); end
        subplot(1,nsub,isub);hold on ; clear h
        for it=1:3
            % inputs
            x = [1];
            xx = x + (it-2)*0.25;
            indtrt = ismember(treatment,trtList(it));
            y = nanmean(vary(indtrt));
            z = sem(vary(indtrt));
            % plot
            [ h(it)  ] = barplot( xx ,y,z, col{it} );
            h(it).BarWidth = 0.25;
%             % scatterplot 1
%             h = plotSpread(Y','distributionColors',col);
%             m = findobj('-property','MarkerFaceColor');
%             set(m,'MarkerSize',18);
%             % scatterplot 2
%             xx = double(unique(sessionList));
%             h(it) = scatter(xx,Y(it,:),40,col{it},'filled');
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend([h(1) h(2) h(3)],trtList);
        end
        ax = gca;
%         ax.YLim = [0.3 1];
        ax.TickLength = [0 0];
        ax.XTick = []; 
        if isub==1 ; ylabel('participation (%)'); end
%         if isub==1 ; ylabel('participation (n° trials)'); end
%         if isub==1 ; ylabel('correct (n° trials)'); end
%         if isub==1 ; ylabel('session length (n° trials)'); end
%         if isub==1 ; ylabel('participation duration (n° trials)'); end
%         if isub==1 ; ylabel('omission duration (n° trials)'); end
        t = title(['monkey ' sub(1) ]); 
        
    end

    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
    
%% effect of option value
    factor = 'R';
%     factor = 'E';
%     factor = 'T';

    % variables
        nbin = [ 7, ntrt , maxSess ];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        selection = selectSubject ;
%         selection = selectSubject & ~selectionRepetition;

        % variables
        X = nan(nbin);
        Y = nan(nbin);
        session = data.session(selection);
        treatment = data.treatment(selection);
        
        participation = selectionParticipate(selection);
        sumE = round(data.ordinalEffortLeft(selection) + data.ordinalEffortRight(selection));
        sumR = round(data.ordinalRewardLeft(selection) + data.ordinalRewardRight(selection));
        nt = data.trialNumber(selection);
        nnt = data.normalizedTrialNumber(selection);
        nnt = quantileranks(nnt,nbin(1));  
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        [~,subsess] = ismember(unique(session),sessionList);
        
            % mean participation
            if factor=='R'; varx = sumR; varx2 = varx;
            elseif factor=='E';  varx = sumE; varx2 = varx;
            elseif factor=='T';  varx = nnt; varx2 = nt;
            end
            xsub = tapply(varx2,{varx,treatment,session},@nanmean);
            ysub = tapply(participation,{varx,treatment,session},@nanmean);
            X(:,subtrt,subsess) =  xsub;
            Y(:,subtrt,subsess) =  ysub;
            

        % second-level stat
         dim=3;
         x = nanmean(X,dim);
         y = nanmean(Y,dim);
         z = sem(Y,dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','participation'); end
        subplot(1,nsub,isub);hold on ; clear h
%         for it=1:3
        for it=1:2
            [~,~,h(it)] = errorscat( x(:,it)  ,y(:,it), z(:,it),col{it});
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend([h(1) h(2) h(3)],trtList);
        end
        ax = gca;
        ax.YLim = [0 1];
        ax.TickLength = [0 0];
        ax.XTick = []; 
        if factor=='R'; xlabel(texlabel('R_left + R_right'));
        elseif factor=='E'; xlabel(texlabel('E_left + E_right')); 
        elseif factor=='T'; xlabel(texlabel('session progression (ntrial)'));  end
        if isub==1 ; ylabel('participation (%)'); end
        t = title(['monkey ' sub(1) ]); 
        
    end

    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
%% GLMM models
% variables
     trt2comp = [1 2];
     vartext = 'participation';
     formula = ' participation ~ 1 + sumR + sumE  + ntrial + part_t_1 ';
     varnames = {'sumR';'R2';'sumE';'E2';'ntrial';'part_t_1';'Rtot';'Etot';'work_t';'pause_t';'cum_part';'participation'};
     effectName = {texlabel('beta_0'),'sumR','sumE','ntrial','part_t_1'};
     nregressor = 5;
     fname= 'logit';  distname = 'binomial';
     dataselect = 'selectSubject';
     predictor_select = '[ sum_r, r2 , sum_e , e2 , nt , p_t*2-1 , rr , ee, work_t, pause_t,cum_part ];';


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
            predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),subject,y,...
                 'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
            model2 = fitglme(predictor2,[ vartext ' ~ clonidine + atomoxetine + (clonidine|subject) + (atomoxetine|subject)']);
            disp(model2.Coefficients);
%             [p(i),F] = coefTest(model2,contrast);
            [p(i),F] = coefTest(model2,[0 1 0]);

            % average effect
%             predictor2 = table(subject,y,...
%                  'VariableNames',{'subject',vartext});
%             model2 = fitglme(predictor2,[ vartext ' ~ 1 + (1|subject)']);
%             [~,~,randomstat] = randomEffects(model2);
%             randomstat.Name = repmat(effectName(i),3,1) ;
%             disp(model2.Coefficients);
%             if i==1
%                 coef = model2.Coefficients;
%                 rcoef = randomstat;
%             else
%                 coef = vertcat(coef,model2.Coefficients);
%                 rcoef = vertcat(rcoef,randomstat);
%             end
%             p(i) = model2.Coefficients.pValue;
        
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
     vartext = 'participation';
     formula = ' participation ~ 1 + sumR + sumE  + ntrial + part_t_1 ';
     varnames = {'sumR';'maxR';'R2';'sumE';'minE';'E2';'ntrial';'part_t_1';
                'Rtot';'Etot';'work_t';'pause_t';'cum_part';'participation'};
     effectName = {texlabel('beta_0'),'sumR','sumE','ntrial','part_t_1'};
     nregressor = 5;
     fname= 'logit';  distname = 'binomial';
     dataselect = 'selectSubject';
%      predictor_select = '[ sum_r, max_r, r2 , sum_e , min_e, e2 , nt , p_t*2-1 , rr , ee, work_t, pause_t,cum_part, p_success , sum_psuccess ];';
     predictor_list = {'sum_r','max_r','r2','sum_e','min_e',...
                        'e2','nt','p_t*2-1','rr','ee',...
                        'work_t','pause_t','cum_part','p_success','sum_psuccess'};
     npredictor = numel(predictor_list);
     model_matrix = [1 4 7 8  ; 
                     2 5 7 8 ;
                     1 14 7 8 ;
                     1 15 7 8];
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
        p_success = data.successRate(selection);
        sum_psuccess = data.successRate_left(selection) + data.successRate_right(selection);
        congruency = double(sign(sdR)==-sign(sdE)) ;
        dimensionality = (sdR~=0) + (sdE~=0) ;
        congruency(dimensionality~=2)=1;
        offertime = data.offerTime(selection);
        offertimebin = quantileranks(offertime,nbin(1)); 
        
        responsetime = data.responseTime(selection)+eps;
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

        subLPP = nan(numel(subsess),nmodel);
        for is = subsess'
            cum_part(session==(is)) = cumsum(participation(session==(is)));
            
            itrt = double(unique(treatment(session==(is))));
            predictor = [];
            for ip=1:npredictor
                eval(['predictor = [ predictor , ' predictor_list{ip} ' ];']);
            end
            predictor = (predictor(session==(is),:));
            predictor = nanzscore(predictor);
            eval(['y = ' vartext '(session==is);']);
            
            
            for iM = 1:nmodel
                inputs = [ ones(numel(y),1) , predictor(:,model_matrix(iM,:)) ];
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