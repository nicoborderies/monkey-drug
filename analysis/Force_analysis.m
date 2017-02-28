%% Force_analysis
%   this script execute the following analysis:
%    - display between-treatments comparisons of force
%    - display between-treatments comparisons of force dependance
%      to experimental factors 
%    - glmm modeling of treatment effect on force
%    - between-treatments comparisons of glmm coefficients

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
        nsess = max(design.sessionNumber);
        nbin = [ ntrt , nsess ];

    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);  
        selection = selectSubject & selectionChosen  ;

        % variables
        Y = nan(nbin);
%         session = data.session(selection);
        session = data.sessionNumber(selection);

        treatment = data.treatment(selection);
        force = data.perf(selection);
%         force = data.time2peak(selection);
%         force = data.cumulativePerf(selection);
        
        % stats
        var = force;
        [~,subtrt] = ismember(unique(treatment),trtList);
%         [~,subsess] = ismember(unique(session),sessionList);
        subsess = unique(session);

        ysub = tapply(var,{treatment,session},@nanmean);
        Y(subtrt,subsess) =  ysub;

        % averaging
         dim=2;
         y = nanmean(Y,dim);     
         z = sem(Y,dim);
%          z = nanstd(Y,[],dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','force'); end
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
            
%             xx = double(unique(sessionList));
%             h(it) = scatter(xx,Y(it,:),40,col{it},'filled');
            
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend([h(1) h(2) h(3)],trtList);
        end
        ax = gca; 
%         ax.YLim = [0 1];
        ax.TickLength = [0 0];
        ax.XTick = []; 
        if isub==1 ; ylabel('force peak (%fmax)'); end
%         if isub==1 ; ylabel('velocity peak (AU)'); end
%         if isub==1 ; ylabel('time to peak (sec)'); end

        t = title(['monkey ' sub(1) ]); 
        
    end
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);



%% option value
% effect of conditions on force
% ******************************

    % variables
%         nbin = [ 4, ntrt , maxSess ];
        nbin = [ 7, ntrt , maxSess ];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        selection = selectSubject & selectionChosen ;
         
        % variables
        Y = nan(nbin);
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);   
        t2p = data.time2peak(selection);   
        t2p = quantileranks(t2p,nbin(1)); 
        yank = data.yankPeak(selection);   
        yank = quantileranks(yank,nbin(1)); 
        yank(yank==0) = 1;
        rchoice = data.chosenCardinalReward(selection);
        echoice = data.chosenCardinalEffort(selection);
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        [~,subsess] = ismember(unique(session),sessionList);
        
        % mean participation
        varcoef  = @(x) nanvar(x)/nanmean(x);
        statistic = @nanmean;
%         statistic = varcoef;
        varx = yank;
        xsub = unique(varx);
        ysub = tapply(force,{varx,treatment,session},statistic);
        Y(:,subtrt,subsess) =  ysub;
            

        % second-level stat
         dim=3;
         y = nanmean(Y,dim);
         z = sem(Y,dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','force'); end
        subplot(1,nsub,isub);hold on ; clear h
        x = [1:nbin(1)];
%         plot(xsub,xsub,'--k');
        for it=1:3
            xx = xsub(x) ;
            [~,~,h(it)] = errorscat( xx ,y(x,it), z(x,it),col{it});
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend([h(1) h(2) h(3)],trtList);
        end
        ax = gca;
%         ax.XLim = [0 1];
%         ax.YLim = [0 1];
        ax.TickLength = [0 0];
        ax.XTick = []; 
%         xlabel(texlabel('instructed force (%fmax)'));
%         xlabel(texlabel('time-to-peak (sec)'));
        xlabel(texlabel('velocity peak (%fmax/sec)'));
        ylabel(texlabel('exerted force (%fmax)'));
        t = title(['monkey ' sub(1) ]); 
        
    end

    % reformat
    setFigProper('FontSize',20,'LineWidth',2);    

% effect of fatigue on force
% **************************

    % variables
        nbin = [ 7, ntrt , maxSess ];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(data.subject, sub);   
        selection = selectSubject & selectionChosen ;
         
        % variables
        Y = nan(nbin);
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);   
        rchoice = data.chosenCardinalReward(selection);
        echoice = data.chosenCardinalEffort(selection);
        nnt = data.normalizedTrialNumber(selection);
        nnt = quantileranks(nnt,nbin(1)); 
        
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        [~,subsess] = ismember(unique(session),sessionList);
        
        % mean participation
        varcoef  = @(x) nanvar(x)/nanmean(x);
        statistic = @nanmean;
%         statistic = varcoef;
        xsub = unique(echoice);
        ysub = tapply(force,{nnt,treatment,session},statistic);
        Y(:,subtrt,subsess) =  ysub;
            

        % second-level stat
         dim=3;
         y = nanmean(Y,dim);
         z = sem(Y,dim);
        
        % display
        if isub==1; fig = figure; set(fig,'Name','force_fatigue'); end
        subplot(1,nsub,isub);hold on ; clear h
        for it=1:3
            x = [1:nbin(1)];
            xx = x ;
            [~,~,h(it)] = errorscat( xx ,y(x,it), z(x,it),col{it});
        end
        
        % legending
        if isa(h,'matlab.graphics.Graphics') && isub==nsub
            legend([h(1) h(2) h(3)],trtList);
        end
        ax = gca;
%         ax.XLim = [0 1];
%         ax.YLim = [0 1];
        ax.TickLength = [0 0];
        ax.XTick = []; 
        xlabel(texlabel('session progression (%)')); 
        ylabel(texlabel('exerted force (%fmax)'));
        t = title(['monkey ' sub(1) ]); 
        
    end

    % reformat
    setFigProper('FontSize',20,'LineWidth',2);


%% GLMM models

 % parameters
     trt2comp = [2];
     vartext = 'force';
     formula = ' force ~ 1 + side + side:echoice + rchoice*echoice + ntrial  ';
     varnames = {'rchoice';'echoice';'side';'ntrial';'Rtot';'Etot';'force'};
     effectName = {'k_0','rchoice','echoice','side','ntrial','rchoice*echoice ','echoice*side'};
     fname= 'identity';  distname = 'normal';
     dataselect = 'selectSubject & selectionChosen';
     predictor_list = {'rchoice';'echoice';'side*2-1';'nt';'rr';'ee'};
     npredictor = numel(predictor_list)+1;
     nmodel=1;
     nfactor = numel(varnames)-1;
     nbin = [ npredictor , ntrt , maxSess ];
     nbin2 = [ npredictor ,ntrt,nsub ];
     
    % prepare variables
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
        force = data.perf(selection);   
        treatment = data.treatment(selection);
        
        rchoice = data.chosenCardinalReward(selection);
        echoice = data.chosenCardinalEffort(selection);
        side = data.sideChoice(selection)-1;   
        nt = data.trialNumber(selection);
        nnt = data.normalizedTrialNumber(selection);
        ee =  data.cumulativeEffort(selection);
        rr =  data.cumulativeReward(selection);
        
        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);

        for is = subsess'
            
            predictor = [ rchoice, echoice , side*2-1 , nt , rr , ee ];
            predictor = predictor(session==(is),:);
            predictor = nanzscore(predictor);
            y = force(session==(is));
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

