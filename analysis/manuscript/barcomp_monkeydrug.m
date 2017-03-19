function [stat,model] = barcomp_monkeydrug(dataset,yvar,yvartext,trt2comp,modelname)

    
    

    % variables
    vartext = yvar;
    ytext = yvartext;
    ylimits = [];
    
    % parameters
    subjectList = {'aliosha','bob','esmeralda'};
    nsub = numel(subjectList);
    trtList = {'clonidine','placebo','atomoxetine'};
    ntrt = numel(trtList);
    col = {[0.54 0.27 0.07],...
       [1 1 1]*0.5,...
       [1 0.5 0]};
    
        
    nbin = [ntrt,nsub];
    Y = nan(nbin);
    Z = nan(nbin);

    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(dataset.subject, sub);   
        selection = selectSubject ;

        % variables
        tab = dataset;
        session = tab.session_number(selection);
        treatment = tab.treatment(selection);
        participation_rate = tab.participation_rate(selection);
        accuracy_R = tab.accuracy_R(selection);
        accuracy_E = 1 - tab.accuracy_E(selection);
        choice_HRHE = tab.choice_HRHE(selection);
        accuracy_RE = tab.accuracy_RE(selection);
        session_reward = tab.session_reward(selection);
        varcoef_force = tab.varcoef_force(selection);
        varcoef_responsetime = tab.varcoef_responsetime(selection);
        peak_force = tab.peak_force(selection);
        responsetime = tab.responsetime(selection);
        for ivar = 1:size(dataset,2)
            varname = dataset.Properties.VariableNames{ivar};
            eval([ varname '= dataset.' varname ]);
        end

        eval(['vary = ' vartext ]);

        
        % subject-display
        if isub==1; fig = figure; set(fig,'Name','barcomp_subjects'); end
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
        fig = figure; set(fig,'Name','barcomp_group');
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
        stat=struct;
        switch modelname
            case 'FFX'
                model = fitglm(predictor,y,'linear','Intercept',false);
                disp(model.Coefficients);
                [stat.p,stat.F,stat.df1,stat.df2] = coefTest(model,contrast);
            case 'RFX'
                predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),tab.subject,y,...
                     'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
                model = fitglme(predictor2,[ vartext ' ~ -1 + clonidine + placebo +  atomoxetine + ( -1 + clonidine + placebo +  atomoxetine | subject)']);
                disp(model.Coefficients);
                [stat.p,stat.F,stat.df1,stat.df2] = coefTest(model,contrast);
        end
         [~,~,stat.coef] = fixedEffects(model);


                xx = x + (trt2comp(trt2comp~=2)-2)*0.25;
        [s] = sigstar( num2cell(xx),stat.p);
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);