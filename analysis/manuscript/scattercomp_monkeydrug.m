function [stat] = scattercomp_monkeydrug(designtab,dataset,xvar,xvarname,yvar,yvarname,modelname,xlimits,ylimits,trt2comp)  
 

% variables
varytext = yvar;
varxtext = xvar;
xtext = xvarname;
ytext = yvarname;
    

% parameters
 %%% display option
 col = {[0.54 0.27 0.07],...
    [1 1 1]*0.5,...
    [1 0.5 0]};
%%% design information
 subjectList = {'aliosha','bob','esmeralda'};
 nsub = numel(subjectList);
 trtList = {'clonidine','placebo','atomoxetine'};
 ntrt = numel(trtList);
 sessionList = unique(dataset.session);
 maxSess = max(designtab.sessionNumber);


xticks = [];
xlabels = [];
BETA = []; RHO = []; T = [];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(dataset.subject, sub);   
        selection = selectSubject ;

        % variables
        session = dataset.session_number(selection);
        treatment = dataset.treatment(selection);
        participation_rate = dataset.participation_rate(selection);
        accuracy_R = dataset.accuracy_R(selection);
        accuracy_E = 1 - dataset.accuracy_E(selection);
        choice_HRHE = dataset.choice_HRHE(selection);
        accuracy_RE = dataset.accuracy_RE(selection);
        varcoef_force = dataset.varcoef_force(selection);
        varcoef_responsetime = dataset.varcoef_responsetime(selection);
        peak_force = dataset.peak_force(selection);
        responsetime = dataset.responsetime(selection);
        success_force = dataset.success_force(selection);
        choice_E = dataset.choice_E(selection);
        for ivar = 1:size(dataset,2)
            varname = dataset.Properties.VariableNames{ivar};
            eval([ varname '= dataset.' varname ]);
        end
        
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
        if ~isempty(xlimits); ax.XLim = xlimits; end
        if ~isempty(ylimits); ax.YLim = ylimits; end
        ax.TickLength = [0 0];
        xlabel(xtext); 
        ylabel(ytext); 
        
        % stat
        predictor = [];
        predictor = [predictor ,(dataset.treatment=='clonidine')];
        predictor = [predictor ,(dataset.treatment=='placebo')];
        predictor = [predictor , (dataset.treatment=='atomoxetine')];
        predictor = [predictor , zscore(dataset.(varxtext))];
        y = zscore(dataset.(varytext));
       
        stat=struct;
        switch modelname
            case 'betweenConditionConfound'
                contrast = [0 0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
                predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),predictor(:,4),dataset.subject,y,...
                     'VariableNames',{'clonidine','placebo','atomoxetine','confound','subject',varytext});
                model = fitglme(predictor2,[ varytext ' ~ -1 + confound + clonidine + placebo +  atomoxetine + ( -1 + confound + clonidine + placebo +  atomoxetine | subject)']);
            case 'withinConditionConfound'
                contrast = [0 1 0 0 ; 0 0 1 0 ; 0 0 0 1 ];
                predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),predictor(:,4),dataset.subject,y,...
                     'VariableNames',{'clonidine','placebo','atomoxetine','confound','subject',varytext});
                model = fitglme(predictor2,[ varytext ' ~ 1 + confound:clonidine + confound:placebo +  confound:atomoxetine + ( 1 + confound:clonidine + confound:placebo +  confound:atomoxetine | subject)']);
        end
        disp(model.Coefficients);
        [stat.p,stat.F,stat.df1,stat.df2] = coefTest(model,contrast);
        [~,~,stat.coef] = fixedEffects(model);
        [~,~,stat.rcoef] = randomEffects(model);
        xx = x + (trt2comp(trt2comp~=2)-2)*0.25;
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
end