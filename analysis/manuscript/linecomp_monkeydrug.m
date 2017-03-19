function [stat] = scattercomp_monkeydrug(designtab,dataset,xvar,xvarname,yvar,yvarname,condvar,modelname,nqtle,xxplot,xlimits,ylimits,dataselection,trt2comp)  
 

% variables
vartext = yvar;
factor = xvar;
binfactor = condvar;
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


% xlimits = [0.5 4.5];
% ylimits = [0 1];
% xxplot = 1;
nbin = [ nqtle , ntrt , maxSess ];
nbin2 = [ nqtle ,ntrt,nsub ];
xbin = [1:nbin(1)];
xticks = [];
xlabels = [];
fstat = @nanmean;

X = nan(nbin2);
Y = nan(nbin2);
Z = nan(nbin2);
BETA = []; RHO = [];


    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(dataset.subject, sub);   
        eval(['selection = ' dataselection ';']);
        
        % variables
        xsub = nan(nbin);
        ysub = nan(nbin);
        
        session = dataset.sessionNumber(selection);
        treatment = dataset.treatment(selection);
        
        participation_rate = selectionParticipate(selection);
%         participation_prediction = dataset.participation_prediction(selection);
%         participation_predictionbin = quantileranks(participation_prediction,nbin(1)); 
%         participation_predictionbin(participation_predictionbin==0)=NaN;
        
        choice = (dataset.sideChoice(selection)-1)*2-1;
        sdR = (round(dataset.deltaOrdinalReward(selection)));
        choiceR = (sign(choice)==sign(sdR));
        dR = abs(sdR);
        sdE = (round(dataset.deltaOrdinalEffort(selection)));
        choiceE = (sign(choice)==sign(sdE));
        dE = abs(sdE);
        r1 = round(dataset.ordinalRewardLeft(selection));
        r2 = round(dataset.ordinalRewardRight(selection));
        e1 = round(dataset.ordinalEffortLeft(selection));
        e2 = round(dataset.ordinalEffortRight(selection));
        sumR = r1+r2;
        sumE = e1+e2;
        congruency = double(sign(sdR)==-sign(sdE)) ;
        dimensionality = (sdR~=0) + (sdE~=0) ;
        congruency(dimensionality~=2)=NaN;
        
        offertime = dataset.intertrialTime(selection);
        offertimebin = quantileranks(offertime,nbin(1)); 

        
        accuracy = (choiceR | ~choiceE) ;
        responsetime = dataset.responseTime(selection);

        responsetimebin = quantileranks(responsetime,nbin(1)); 
        responsetimebin(responsetimebin==0)=NaN;
        echoice = dataset.chosenCardinalEffort(selection);
        rchoice = dataset.chosenCardinalReward(selection);
        peak_force = dataset.perf(selection);   
        peak_force_bin = quantileranks(peak_force,nbin(1)); 
        peak_velocity = dataset.yankPeak(selection);   
        time2peak = dataset.time2peak(selection);   
        nnt = dataset.normalizedTrialNumber(selection);
        nnt = quantileranks(nnt,nbin(1)); 
        nt = dataset.trialNumber(selection);
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
        
        for is = subsess'
            x = varx(session==is);
            y = vary(session==is);
            beta = glmfit(x,y,'normal');
            rho = corr(x,y,'type','Pearson','rows','pairwise');
            BETA = [BETA;beta(2)];
            RHO = [RHO;rho];
        end
            

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
            [~,s,h(it)] = errorscat( x(xx,it) ,y(xx,it), z(xx,it),col{it});
            s.SizeData = 20;
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
        [~,s,h(it)] = errorscat( x(xx) ,y(xx), z(xx),col{it});
        s.SizeData = 20;
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
    
     % stat
    predictor = [];
    predictor = [predictor ,(designtab.treatment=='clonidine')];
    predictor = [predictor ,(designtab.treatment=='placebo')];
    predictor = [predictor , (designtab.treatment=='atomoxetine')];
    switch modelname
        case 'glm'
            y = BETA;
        case 'corr'
            y = RHO;
    end

    contrast = [0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
    stat=struct;
    predictor2 = table(predictor(:,1),predictor(:,2),predictor(:,3),designtab.subject,y,...
         'VariableNames',{'clonidine','placebo','atomoxetine','subject',vartext});
    model = fitglme(predictor2,[ vartext ' ~ -1 + clonidine + placebo +  atomoxetine + ( -1 + clonidine + placebo +  atomoxetine | subject)']);
    disp(model.Coefficients);
    [stat.p,stat.F,stat.df1,stat.df2] = coefTest(model,contrast);
     [~,~,stat.coef] = fixedEffects(model);
    xx = x + (trt2comp(trt2comp~=2)-2)*0.25;
    [s] = sigstar(num2cell(median(xx)),stat.p);
    
    
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);  
    
end