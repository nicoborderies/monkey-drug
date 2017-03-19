function [dataset,stat] = multcomp_monkeydrug(designtab,dataset,yvar,paramNames,trt2comp,modelname)  

% variables
 vartext = 'y';
     
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
     
    
    % stat
    if ischar(yvar)
        Y = dataset.(yvar);  
    else
        Y = dataset{:,yvar};
    end
    nregressor = size(Y,2);


    T = dataset.treatment;
    predictor = [];
    predictor = [predictor ,(T=='clonidine')];
    predictor = [predictor ,(T=='placebo')];
    predictor = [predictor , (T=='atomoxetine')];
    subject = dataset.subject;
    stat=struct;
    stat.p = nan(1,nregressor);    stat.F = nan(1,nregressor);
    stat.pRFX = nan(3,nregressor); stat.F_RFX = nan(3,nregressor); 
    for i = 1:nregressor
        contrast = [0 0 0]; contrast(trt2comp) = 1; contrast(2) = -1;
        y = Y(:,i);

        switch modelname
            case 'FFX'
                % fixed-effect model
                model = fitglm(predictor,y,'linear','Intercept',false);
                [p(i),F(i)] = coefTest(model,contrast);
                [~,~,fixedstat] = fixedEffects(model);
                nffx = numel(fixedstat.Estimate);
                fixedstat.Name = repmat(paramNames(i),nffx,1) ;
            case 'RFX'
                if ~isempty(find(y~=0 & ~isnan(y)));
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
                else
                    fixedstat = mat2dataset(nan(3,8),'VarNames',{'Name','Estimate','SE','tStat','DF','pValue','Lower','Upper'});
                    fixedstat.Name = repmat(paramNames(i),3,1);
                    fixedstat.Estimate = repmat(0,3,1);
                    fixedstat.SE = repmat(0,3,1);
                    coef = vertcat(coef,fixedstat);

                end
        end
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


         % stores
         stat.coef = coef;
         stat.paramNames = paramNames;
         stat.p = p;
         stat.F = F;
         if exist('randomstat'); stat.rcoef = rcoef;end
end
     
end

         
        
