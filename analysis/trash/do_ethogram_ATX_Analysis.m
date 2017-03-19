%% do_ethogram_ATX_Analysis

%     clear all;
%     close all;
%     clear h ;
%     clear b;

%% renaming
    data.trt = nominal(data.trt);
    data.trt(data.trt=='free')='placebo';
    data.trt(data.trt=='atomoxetine')='atx';
    data.trt = removecats(data.trt);
    
    data.trt = nominal(data.trt);
    data.trt(data.trt=='placebo')='1_placebo';
    data.trt(data.trt=='atx')='2_atx';
    data.trt = removecats(data.trt);
%%
    subjectList = unique(data.subject);
    sessionList = unique(data.session);
    behaviorList = getlevels(data.behavior);
    trtList = unique(data.trt);
    trtSession=cell(2,1);
    for i=1:2
    trtSession{i} = find(ismember(sessionList,unique(data.session(data.trt == trtList(i)))));
    end
    col = {[1 1 1]*0.5,[1 0.5 0]};

%%

fig = figure; set(fig,'Name','ethogram');
    for iSub = 1:numel(subjectList)
        sub = subjectList(iSub);
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject  ;
        
        session = data.session(selection);
        treatment = data.trt(selection);
        behavior = double(data.behavior(selection));
        time = data.time(selection);
        
        freq = nan(numel(sessionList),numel(behaviorList));
        cumtime = nan(numel(sessionList),numel(behaviorList));

        for is = 1:numel(sessionList)
            ind = tapply(behavior,{behavior},@unique);

            freq(is,:) = histcounts(behavior(session==sessionList(is)),[1:13],'Normalization','Probability');
            cumtime(is,ind) = tapply(time,{behavior},@nansum);
        end
        
        for it=1:2
            
%             subplot(2,numel(subjectList),(it-1)*2+iSub);hold on

            
            x = [1:12];
            y = nanmean(freq(trtSession{it},:),1);
            z = sem(freq(trtSession{it},:),1);
            
           [h,~,l] = errorscat(x ,y, z,col{it});
%             l.LineStyle = 'none';
            
            ax = gca;    
            ax.XTick = x;
            ax.XTickLabel = cellstr(behaviorList);
            ax.XTickLabelRotation = 45;
            ax.YLim = [0 0.5];
            xlabel('behavioral categories'); 
            ylabel('frequency'); 
            yy =ylim;
            t = title(['monkey ' char(sub(1)) ]); 
            
        end
       

    end
    
%%

fig = figure; set(fig,'Name','ethogram2');
    for iSub = 1:numel(subjectList)
        sub = subjectList(iSub);
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
        selection = selectSubject  ;
        
        session = data.session(selection);
        treatment = data.trt(selection);
        behavior = double(data.behavior(selection));
        time = data.time(selection);
        
        freq = nan(numel(sessionList),numel(behaviorList));
        cumtime = zeros(numel(sessionList),numel(behaviorList));

        for is = 1:numel(sessionList)
            ind = tapply(behavior(session==sessionList(is)),{behavior(session==sessionList(is))},@unique);

            freq(is,:) = histcounts(behavior(session==sessionList(is)),[1:13],'Normalization','Probability');
            cumtime(is,ind) = tapply(time(session==sessionList(is)),{behavior(session==sessionList(is))},@nansum);
        end
        
        for it=1:2
            

            
            x = [1:12];
            y = nanmean(cumtime(trtSession{it},:),1);
            z = sem(cumtime(trtSession{it},:),1);
            
           [h,~,l] = errorscat(x ,y, z,col{it});
%             l.LineStyle = 'none';
            
            ax = gca;    
            ax.XTick = x;
            ax.XTickLabel = cellstr(behaviorList);
            ax.XTickLabelRotation = 45;
%             ax.YLim = [0 0.5];
            xlabel('behavioral categories'); 
            ylabel('mean time (sec)'); 
            yy =ylim;
            t = title(['monkey ' char(sub(1)) ]); 
            
        end
       

    end